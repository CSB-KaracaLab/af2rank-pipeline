from pathlib import Path

import numpy as np
import pytest

from af2rank_pipeline import cli
from af2rank_pipeline.actifptm_runner import (
    LINKED_CHAIN,
    drop_linker_outputs,
    linked_asym_id,
    read_pdb_residues,
    unlink_pdb,
    write_linked_pdb,
)
from af2rank_pipeline.af2rank_runner import AF2RankConfig, model_subdir

REPO = Path(__file__).resolve().parents[1]


def residue_line(serial, atom, resname, chain, resseq, x):
    return (f"ATOM  {serial:5d}  {atom:<3s} {resname} {chain}{resseq:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           {atom[0]}\n")


def two_chain_pdb(path):
    lines, serial = [], 1
    # Chain A has a gap (residue 4 missing); chain B starts at 2.
    for chain, numbers, start_x in (("A", (3, 5, 6), 0.0), ("B", (2, 3), 30.0)):
        for i, resseq in enumerate(numbers):
            for atom in ("N", "CA", "C", "O"):
                lines.append(residue_line(serial, atom, "ALA", chain, resseq, start_x + 3.8 * i))
                serial += 1
        lines.append("TER\n")
    path.write_text("".join(lines) + "END\n")
    return path


def test_linked_pdb_joins_chains_with_glycines_and_keeps_gaps(tmp_path):
    source = two_chain_pdb(tmp_path / "model.pdb")
    position_map = write_linked_pdb(source, ["A", "B"], 20, tmp_path / "linked.pdb")

    assert position_map == [("A", 3), ("A", 5), ("A", 6), *[None] * 20, ("B", 2), ("B", 3)]
    residues = read_pdb_residues(tmp_path / "linked.pdb")
    assert {chain for chain, _, _ in residues} == {LINKED_CHAIN}
    numbers = [resseq for _, resseq, _ in residues]
    # A keeps its internal gap, the linker follows A, and B follows the linker.
    assert numbers == [1, 3, 4, *range(5, 25), 25, 26]
    linker = residues[3:23]
    assert all(lines[0][17:20] == "GLY" and len(lines) == 4 for _, _, lines in linker)
    # Linker CA atoms lie between A's last CA and B's first CA.
    linker_x = [float(lines[1][30:38]) for _, _, lines in linker]
    assert 7.6 < min(linker_x) and max(linker_x) < 30.0


@pytest.mark.parametrize("renumbered", [False, True])
def test_unlinking_restores_chains_numbering_and_drops_linker(tmp_path, renumbered):
    source = two_chain_pdb(tmp_path / "model.pdb")
    position_map = write_linked_pdb(source, ["A", "B"], 3, tmp_path / "linked.pdb")
    # ColabDesign keeps the linked numbering; also accept an output renumbered 1..L.
    scored, serial = [], 1
    for position, (_, resseq, lines) in enumerate(read_pdb_residues(tmp_path / "linked.pdb"), start=1):
        for line in lines:
            number = position if renumbered else resseq
            scored.append(f"ATOM  {serial:5d}{line[11:22]}{number:4d}{line[26:]}\n")
            serial += 1
    (tmp_path / "scored_linked.pdb").write_text("".join(scored) + "END\n")

    unlink_pdb(tmp_path / "scored_linked.pdb", position_map, tmp_path / "scored.pdb")
    text = (tmp_path / "scored.pdb").read_text()
    residues = read_pdb_residues(tmp_path / "scored.pdb")
    assert [(chain, resseq) for chain, resseq, _ in residues] == [("A", 3), ("A", 5), ("A", 6), ("B", 2), ("B", 3)]
    assert "GLY" not in text
    assert text.count("TER") == 2 and text.rstrip().endswith("END")


def test_linker_is_dropped_from_scores_and_chain_ids(tmp_path):
    position_map = [("A", 1), ("A", 2), None, None, ("B", 1)]
    keep = np.flatnonzero([entry is not None for entry in position_map])
    logits = np.arange(5 * 5 * 2, dtype=float).reshape(5, 5, 2)
    outputs = {
        "predicted_aligned_error": {"logits": logits, "breaks": np.array([1.0])},
        "distogram": {"logits": logits, "bin_edges": np.array([2.0])},
    }
    trimmed = drop_linker_outputs(outputs, keep)
    assert trimmed["predicted_aligned_error"]["logits"].shape == (3, 3, 2)
    np.testing.assert_array_equal(trimmed["distogram"]["logits"][2, 0], logits[4, 0])
    np.testing.assert_array_equal(trimmed["predicted_aligned_error"]["breaks"], [1.0])
    np.testing.assert_array_equal(linked_asym_id(position_map), [0, 0, 1])


def test_linker_metrics_match_unlinked_metrics_on_the_same_residues():
    pytest.importorskip("jax")
    from af2rank_pipeline.actifptm_runner import compute_pair_metrics, linked_confidence

    rng = np.random.default_rng(0)
    full = 30 + 5
    keep = np.r_[0:12, 17:35]
    asym_id = np.r_[np.zeros(12), np.ones(18)].astype(int)
    pae_logits = rng.normal(size=(full, full, 64))
    dist_logits = rng.normal(size=(full, full, 64))
    outputs = {
        "predicted_aligned_error": {"logits": pae_logits, "breaks": np.linspace(0, 31, 63)},
        "distogram": {"logits": dist_logits, "bin_edges": np.linspace(2.3125, 21.6875, 63)},
    }
    aux = {"plddt": rng.uniform(size=full), "pae": rng.uniform(0, 30, (full, full))}
    trimmed = drop_linker_outputs(outputs, keep)
    confidence = linked_confidence(aux, trimmed, keep, asym_id)
    reference = {name: {**outputs[name], "logits": outputs[name]["logits"][keep][:, keep]} for name in outputs}

    assert confidence["plddt"] == pytest.approx(aux["plddt"][keep].mean())
    trimmed_pae = aux["pae"][keep][:, keep]
    np.testing.assert_array_equal(confidence["pae"], trimmed_pae)
    # Like ColabDesign's logged PAE: the mean over within-chain pairs only.
    same_chain = asym_id[:, None] == asym_id[None, :]
    assert confidence["pae_mean"] == pytest.approx(trimmed_pae[same_chain].mean())
    assert 0 < confidence["ptm"] < 1 and 0 < confidence["i_ptm"] < 1
    metrics = compute_pair_metrics(trimmed, asym_id, ["A", "B"])
    assert metrics == compute_pair_metrics(reference, asym_id, ["A", "B"])


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["--model-mode", "alphafold-multimer"], 20),
        (["--model-mode", "alphafold-multimer", "--multimer-linker", "0"], 0),
        (["--model-mode", "alphafold-multimer", "--multimer-linker", "12"], 12),
        ([], 0),
    ],
)
def test_multimer_scoring_links_chains_by_default(argv, expected):
    args = cli.build_parser().parse_args(
        ["af2rank", "--target-spec", "t.json", "--out", "o", "--params", "p", *argv])
    assert cli._af2rank_config(args).linker_length == expected


def test_linker_is_rejected_for_ptm_weights_before_cleaning(tmp_path, capsys):
    status = cli.main([
        "run-target", "--target", "1brs", "--fasta", str(REPO / "examples/1brs/target.fasta"),
        "--models", str(REPO / "examples/1brs/model.pdb"), "--out", str(tmp_path / "out"),
        "--params", str(tmp_path / "params"), "--multimer-linker", "20",
    ])
    assert status == 1
    assert "--multimer-linker requires --model-mode alphafold-multimer" in capsys.readouterr().err
    assert not (tmp_path / "out" / "cleaned").exists()


def test_linked_runs_get_their_own_output_folder():
    linked = AF2RankConfig(model_mode="alphafold-multimer", version="v3", model_num=1, linker_length=20)
    separate = AF2RankConfig(model_mode="alphafold-multimer", version="v3", model_num=1)
    assert model_subdir(linked) == "alphafold-multimer_v3_model1_rec1_link20"
    assert model_subdir(separate) == "alphafold-multimer_v3_model1_rec1"
