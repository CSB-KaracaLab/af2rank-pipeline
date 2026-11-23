import pickle

import numpy as np
import pytest

from af2rank_pipeline.msa_input import GAP, RESTYPES, build_msa_inputs, column_map, load_msa_features


def _features(tmp_path, seq_a="ACDEFG", seq_b="KLMN", rows=None):
    """A small AF-multimer-style features.pkl: query row plus `rows` (strings over both chains, '-' = gap)."""
    aa = {c: i for i, c in enumerate(RESTYPES)}
    query = seq_a + seq_b
    rows = [query] + (rows or [])
    msa = np.array([[aa.get(c, GAP) for c in r] for r in rows])
    deletion = np.zeros(msa.shape)
    if len(rows) > 1:
        deletion[1, 0] = 3
    feats = {"msa": msa, "deletion_matrix": deletion, "aatype": msa[0],
             "asym_id": np.array([1] * len(seq_a) + [2] * len(seq_b))}
    path = tmp_path / "features.pkl"
    path.write_bytes(pickle.dumps(feats))
    return load_msa_features(path)


def test_chains_split_by_asym_id(tmp_path):
    f = _features(tmp_path)
    assert [c["sequence"] for c in f["chains"]] == ["ACDEFG", "KLMN"]


def test_column_map_with_linker_and_missing_residue(tmp_path):
    f = _features(tmp_path)
    # chain A misses its first residue; two linker positions; chain B complete
    positions = [("A", c) for c in "CDEFG"] + [None, None] + [("B", c) for c in "KLMN"]
    cols = column_map(positions, ["A", "B"], f)
    assert cols.tolist() == [1, 2, 3, 4, 5, -1, -1, 6, 7, 8, 9]


def test_column_map_rejects_wrong_sequence(tmp_path):
    f = _features(tmp_path)
    with pytest.raises(ValueError):
        column_map([("A", c) for c in "WWWWWW"] + [("B", c) for c in "KLMN"], ["A", "B"], f)


def test_build_msa_inputs_shapes_and_features(tmp_path):
    f = _features(tmp_path, rows=["ACDEFGKLMN", "ACDEFG----", "----------", "ACDEWGKLMW"])
    positions = [("A", c) for c in "ACDEFG"] + [None] + [("B", c) for c in "KLMN"]
    cols = column_map(positions, ["A", "B"], f)
    out = build_msa_inputs(f, cols, n_clusters=2, n_extra=10, multimer=True)
    assert out["n_rows"] == 4  # the all-gap row is dropped
    assert out["msa_feat"].shape == (2, 11, 49)
    assert out["extra_msa"].shape == (2, 11)
    # linker column is a gap in every row
    assert out["msa_feat"][:, 6, GAP].tolist() == [1.0, 1.0]
    assert (out["extra_msa"][:, 6] == GAP).all()
    # row 1 has one deletion at column 0: has_deletion = 1, transformed value = arctan(1) * 2/pi = 0.5
    assert out["msa_feat"][1, 0, 23] == 1.0
    assert out["msa_feat"][1, 0, 24] == pytest.approx(0.5)
    # profiles are probability vectors
    assert np.allclose(out["msa_feat"][:, :, 25:48].sum(-1), 1.0)
    # multimer keeps raw extra deletions; ptm transforms them
    ptm = build_msa_inputs(f, cols, n_clusters=1, n_extra=10, multimer=False)
    assert ptm["extra_deletion_value"].max() == pytest.approx(0.5)
    assert out["extra_deletion_value"].max() == 0.0 or out["extra_deletion_value"].max() == 3.0
