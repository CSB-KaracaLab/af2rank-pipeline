import hashlib
import importlib.util
import io
import tarfile
import zipfile
from pathlib import Path

import pytest

from af2rank_pipeline import cli
from af2rank_pipeline.actifptm_runner import model_name_from_args
from af2rank_pipeline.af2rank_runner import AF2RankConfig, model_subdir, run_af2rank_stage

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("version", "name"),
    [("v1", "model_3_multimer"), ("v2", "model_3_multimer_v2"), ("v3", "model_3_multimer_v3")],
)
def test_multimer_versions_map_to_release_weight_names(version, name):
    assert model_name_from_args("alphafold-multimer", 3, version) == name
    assert model_name_from_args("alphafold", 2, version) == "model_2_ptm"


def test_cli_selects_multimer_version_and_keeps_runs_apart():
    args = cli.build_parser().parse_args([
        "af2rank", "--target-spec", "t.json", "--out", "o", "--params", "p",
        "--model-mode", "alphafold-multimer", "--multimer-version", "v1", "--model-num", "4",
    ])
    cfg = cli._af2rank_config(args)
    assert (cfg.version, cfg.model_num, cfg.linker_length) == ("v1", 4, 20)
    assert model_subdir(cfg) == "alphafold-multimer_v1_model4_rec1_link20"
    default = cli._af2rank_config(cli.build_parser().parse_args(
        ["af2rank", "--target-spec", "t.json", "--out", "o", "--params", "p"]))
    assert default.version == "ptm"


def test_stage_rejects_unknown_multimer_version(tmp_path):
    from af2rank_pipeline.target import load_target_spec

    spec = load_target_spec("1brs", REPO / "examples/1brs/target.fasta")
    cfg = AF2RankConfig(model_mode="alphafold-multimer", version="v4", params=str(tmp_path))
    with pytest.raises(ValueError, match="v1 or v2 or v3"):
        run_af2rank_stage(spec, tmp_path / "out", cfg)


def load_helper():
    spec = importlib.util.spec_from_file_location("weights", REPO / "env/download_default_params.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("member_prefix", ["", "./"])
def test_helper_fetches_named_weight_from_either_archive_layout(tmp_path, monkeypatch, member_prefix):
    weight = io.BytesIO()
    with zipfile.ZipFile(weight, "w") as archive:
        archive.writestr("w.npy", b"weights")
    archive_bytes = io.BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w", format=tarfile.USTAR_FORMAT) as tar:
        for name, data in (("LICENSE", b"license"), ("params_model_1_multimer.npz", weight.getvalue())):
            info = tarfile.TarInfo(member_prefix + name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    data = archive_bytes.getvalue()

    helper = load_helper()
    requested = []
    monkeypatch.setattr(helper, "get_range", lambda start, end, url: (
        requested.append(url) or io.BytesIO(data[start:end + 1])))
    url = helper.archive_url("2021-10-27")
    sha = hashlib.sha256(weight.getvalue()).hexdigest()
    path = helper.download_weight(tmp_path, "params_model_1_multimer.npz", sha, url)

    assert path.read_bytes() == weight.getvalue()
    assert set(requested) == {url}
    with pytest.raises(RuntimeError, match="Unexpected checksum"):
        helper.download_weight(tmp_path, "params_model_1_multimer.npz", "0" * 64, url)


def test_seed_and_dropout_reach_the_inference_script(tmp_path, monkeypatch):
    from af2rank_pipeline import af2rank_runner
    from af2rank_pipeline.target import load_target_spec

    args = cli.build_parser().parse_args([
        "af2rank", "--target-spec", "t.json", "--out", "o", "--params", "p", "--seed", "7", "--dropout"])
    cfg = cli._af2rank_config(args)
    assert (cfg.seed, cfg.dropout) == (7, True)
    assert model_subdir(cfg) == "alphafold_ptm_model2_rec1_seed7_dropout"
    assert model_subdir(AF2RankConfig()) == "alphafold_ptm_model2_rec1"

    commands = []
    monkeypatch.setattr(af2rank_runner.subprocess, "run", lambda cmd, check: commands.append(cmd))
    spec = load_target_spec("1brs", REPO / "examples/1brs/target.fasta")
    af2rank_runner._run_actifptm_script(spec, tmp_path, tmp_path, tmp_path, cfg)
    af2rank_runner._run_actifptm_script(spec, tmp_path, tmp_path, tmp_path, AF2RankConfig(params="p"))
    seeded, default = commands
    assert seeded[seeded.index("--seed") + 1] == "7" and "--dropout" in seeded
    assert "--seed" not in default and "--dropout" not in default
