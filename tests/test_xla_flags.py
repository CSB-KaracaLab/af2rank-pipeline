import importlib.util
import os
import subprocess
import sys

import pytest

from af2rank_pipeline.actifptm_runner import merge_xla_flags


def test_merge_keeps_every_flag_and_last_value_wins():
    assert merge_xla_flags("--a=1 --b=2", "--b=3 --c") == "--a=1 --b=3 --c"
    assert merge_xla_flags("", "") == ""


@pytest.mark.skipif(importlib.util.find_spec("colabdesign") is None, reason="needs the inference extra")
def test_caller_xla_flags_survive_colabdesign_import():
    # ColabDesign overwrites XLA_FLAGS when first imported; H100 runs need their own flag kept.
    # A fresh interpreter makes sure the import (and its overwrite) really happens.
    code = ("import os; from af2rank_pipeline.actifptm_runner import setup_imports; "
            "setup_imports(None); print(os.environ['XLA_FLAGS'])")
    env = {**os.environ, "XLA_FLAGS": "--xla_gpu_force_compilation_parallelism=1"}
    flags = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True,
                           check=True).stdout.split()
    assert "--xla_gpu_force_compilation_parallelism=1" in flags
    assert "--xla_gpu_enable_triton_gemm=false" in flags
