# Running AF2Rank on MareNostrum 5 (H100)

These scripts ran the E2A–HPr example (`examples/e2a_hpr`, 3 NMR models × 17 weights) in 16 min on one H100 of the MareNostrum 5 accelerated partition. They should carry over to other Slurm clusters with small changes to the module lines.

| File | What it does |
|---|---|
| `requirements-lock.txt` | The `uv.lock` versions for Linux and Python 3.11, plus the CUDA 12 JAX wheels the Colab notebook uses. `nvidia-cuda-nvcc-cu12` is pinned to 12.3.107. |
| `setup_env.sh` | Builds the virtual environment offline from a folder of wheels (login nodes have no internet). |
| `link_weights.sh` | Links AF2-ptm 1–2 and AF2.3 multimer 1–5 from BSC's AlphaFold module; AF2.1 and AF2.2 multimer weights are copied in separately. |
| `run_target.slurm` | Scores one target with all 17 weights, one after another, then writes a summary. |
| `summarize.py` | Collects every weight's `final_models.csv` into `scores_all_weights.csv` and draws ipTM against DockQ input→output. |

## Steps

1. **Download the wheels** on a machine with internet, from the repository root:

   ```bash
   P=(); for t in 2_28 2_27 2_26 2_25 2_24 2_17 2_12 2_5; do P+=(--platform manylinux_${t}_x86_64); done
   grep -v '^colabdesign' examples/marenostrum5/requirements-lock.txt > req-binary.txt
   python3 -m pip download --no-deps --only-binary=:all: --python-version 3.11 --implementation cp \
     --abi cp311 --abi abi3 --abi none "${P[@]}" --platform manylinux2014_x86_64 --platform any \
     -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html -r req-binary.txt -d wheels
   python3 -m pip wheel --no-deps "git+https://github.com/sokrypton/ColabDesign.git@094e2cb3603dee7d99846e0977736bd943c830c2" -w wheels
   uv build --wheel -o wheels
   ```

2. **Copy** `wheels/` and these scripts to your work folder on MareNostrum (`alogin1.bsc.es`), then run `AF2RANK_ROOT=<work folder> bash setup_env.sh`.
3. **Weights:** `AF2RANK_PARAMS=<params folder> bash link_weights.sh`, then copy `params_model_{1..5}_multimer.npz` and `params_model_{1..5}_multimer_v2.npz` into the same folder. The notebook's weight cell lists their SHA-256 sums.
4. **Submit** (the header of `run_target.slurm` shows the full command), with your project as `--account`.

## H100 notes

- JAX 0.4.20 fails on H100 with `nvlink fatal: cannot link accelerated arch to non-accelerated`. `run_target.slurm` sets `XLA_FLAGS=--xla_gpu_force_compilation_parallelism=1`, which avoids the failing link step. The pipeline keeps this flag even though ColabDesign overwrites `XLA_FLAGS` on import.
- The job starts with a small GPU compile test, so a broken setup fails within seconds instead of after the first weight.
- `acc_debug` allows one job of up to 2 h per user; use `acc_ehpc` for longer batches.
