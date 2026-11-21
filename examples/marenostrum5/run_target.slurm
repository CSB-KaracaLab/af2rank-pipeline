#!/bin/bash
#SBATCH --job-name=af2rank
#SBATCH --qos=acc_debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# AF2Rank on one target with all 17 template-capable AF2 weights, one weight after another.
# Same command per weight as the Colab notebook's scoring cell. Submit with, for example:
#   sbatch --account=<project> --export=ALL,AF2RANK_ROOT=$PWD,AF2RANK_PARAMS=$PWD/../af2rank-weights/params,\
#TARGET=e2a_hpr,FASTA=$PWD/inputs/e2a_hpr/target.fasta,MODELS=$PWD/inputs/e2a_hpr/models run_target.slurm
set -euo pipefail
ROOT=${AF2RANK_ROOT:?}; PARAMS=${AF2RANK_PARAMS:?}; TARGET=${TARGET:?}; FASTA=${FASTA:?}; MODELS=${MODELS:?}
OUT=$ROOT/results/${TARGET}-${SLURM_JOB_ID}
module purge
module load intel impi mkl hdf5 python/3.11.5-gcc
unset PYTHONHOME PYTHONPATH
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
# JAX 0.4.20 on H100: XLA's parallel compile links with nvlink, which fails for sm_90a
# ("cannot link accelerated arch to non-accelerated"). One compile thread avoids the link step.
export XLA_FLAGS=--xla_gpu_force_compilation_parallelism=1
CLI=$ROOT/venv/bin/af2rank-pipeline
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
# Fail fast if XLA cannot compile for this GPU.
$ROOT/venv/bin/python -c "import jax, jax.numpy as jnp; x = jnp.ones((256, 256)); print('JAX', jax.__version__, jax.devices(), 'GPU compile test:', float(jnp.concatenate([x @ x, x]).sum()))"

WEIGHTS=()
for n in 1 2; do WEIGHTS+=("model_${n}_ptm|--model-num $n"); done
for spec in "v1:multimer" "v2:multimer_v2" "v3:multimer_v3"; do
  version=${spec%%:*}; suffix=${spec#*:}
  for n in 1 2 3 4 5; do
    WEIGHTS+=("model_${n}_${suffix}|--model-mode alphafold-multimer --multimer-version $version --model-num $n --multimer-linker 20")
  done
done
for entry in "${WEIGHTS[@]}"; do
  name=${entry%%|*}; args=${entry#*|}
  start=$(date +%s)
  echo "== $name"
  $CLI run-target --target "$TARGET" --fasta "$FASTA" --models "$MODELS" --params "$PARAMS" \
       --workers 1 --strict-clean $args --out "$OUT/$name"
  echo "== $name done in $(( $(date +%s) - start )) s"
done
$ROOT/venv/bin/python $ROOT/summarize.py "$OUT"
echo "results: $OUT"
