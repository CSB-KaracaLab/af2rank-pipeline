#!/bin/bash
# Build the pipeline environment on a MareNostrum 5 login node, which has no internet.
# 1. On a machine with internet, download the Linux wheels listed in requirements-lock.txt
#    (README.md shows the command) and copy them to $AF2RANK_ROOT/wheels together with the
#    ColabDesign and af2rank_pipeline wheels.
# 2. Run this script on alogin1/alogin2:  AF2RANK_ROOT=/gpfs/projects/<project>/<user>/af2rank bash setup_env.sh
set -euo pipefail
ROOT=${AF2RANK_ROOT:?set AF2RANK_ROOT to your work folder}
module purge
module load intel impi mkl hdf5 python/3.11.5-gcc
# The module sets PYTHONHOME/PYTHONPATH to its own packages; the venv must not see them.
unset PYTHONHOME PYTHONPATH
cd "$ROOT"
python3 -m venv venv
venv/bin/pip install -q --no-index --find-links wheels -r requirements-lock.txt
venv/bin/pip install -q --no-index --find-links wheels --no-deps colabdesign af2rank_pipeline
venv/bin/python -c "import jax, colabdesign, af2rank_pipeline, DockQ; print('jax', jax.__version__)"
