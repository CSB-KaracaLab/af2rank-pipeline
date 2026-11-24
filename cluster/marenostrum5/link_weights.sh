#!/bin/bash
# AF2-ptm 1-2 and AF2.3 multimer 1-5 are provided by BSC's AlphaFold 2.3.2 module (their SHA-256
# sums match the notebook's). Copy the AF2.1 and AF2.2 multimer weights into $AF2RANK_PARAMS yourself.
set -euo pipefail
SRC=/gpfs/apps/MN5/ACC/ALPHAFOLD/SRC/database/Alphafold/v2.3.2/params
DEST=${AF2RANK_PARAMS:?set AF2RANK_PARAMS to your params folder}
mkdir -p "$DEST"
for w in model_1_ptm model_2_ptm model_{1,2,3,4,5}_multimer_v3; do
  ln -sfn "$SRC/params_$w.npz" "$DEST/params_$w.npz"
done
ls "$DEST"
