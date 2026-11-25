# AF2Rank Pipeline

AF2Rank Pipeline scores existing protein-complex structures. It standardizes chains and residues against a target FASTA, runs AF2Rank through ColabDesign, calculates actifpTM and interface metrics, and combines the results in a CSV for model comparison.

The [report](docs/af2rank_all_weights.md) explains how the 17 template-capable AlphaFold2 weights and the multimer glycine linker are used. Worked examples and benchmarks are listed under [Examples and benchmarks](#examples-and-benchmarks); the scores they report are defined under [Scores and settings](#scores-and-settings).

## Try it in Colab

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CSB-KaracaLab/af2rank-pipeline/blob/main/notebooks/colab_demo.ipynb)

The notebook opens on a GPU runtime; if Colab gives you a CPU runtime, choose *Runtime → Change runtime type →* a GPU (for example T4) before running. Run the included 1BRS example or upload your own FASTA and structures.

Tick any of these template-capable AlphaFold2 weight sets to score the models:

| Weight set | Weights |
|---|---|
| AF2-ptm 1, 2 (default) | `model_1_ptm`, `model_2_ptm` |
| AF2.1 multimer v1 | `model_1_multimer` … `model_5_multimer` |
| AF2.2 multimer v2 | `model_1_multimer_v2` … `model_5_multimer_v2` |
| AF2.3 multimer v3 | `model_1_multimer_v3` … `model_5_multimer_v3` |

Multimer runs join the chains with a 20-glycine linker. The linker's template is masked, and the linker is left out of every score.

Every model is scored with each selected weight. The notebook then shows the scores for each model and weight, the rescored structures and PAE maps, and a ZIP of all results.

## Run locally

Install Git and [uv](https://docs.astral.sh/uv/getting-started/installation/), then run these commands from a clone of this repository. The project uses Python 3.10 and a committed lockfile:

```bash
uv python install 3.10
uv sync --locked --all-extras
source .venv/bin/activate
```

The lockfile installs CPU JAX. A GPU is **not required** for the bundled single-model example: the full 199-residue 1BRS run has been tested on CPU. Larger structures or many models may take much longer; use a GPU for those workloads.

Full scoring needs AlphaFold model parameters, but not AlphaFold's sequence databases. Download the default model-2 pTM weight (373 MB) from the [official parameter archive](https://github.com/google-deepmind/alphafold#model-parameters):

```bash
python env/download_default_params.py params
```

The helper downloads just the required file and verifies its checksum. `params/` is excluded from version control. If you already have a parameter directory, use it instead. Other weights are fetched by name, for example `python env/download_default_params.py params --weight model_1_ptm --weight model_1_multimer_v2`; the [official download script](https://github.com/google-deepmind/alphafold/blob/main/scripts/download_alphafold_params.sh) provides the complete parameter bundle.

AF2Rank passes each input structure to AlphaFold as a template, so only template-capable weights apply:

- **Weights:** AF2-ptm models 1 and 2 (`--model-num 1` or `2`), and AF2-multimer models 1–5 from AF2.1, AF2.2 and AF2.3 (`--model-mode alphafold-multimer --multimer-version v1`, `v2` or `v3`, default v3). Each needs its weight file (`params_model_N_multimer.npz`, `_multimer_v2.npz` or `_multimer_v3.npz`), which `--weight` above fetches.
- **Glycine linker:** AF2-multimer was trained with per-chain templates only, so multimer runs join the chains into one chain with a 20-glycine linker and the whole complex becomes a single intra-chain template. The linker residues carry no template information. They are removed from the scored PDB, the PAE matrix and every metric, and chains and residue numbers are restored.
- **Separate chains instead:** `--multimer-linker N` sets the linker length; `--multimer-linker 0` gives the chains as separate templates, as in the original af2rank-pipeline.
- **Linker outputs:** linked runs keep the joined input (`<model>.pdb`), AlphaFold's output with the linker (`<model>_scored_linked.pdb`) and the linker residue numbers (`<model>_linker.json`) under `af2rank/<run>/linked_inputs/`. The Colab structure viewer can show that output with the linker highlighted.

Run the included 1BRS complex:

```bash
af2rank-pipeline run-target \
  --target 1brs \
  --fasta examples/1brs/target.fasta \
  --models examples/1brs/model.pdb \
  --out runs/1brs \
  --params params
```

The result table is `runs/1brs/manifests/final_models.csv`; the scored structure and PAE data are under `runs/1brs/af2rank/`. The [example notes](examples/1brs/README.md) describe the source of the bundled data. Use a new output directory for each independent run.


### NVIDIA GPU on Linux

For an NVIDIA GPU with a compatible CUDA 12 driver, install the matching JAX wheel **after** `uv sync`:

```bash
uv pip install --python .venv/bin/python \
  --constraint env/cuda12-constraints.txt \
  --find-links https://storage.googleapis.com/jax-releases/jax_cuda_releases.html \
  'jax[cuda12_pip]==0.4.20'
python -c "import jax; print(jax.devices())"
```

Check that JAX reports a GPU, then run `af2rank-pipeline` from the activated environment as above. Running `uv sync` or `uv run` afterward can restore the locked CPU wheel. This CUDA overlay has completed the bundled example in hosted Colab; GPU drivers and libraries vary across machines. Follow [JAX's accelerator guidance](https://docs.jax.dev/en/latest/installation.html) if device detection fails. On a shared cluster, use the site's approved environment and run inference inside a scheduled GPU allocation.

## Use your own models

Provide one FASTA record per physical protein chain. Homomer copies need separate records with unique headers; record order determines the canonical chain order. Use standard amino-acid letters.

```fasta
>chain_A
ACDEFGHIKLMNPQRSTVWY
>chain_B
ACDEFGHIKLMNPQRSTVWY
```

`--models` accepts a PDB or mmCIF file, a gzipped PDB, a multi-model PDB or mmCIF ensemble, or a directory of structures. Each ensemble member is scored separately. Nonprotein components are not retained during cleaning. To score your own models:

```bash
af2rank-pipeline run-target \
  --target my_complex \
  --fasta /path/to/target.fasta \
  --models /path/to/models \
  --out runs/my_complex \
  --params params
```

Models that fail cleaning are logged in `logs/failures.jsonl` and skipped; the remaining models are scored. Use `--strict-clean` with `run-target` or `clean` to stop if any model fails cleaning. A run stops if no models can be cleaned, or if scoring is incomplete for the cleaned models. The Colab notebook uses `--strict-clean`. Malformed ensemble markers are rejected before cleaning.

The default run uses pTM model 2, one recycle, and one inference iteration. It masks template sequence and side chains. Use `af2rank-pipeline run-target --help` for other model and scoring settings. `--tm-exec` optionally points to a working [TMscore](https://zhanggroup.org/TM-score/) executable; without it, TM-score fields and the composite score are empty.

## Results

| Path in the output directory | Contents |
| --- | --- |
| `target_spec.json` | Target chains and sequences |
| `cleaned/af2rank_input/` | Standardized input structures |
| `manifests/cleaning.jsonl` | Chain alignment, coverage, and source provenance |
| `af2rank/` | AF2Rank scores, scored PDBs, and PAE data |
| `additional_metrics/` | Input–output DockQ and interface PAE metrics |
| `manifests/final_models.csv` | Combined scores and provenance for each model |
| `logs/failures.jsonl` | Stage errors, if any |

DockQ uses the cleaned input as its reference, so a high value means the interface remained similar after AF2Rank inference; it does not establish docking accuracy (see [Scores and settings](#scores-and-settings)). Check the output row count and any failure log before comparing models.

## Examples and benchmarks

Each folder except `examples/1brs` (the quick-start run, default weight only) reports scores for all 17 template-capable weights. [Scores and settings](#scores-and-settings) defines the two settings.

| Folder | What | Setting | Size |
| --- | --- | --- | --- |
| [`examples/1brs`](examples/1brs) | barnase–barstar, the quick-start example | masked | 1 model |
| [`examples/2ptc`](examples/2ptc) | trypsin–BPTI crystal structure with waters and ions to clean | masked | 1 model |
| [`examples/e2a_hpr`](examples/e2a_hpr) | E2A–HPr NMR ensembles | masked | 3 and 30 models |
| [`benchmarks/nmr`](benchmarks/nmr) | curated NMR dimers: homodimer, heterodimer and flexible sets | full resolution | 65 dimers |
| [`benchmarks/xray`](benchmarks/xray) | unseen X-ray dimers matched to the NMR set | heterodimers: both settings; homodimers: not scored yet | 53 hetero, 39 homo |
| [`benchmarks/capri`](benchmarks/capri) | CAPRI docking decoys: which score and ranking strategy separate good models from bad ones (clustering, linker vs separate chains) | both settings; masked also with separate chains | 9 targets, 1,340 decoys |
| [`cluster/marenostrum5`](cluster/marenostrum5) | Slurm scripts for H100 clusters | – | – |

## Scores and settings

**Settings.**
- **Masked** (the pipeline default and the AF2Rank paper's setting): the template's sequence and side chains are hidden from AlphaFold, 1 recycle. AlphaFold sees only the backbone, so it can move a wrong interface away.
- **Full resolution:** sequence and side chains kept, 0 recycles (`--no-mask-sequence --no-mask-sidechains --recycles 0`). AlphaFold tends to copy the template, right or wrong.

**AlphaFold confidence** (higher is better): pLDDT (per residue), pTM (whole complex), ipTM (interface), actifpTM (ipTM restricted to the residues actually at the interface; Varga et al., 2025). PAE is AlphaFold's predicted aligned error in Å (lower is better); *PAE at input contacts* averages it over the residue pairs in contact in the input.

**How much AlphaFold kept the input** (the input is the reference, not the true structure): *DockQ input→output* (DockQ_io in tables) is the DockQ of AlphaFold's output against the cleaned input; below 0.23 the input interface counts as lost. *fnat input→output* is the fraction of input interface contacts kept in the output, and *iRMSD input→output* the interface RMSD in Å. When the input is the native structure, as in the X-ray benchmark, DockQ input→output is the true DockQ of the output.

**Other interface scores** used in the CAPRI benchmark: ipSAE, LIS and minimum inter-chain PAE (built from the inter-chain PAE), and pDockQ (built from pLDDT and contact count).

**Decoy benchmarks.** *CAPRI classes* are CAPRI's quality classes from fnat, ligand RMSD and interface RMSD to the native: incorrect, acceptable, medium, high. A *good* decoy is acceptable or better. *AUC* is the ROC AUC for good vs incorrect decoys within one target and weight (0.5 = random). *Top-1* or *top-5 success* means a good decoy is ranked first or within the first five.

## Other workflows

`af2rank-pipeline clean` runs the input-cleaning stage without model parameters or inference. The CLI also provides staged commands for batching, AF2Rank inference, DockQ, aggregation, and cleanup; use `af2rank-pipeline --help` and each subcommand's `--help`. `--resume` can reuse completed cleaning and AF2Rank results when the model set and scoring configuration match; start a new output directory if inputs or weights change, because resume does not hash their contents.

To check a local installation or contribute changes:

```bash
uv run --locked --all-extras python -m pytest -q
```

## License

The project-owned code is [MIT licensed](LICENSE). See [third-party components and citations](THIRD_PARTY.md) for bundled code, example data, and external dependencies.
