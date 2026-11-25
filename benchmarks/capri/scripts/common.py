"""Shared loader for the CAPRI benchmark tables (run the scripts from benchmarks/capri)."""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]


def load(setting="fullres"):
    """One row per decoy x weight. setting: 'fullres' (all 153 runs), 'masked' (all 153 runs, the paper's setting,
    20-glycine linker for multimer weights) or 'separate' (masked, multimer weights given as separate chains, 135 runs).

    Adds: good (CAPRI acceptable or better), dockq / main_classification (CAPRI's own values), neg_pae_contacts.
    Decoys whose input has no interface (11 of T041/T131, all CAPRI-incorrect) get DockQ_io and fnat_io 0 and the
    maximum PAE (31.75 A).
    """
    f = {"fullres": "scores_capri_fullres_all_weights.csv", "masked": "scores_capri_masked_all_weights.csv",
         "separate": "scores_capri_masked_separate_chains.csv"}[setting]
    d = pd.read_csv(HERE / "results" / f)
    d["target"] = d.target.str.replace(".", "_", regex=False)
    d = d.rename(columns={"capri_dockq": "dockq", "capri_class": "main_classification"})
    d["good"] = d.main_classification.ne("incorrect")
    d["dockq_io"] = d.dockq_io.fillna(0.0); d["fnat_io"] = d.fnat_io.fillna(0.0)
    d["irmsd_io"] = d.irmsd_io.fillna(d.irmsd_io.max())
    d["pae_contacts_mean"] = d.pae_contacts_mean.fillna(31.75)
    d["neg_pae_contacts"] = -d.pae_contacts_mean
    d["family"] = d.weight.str.replace(r"model_\d_", "", regex=True)
    return d


VERSION = {"ptm": "AF2-ptm", "multimer": "AF2.1", "multimer_v2": "AF2.2", "multimer_v3": "AF2.3"}


def load_pair():
    """Masked multimer runs in both input modes, one row per decoy x weight (135 runs): actifpTM with the 20-glycine
    linker (actifptm_linker), with separate chains (actifptm_separate), their product, DockQ input->output in both
    modes, CAPRI DockQ and class, good, version."""
    keep = ["target", "weight", "model_id", "actifptm", "i_ptm", "dockq_io"]
    lk = load("masked"); lk = lk[lk.weight.str.contains("multimer")]
    sp = load("separate")
    d = lk[keep + ["dockq", "main_classification", "good", "family"]].merge(
        sp[keep], on=["target", "weight", "model_id"], suffixes=("_linker", "_separate"))
    d["product"] = d.actifptm_linker * d.actifptm_separate
    d["version"] = d.family.map(VERSION)
    return d


def targets():
    return pd.read_csv(HERE / "data/capri_targets.csv").set_index("tid")
