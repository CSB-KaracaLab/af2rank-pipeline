"""Shared loader for the CAPRI benchmark tables (run the scripts from benchmarks/capri)."""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]


def load(setting="fullres"):
    """One row per decoy x weight. setting: 'fullres' (all 153 runs) or 'masked' (all 153 runs, the paper's setting).

    Adds: good (CAPRI acceptable or better), dockq / main_classification (CAPRI's own values), neg_pae_contacts.
    Decoys whose input has no interface (11 of T041/T131, all CAPRI-incorrect) get DockQ_io and fnat_io 0 and the
    maximum PAE (31.75 A).
    """
    f = {"fullres": "scores_capri_fullres_all_weights.csv", "masked": "scores_capri_masked_all_weights.csv"}[setting]
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


def targets():
    return pd.read_csv(HERE / "data/capri_targets.csv").set_index("tid")
