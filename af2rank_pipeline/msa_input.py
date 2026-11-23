"""MSA input for AF2Rank: turn AlphaFold's own MSA features into ColabDesign inputs.

AF2Rank normally runs on a single sequence. With an MSA, the scored model still enters as a template, and the
MSA comes from AlphaFold 2.3's data pipeline (`features.pkl` from pipeline_multimer: paired + unpaired rows,
already in AlphaFold's residue order, 0-19 amino acids, 20 unknown, 21 gap).

The MSA columns are mapped onto the residues ColabDesign scores (residues present in the model, plus any linker
positions, which get gaps). Rows are split as in AlphaFold: the first `n_clusters` rows are cluster centres
(row 0 is the query), the next `n_extra` rows form the extra MSA, and each extra row adds to the profile of the
centre it agrees with most. Unlike AlphaFold, rows are taken in order rather than sampled, and no positions
are masked, so a run is deterministic.
"""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

RESTYPES = "ARNDCQEGHILKMFPSTWYV"
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
GAP = 21


def load_msa_features(path: str | Path) -> dict:
    """Load features.pkl and return the MSA, deletions and the per-chain query sequences."""
    with open(path, "rb") as handle:
        feats = pickle.load(handle)
    msa = np.asarray(feats["msa"], dtype=np.int32)
    deletion = np.asarray(feats["deletion_matrix"], dtype=np.float32)
    asym = np.asarray(feats["asym_id"], dtype=np.int32)
    aatype = np.asarray(feats["aatype"], dtype=np.int32)
    if aatype.ndim == 2:  # one-hot in some AlphaFold versions
        aatype = aatype.argmax(-1)
    chains = []
    for a in sorted(set(asym.tolist())):
        cols = np.flatnonzero(asym == a)
        chains.append({"columns": cols, "sequence": "".join(RESTYPES[i] if i < 20 else "X" for i in aatype[cols])})
    return {"msa": msa, "deletion": deletion, "chains": chains}


def align_positions(model_seq: str, query_seq: str) -> list[int | None]:
    """For each residue of the model sequence, the index of the aligned query residue (None if unaligned)."""
    if model_seq == query_seq:
        return list(range(len(model_seq)))
    start = query_seq.find(model_seq)
    if start >= 0:
        return list(range(start, start + len(model_seq)))
    from Bio import Align
    aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-5,
                                    extend_gap_score=-0.5)
    aln = aligner.align(model_seq, query_seq)[0]
    mapping: list[int | None] = [None] * len(model_seq)
    for (m0, m1), (q0, q1) in zip(*aln.aligned):
        for k in range(m1 - m0):
            mapping[m0 + k] = q0 + k
    return mapping


def column_map(positions: list[tuple[str, str] | None], chain_order: list[str], feats: dict) -> np.ndarray:
    """MSA column for each scored position.

    `positions` lists the scored residues in model order as (chain_id, one-letter residue), or None for a
    linker position. `chain_order` gives the model chain matching each MSA chain (same order as the FASTA).
    Returns an int array with -1 where the position has no MSA column.
    """
    if len(chain_order) != len(feats["chains"]):
        raise ValueError(f"MSA has {len(feats['chains'])} chains, the model {len(chain_order)}")
    out = np.full(len(positions), -1, dtype=np.int64)
    for chain_id, chain in zip(chain_order, feats["chains"]):
        idx = [i for i, p in enumerate(positions) if p is not None and p[0] == chain_id]
        model_seq = "".join(positions[i][1] for i in idx)
        mapping = align_positions(model_seq, chain["sequence"])
        matched = sum(m is not None and chain["sequence"][m] == model_seq[k] for k, m in enumerate(mapping))
        if matched < 0.9 * len(model_seq):
            raise ValueError(f"chain {chain_id}: only {matched}/{len(model_seq)} residues match the MSA query")
        for i, m in zip(idx, mapping):
            if m is not None:
                out[i] = chain["columns"][m]
    return out


def _deletion_transform(d: np.ndarray) -> np.ndarray:
    return np.arctan(d / 3.0) * (2.0 / np.pi)


def build_msa_inputs(feats: dict, cols: np.ndarray, n_clusters: int, n_extra: int, multimer: bool) -> dict:
    """ColabDesign-shaped MSA arrays for the scored positions.

    Returns msa_feat (C, L, 49) for the cluster rows and extra_* arrays (E, L). For AF-multimer the extra
    deletion values stay raw counts (its model transforms them); for AF2-ptm they are pre-transformed.
    """
    msa = np.where(cols[None, :] >= 0, feats["msa"][:, np.clip(cols, 0, None)], GAP)
    dele = np.where(cols[None, :] >= 0, feats["deletion"][:, np.clip(cols, 0, None)], 0.0).astype(np.float32)
    # keep rows that are not gaps at every scored position (the query always stays)
    keep = np.concatenate([[True], (msa[1:] != GAP).any(1)])
    msa, dele = msa[keep], dele[keep]
    n_c = max(1, min(n_clusters, len(msa)))
    c_msa, c_del = msa[:n_c], dele[:n_c]
    e_msa, e_del = msa[n_c:n_c + n_extra], dele[n_c:n_c + n_extra]
    L = msa.shape[1]

    one_hot = lambda a: np.eye(23, dtype=np.float32)[a]  # noqa: E731
    c_1hot = one_hot(c_msa)
    profile, del_sum, count = c_1hot.copy(), c_del.copy(), np.ones((n_c, 1), np.float32)
    if len(e_msa):
        e_1hot = one_hot(e_msa)
        weights = np.ones(23, np.float32)
        weights[GAP] = 0.0  # gaps do not count as agreement (AlphaFold's gap_agreement_weight)
        agree = np.einsum("elk,clk->ec", e_1hot * weights, c_1hot)
        assign = agree.argmax(1)
        np.add.at(profile, assign, e_1hot)
        np.add.at(del_sum, assign, e_del)
        np.add.at(count, assign, 1.0)
    profile /= count[:, :, None]
    del_mean = del_sum / count
    has_del = np.clip(c_del, 0.0, 1.0)[..., None]
    msa_feat = np.concatenate([c_1hot, has_del, _deletion_transform(c_del)[..., None], profile,
                               _deletion_transform(del_mean)[..., None]], axis=-1).astype(np.float32)
    E = max(1, len(e_msa))
    extra_msa = np.full((E, L), GAP, np.int32)
    extra_del = np.zeros((E, L), np.float32)
    extra_mask = np.zeros((E, L), np.float32)
    if len(e_msa):
        extra_msa[:len(e_msa)], extra_del[:len(e_msa)], extra_mask[:len(e_msa)] = e_msa, e_del, 1.0
    return {
        "msa_feat": msa_feat,
        "msa_mask": np.ones((n_c, L), np.float32),
        "msa_row_mask": np.ones(n_c, np.float32),
        "extra_msa": extra_msa,
        "extra_has_deletion": np.clip(extra_del, 0.0, 1.0),
        "extra_deletion_value": extra_del if multimer else _deletion_transform(extra_del),
        "extra_msa_mask": extra_mask,
        "extra_msa_row_mask": (extra_mask.sum(1) > 0).astype(np.float32),
        "n_rows": int(len(msa)),
    }
