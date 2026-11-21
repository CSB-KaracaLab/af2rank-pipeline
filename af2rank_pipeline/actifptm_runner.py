#!/usr/bin/env python3
"""AF2Rank inference and actifpTM metrics for cleaned protein-complex models.

Uses installed ColabDesign and AlphaFold parameters supplied at runtime. Writes confidence
metrics, scored structures and PAE data. Raw debug arrays are optional.
"""

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
from copy import deepcopy
from pathlib import Path

import numpy as np


ACTIFPTM_VARIANT = "both_probability_weighted_and_binary_contacts"


def positive_iterations(value):
    try:
        iterations = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("iterations must be an integer of at least 1") from exc
    if iterations < 1:
        raise argparse.ArgumentTypeError("iterations must be at least 1")
    return iterations


def parse_args():
    parser = argparse.ArgumentParser(description="Run AF2Rank and ActifPTM cleanly.")
    parser.add_argument("--decoy_dir", required=True,
                        help="Directory containing input/decoy PDB files.")
    parser.add_argument("--chain", default=None, help="Comma-separated chain IDs, for example A,B,C.")
    parser.add_argument("--model_mode", default="alphafold", choices=["alphafold", "alphafold-multimer"])
    parser.add_argument("--model_num", type=int, default=2)
    parser.add_argument("--multimer_version", default="v3", choices=sorted(MULTIMER_SUFFIX),
                        help="AF2-multimer weights: v1 (AF2.1), v2 (AF2.2) or v3 (AF2.3).")
    parser.add_argument("--recycle", type=int, default=1)
    parser.add_argument("--iterations", type=positive_iterations, default=1)
    parser.add_argument("--params", required=True, help="AlphaFold params dir or parent data dir.")
    parser.add_argument("--colabdesign", default=None, help="Optional ColabDesign checkout override.")
    parser.add_argument("--tm", default=None, help="Optional TMscore executable; if absent, TM-score metrics are skipped.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where final CSV/JSON/PDB outputs are written.")
    parser.add_argument("--pae_json_dir", default=None,
                        help="Directory where per-model PAE JSON files are written.")
    parser.add_argument("--no_pae_json_gzip", action="store_false", dest="pae_json_gzip",
                        help="Write PAE JSON without gzip compression.")
    parser.add_argument("--target", default="target")
    parser.add_argument("--mask_sequence", action="store_true")
    parser.add_argument("--mask_sidechains", action="store_true")
    parser.add_argument("--mask_interchain", action="store_true")
    parser.add_argument("--keep_raw", action="store_true", help="Save raw PAE/distogram NPZ files for debugging.")
    parser.add_argument("--linker_length", type=int, default=0,
                        help="Join chains with this many glycines before multimer inference (0 = off).")
    parser.set_defaults(pae_json_gzip=True)
    args = parser.parse_args()
    if args.linker_length < 0:
        parser.error("--linker_length must be 0 or more")
    if args.linker_length and args.model_mode != "alphafold-multimer":
        parser.error("--linker_length applies only to --model_mode alphafold-multimer")
    return args


def setup_imports(colabdesign_dir):
    if colabdesign_dir:
        checkout = Path(colabdesign_dir).expanduser().resolve()
        if not (checkout / "colabdesign").is_dir():
            raise FileNotFoundError(f"ColabDesign checkout not found: {checkout}")
        if str(checkout) not in sys.path:
            sys.path.insert(0, str(checkout))
    # ColabDesign replaces XLA_FLAGS on import; keep the caller's flags too
    # (for example --xla_gpu_force_compilation_parallelism=1 on H100 GPUs).
    caller_flags = os.environ.get("XLA_FLAGS", "")
    try:
        import colabdesign  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "ColabDesign is unavailable. Install the inference extra with "
            "'uv sync --locked --all-extras', or pass --colabdesign to a checkout."
        ) from exc
    os.environ["XLA_FLAGS"] = merge_xla_flags(os.environ.get("XLA_FLAGS", ""), caller_flags)


def merge_xla_flags(*flag_strings):
    """Combine XLA flag strings; for a flag given more than once, the last value wins."""
    merged = {}
    for flags in flag_strings:
        for flag in flags.split():
            merged[flag.split("=", 1)[0]] = flag
    return " ".join(merged.values())


def get_data_dir(params_path):
    """ColabDesign expects data_dir/params/*.npz; users often pass the params dir itself."""
    p = Path(params_path).resolve()
    if p.name == "params":
        return str(p.parent)
    return str(p)


# Weight file suffix per AF2-multimer release: AF2.1, AF2.2, AF2.3.
MULTIMER_SUFFIX = {"v1": "multimer", "v2": "multimer_v2", "v3": "multimer_v3"}


def model_name_from_args(model_mode, model_num, multimer_version="v3"):
    if model_mode == "alphafold":
        return f"model_{model_num}_ptm"
    return f"model_{model_num}_{MULTIMER_SUFFIX[multimer_version]}"


def validate_params(params_path, model_mode, model_num, multimer_version="v3"):
    """Check that ColabDesign can read the requested params directly from shared storage."""
    params_path = Path(params_path).resolve()
    params_dir = params_path if params_path.name == "params" else params_path / "params"
    if not params_dir.exists():
        raise FileNotFoundError(f"Parameters directory not found: {params_dir}")

    if model_mode == "alphafold":
        expected = params_dir / f"params_model_{model_num}_ptm.npz"
    else:
        expected = params_dir / f"params_{model_name_from_args(model_mode, model_num, multimer_version)}.npz"

    if not expected.exists():
        raise FileNotFoundError(f"Expected model parameters not found: {expected}")


def parse_chain_labels(chain):
    if chain is None or chain.strip() == "":
        return None
    return [x.strip() for x in chain.split(",") if x.strip()]


def get_chain_lengths_from_pdb(pdb_path, chain_labels=None):
    """Return chain lengths from ATOM records, using unique residues per chain."""
    chain_order = []
    residues_by_chain = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip() or "_"
            if chain_labels is not None and chain_id not in chain_labels:
                continue
            res_key = (line[22:26].strip(), line[26].strip())
            if chain_id not in residues_by_chain:
                residues_by_chain[chain_id] = set()
                chain_order.append(chain_id)
            residues_by_chain[chain_id].add(res_key)

    if chain_labels is not None:
        chain_order = [c for c in chain_labels if c in residues_by_chain]

    return [len(residues_by_chain[c]) for c in chain_order], chain_order


def build_asym_id_from_lengths(lengths):
    if not lengths:
        return None
    return np.concatenate([np.full(length, i, dtype=np.int32) for i, length in enumerate(lengths)])


def get_asym_id(af_model, pdb_path, chain):
    inputs = af_model.model._inputs
    batch = inputs.get("batch", {})
    # Chain IDs must follow prepared residues rather than raw PDB counts.
    for features in (inputs, batch):
        if "asym_id" in features:
            asym_id = np.array(features["asym_id"])
            return asym_id[0] if asym_id.ndim > 1 else asym_id

    if "res_mask_per_chain" in af_model.model._inputs:
        masks = np.array(af_model.model._inputs["res_mask_per_chain"])
        asym_id = np.zeros(masks.shape[1], dtype=np.int32)
        for i in range(masks.shape[0]):
            asym_id[masks[i].astype(bool)] = i
        return asym_id

    labels = parse_chain_labels(chain)
    lengths, _ = get_chain_lengths_from_pdb(pdb_path, labels)
    asym_id = build_asym_id_from_lengths(lengths)
    if asym_id is not None:
        return asym_id

    L = len(batch["aatype"])
    return np.zeros(L, dtype=np.int32)


LINKED_CHAIN = "A"


def read_pdb_residues(pdb_path, chain_labels=None):
    """Return ATOM residues in file order as (chain, resseq, atom lines)."""
    residues = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21]
            if chain_labels is not None and chain_id not in chain_labels:
                continue
            if line[26].strip():
                raise ValueError(f"{Path(pdb_path).name}: insertion codes cannot be linked into one chain")
            resseq = int(line[22:26])
            if not residues or residues[-1][:2] != (chain_id, resseq):
                residues.append((chain_id, resseq, []))
            residues[-1][2].append(line.rstrip("\n"))
    if chain_labels is not None:
        order = {label: i for i, label in enumerate(chain_labels)}
        residues.sort(key=lambda residue: order[residue[0]])
    return residues


def _atom_coord(lines, atom_name):
    for line in lines:
        if line[12:16].strip() == atom_name:
            return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


def _atom_line(serial, atom_name, resname, chain_id, resseq, xyz, element):
    return (f"ATOM  {serial:5d}  {atom_name:<3s} {resname:>3s} {chain_id}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           {element:>1s}")


def write_linked_pdb(pdb_path, chain_labels, linker_length, out_path):
    """Join the chains into one chain separated by poly-glycine linkers.

    AlphaFold-multimer was trained with per-chain templates only. A single
    linked chain lets the whole complex act as one intra-chain template. Each
    chain keeps its internal residue-number gaps; linker residues get
    placeholder backbone atoms whose template information is removed before
    inference. Returns one (chain, resseq) entry per linked position, with
    None for linker positions.
    """
    residues = read_pdb_residues(pdb_path, chain_labels)
    chains = list(dict.fromkeys(chain_id for chain_id, _, _ in residues))
    position_map = []
    out_lines = []
    serial = 1
    next_number = 1
    previous_residues = None
    for chain_id in chains:
        chain_residues = [residue for residue in residues if residue[0] == chain_id]
        if previous_residues is not None:
            previous_ca = _atom_coord(previous_residues[-1][2], "CA")
            next_ca = _atom_coord(chain_residues[0][2], "CA")
            if previous_ca is None or next_ca is None:
                raise ValueError(f"{Path(pdb_path).name}: chain ends need CA atoms to place the linker")
            for k in range(1, linker_length + 1):
                ca = previous_ca + (next_ca - previous_ca) * k / (linker_length + 1)
                backbone = (("N", ca + (-1.2, 0.5, 0.0), "N"), ("CA", ca, "C"),
                            ("C", ca + (1.2, 0.5, 0.0), "C"), ("O", ca + (1.2, 1.7, 0.0), "O"))
                for atom_name, xyz, element in backbone:
                    out_lines.append(_atom_line(serial, atom_name, "GLY", LINKED_CHAIN, next_number, xyz, element))
                    serial += 1
                position_map.append(None)
                next_number += 1
        shift = next_number - chain_residues[0][1]
        for _, resseq, lines in chain_residues:
            for line in lines:
                out_lines.append(f"ATOM  {serial:5d}{line[11:21]}{LINKED_CHAIN}{resseq + shift:4d}{line[26:]}")
                serial += 1
            position_map.append((chain_id, resseq))
        next_number = chain_residues[-1][1] + shift + 1
        previous_residues = chain_residues
    if next_number > 10000:
        raise ValueError(f"{Path(pdb_path).name}: linked chain is too long for PDB residue numbering")
    Path(out_path).write_text("\n".join(out_lines) + "\nTER\nEND\n")
    return position_map


def unlink_pdb(linked_pdb, position_map, out_path):
    """Drop linker residues from a ColabDesign PDB and restore chains and numbering.

    Residues are matched by their order in the file, so this works whether
    ColabDesign keeps the linked numbering or renumbers the chain.
    """
    out_lines = []
    serial = 1
    previous_chain = None
    last_line = None
    position = -1
    residue_key = None
    with open(linked_pdb, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21:27] != residue_key:
                residue_key = line[21:27]
                position += 1
            if position >= len(position_map):
                raise ValueError(f"{Path(linked_pdb).name} has more residues than the linked input")
            original = position_map[position]
            if original is None:
                continue
            chain_id, resseq = original
            if previous_chain is not None and chain_id != previous_chain:
                out_lines.append(f"TER   {serial:5d}      {last_line[17:20]} {previous_chain}{last_line[22:26]}")
                serial += 1
            last_line = f"ATOM  {serial:5d}{line[11:21]}{chain_id}{resseq:4d}{line[26:].rstrip()}"
            out_lines.append(last_line)
            serial += 1
            previous_chain = chain_id
    if position + 1 != len(position_map):
        raise ValueError(
            f"{Path(linked_pdb).name} has {position + 1} residues, the linked input has {len(position_map)}")
    out_lines.append(f"TER   {serial:5d}      {last_line[17:20]} {previous_chain}{last_line[22:26]}")
    Path(out_path).write_text("\n".join(out_lines) + "\nEND\n")


def linked_asym_id(position_map):
    """Chain index per non-linker position, in linked order."""
    chains = list(dict.fromkeys(entry[0] for entry in position_map if entry is not None))
    return np.array([chains.index(entry[0]) for entry in position_map if entry is not None], dtype=np.int32)


def drop_linker_outputs(outputs, keep):
    """Remove linker rows and columns from the PAE and distogram outputs."""
    trimmed = {}
    for name in ("predicted_aligned_error", "distogram"):
        head = dict(outputs[name])
        head["logits"] = np.asarray(head["logits"])[keep][:, keep]
        trimmed[name] = head
    return trimmed


def alphafold_tm_score(logits, breaks, asym_id=None):
    """AlphaFold's predicted TM-score, as ColabDesign logs it: pTM, or ipTM with asym_id."""
    logits = np.asarray(logits, dtype=float)
    breaks = np.asarray(breaks, dtype=float)
    step = breaks[1] - breaks[0]
    bin_centers = np.append(breaks + step / 2, breaks[-1] + 1.5 * step)
    d0 = 1.24 * (max(logits.shape[0], 19) - 15) ** (1.0 / 3) - 1.8
    tm_per_pair = (softmax(logits, axis=-1) / (1 + np.square(bin_centers / d0))).sum(-1)
    if asym_id is None:
        pair_mask = np.ones(tm_per_pair.shape)
    else:
        asym_id = np.asarray(asym_id)
        pair_mask = (asym_id[:, None] != asym_id[None, :]).astype(float)
    weights = pair_mask / (1e-8 + pair_mask.sum(-1, keepdims=True))
    return float(np.max((tm_per_pair * pair_mask * weights).sum(-1)))


def linked_confidence(aux, outputs, keep, asym_id):
    """Recompute pLDDT, pTM, ipTM and PAE without the linker residues.

    Values follow ColabDesign's logged definitions so linked and unlinked runs
    compare directly: mean pLDDT (0-1), AlphaFold pTM/ipTM, and the mean PAE
    over within-chain residue pairs. `outputs` must already be trimmed with
    drop_linker_outputs.
    """
    pae_head = outputs["predicted_aligned_error"]
    logits, breaks = pae_head["logits"], pae_head["breaks"]
    pae = np.asarray(aux["pae"])[keep][:, keep]
    asym_id = np.asarray(asym_id)
    same_chain = asym_id[:, None] == asym_id[None, :]
    return {
        "plddt": float(np.mean(np.asarray(aux["plddt"])[keep])),
        "ptm": alphafold_tm_score(logits, breaks),
        "i_ptm": alphafold_tm_score(logits, breaks, asym_id) if len(np.unique(asym_id)) > 1 else np.nan,
        "pae": pae,
        "pae_mean": float(pae[same_chain].mean()),
    }


def chain_ranges_from_asym_id(asym_id, chain_labels=None):
    asym_id = np.array(asym_id)
    if asym_id.ndim > 1:
        asym_id = asym_id[0]
    unique = list(np.unique(asym_id))
    if chain_labels is None:
        chain_labels = [chr(ord("A") + i) for i in range(len(unique))]
    ranges = {}
    for i, chain_value in enumerate(unique):
        positions = np.where(asym_id == chain_value)[0]
        if len(positions) == 0:
            continue
        label = chain_labels[i] if i < len(chain_labels) else f"chain{i + 1}"
        ranges[label] = (int(positions[0]), int(positions[-1]))
    return ranges


def calculate_bin_centers(breaks, use_jnp=False):
    import jax.numpy as jnp

    xp = jnp if use_jnp else np
    breaks = xp.array(breaks)
    return 0.5 * (breaks[:-1] + breaks[1:])


def softmax(x, axis=-1):
    import scipy.special

    return scipy.special.softmax(x, axis=axis)


def predicted_tm_score_modified(logits, breaks, residue_weights=None, asym_id=None,
                                pair_residue_weights=None, use_jnp=False):
    import jax.numpy as jnp

    xp = jnp if use_jnp else np
    logits = xp.array(logits)
    breaks = xp.array(breaks)

    if residue_weights is None:
        residue_weights = xp.ones(logits.shape[0])
    else:
        residue_weights = xp.array(residue_weights)

    bin_centers = calculate_bin_centers(breaks, use_jnp=use_jnp)
    if bin_centers.shape[0] != logits.shape[-1]:
        recreated_breaks = xp.linspace(0, 31, logits.shape[-1] + 1)
        bin_centers = calculate_bin_centers(recreated_breaks, use_jnp=use_jnp)

    num_res = residue_weights.shape[0]
    clipped_num_res = xp.maximum(num_res, 19)
    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8
    tm_per_bin = 1.0 / (1 + xp.square(bin_centers) / xp.square(d0))

    probs = jnp.asarray(softmax(np.asarray(logits), axis=-1)) if use_jnp else softmax(logits, axis=-1)
    predicted_tm_term = (probs * tm_per_bin).sum(-1)

    if asym_id is None:
        pair_mask = xp.ones(predicted_tm_term.shape, dtype=bool)
    else:
        asym_id = xp.array(asym_id)
        pair_mask = asym_id[:, None] != asym_id[None, :]

    predicted_tm_term *= pair_mask

    if pair_residue_weights is None:
        pair_residue_weights = pair_mask * (residue_weights[None, :] * residue_weights[:, None])
    else:
        pair_residue_weights = xp.array(pair_residue_weights)

    normed_residue_mask = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1, keepdims=True))
    per_alignment = (predicted_tm_term * normed_residue_mask).sum(-1)
    return per_alignment * residue_weights


def distogram_bins(distogram):
    logits = np.array(distogram["logits"])
    bin_edges = np.array(distogram["bin_edges"])

    if logits.shape[-1] == 64:
        if bin_edges.shape[0] == 63:
            return np.append(0, bin_edges)
        return np.append(0, np.linspace(2.3125, 21.6875, 63))
    if logits.shape[-1] == 39:
        return np.linspace(3.25, 50.75, 39) + 1.25
    if bin_edges.shape[0] == logits.shape[-1]:
        return bin_edges

    return np.linspace(0, 31, logits.shape[-1])


def contact_map(distogram, dist=8.0):
    logits = np.array(distogram["logits"])
    bins = distogram_bins(distogram)
    return (softmax(logits, axis=-1) * (bins < dist)).sum(-1)


def compute_pair_metrics(outputs, asym_id, chain_labels=None):
    import jax.numpy as jnp

    asym_id = np.array(asym_id)
    if asym_id.ndim > 1:
        asym_id = asym_id[0]

    chain_ranges = chain_ranges_from_asym_id(asym_id, chain_labels)
    distogram = outputs["distogram"]
    pae = outputs["predicted_aligned_error"]
    cmap = contact_map(distogram)

    metrics = {
        "pairwise_actifptm": {},
        "pairwise_actifptm_binary": {},
        "pairwise_iptm": {},
        "per_chain_ptm": {},
        "actifptm": None,
        "actifptm_prob": None,
        "actifptm_binary": None,
    }

    full_length = len(asym_id)
    breaks = pae["breaks"]
    logits = pae["logits"]

    for label, (start, end) in chain_ranges.items():
        sliced_logits = jnp.array(logits)[start:end + 1, start:end + 1, :]
        residue_weights = jnp.ones(sliced_logits.shape[0], dtype=float)
        score = predicted_tm_score_modified(
            sliced_logits, breaks, residue_weights=residue_weights, use_jnp=True
        )
        metrics["per_chain_ptm"][label] = round(float(score.max()), 3)

    labels = list(chain_ranges.keys())
    full_prob_weights = np.zeros((full_length, full_length), dtype=float)
    full_binary_weights = np.zeros((full_length, full_length), dtype=float)

    for i, label_i in enumerate(labels):
        start_i, end_i = chain_ranges[label_i]
        for j, label_j in enumerate(labels):
            if j <= i:
                continue
            start_j, end_j = chain_ranges[label_j]
            key = f"{label_i}_{label_j}"

            seq_mask = np.zeros(full_length, dtype=float)
            seq_mask[start_i:end_i + 1] = 1.0
            seq_mask[start_j:end_j + 1] = 1.0

            pair_prob_weights = np.zeros((full_length, full_length), dtype=float)
            pair_prob_weights[start_i:end_i + 1, start_j:end_j + 1] = cmap[start_i:end_i + 1, start_j:end_j + 1]
            pair_prob_weights[start_j:end_j + 1, start_i:end_i + 1] = cmap[start_j:end_j + 1, start_i:end_i + 1]
            full_prob_weights += pair_prob_weights

            actif_prob = predicted_tm_score_modified(
                logits, breaks, residue_weights=seq_mask, asym_id=asym_id,
                pair_residue_weights=pair_prob_weights, use_jnp=True
            )
            metrics["pairwise_actifptm"][key] = round(float(actif_prob.max()), 3)

            contacts = np.where(cmap[start_i:end_i + 1, start_j:end_j + 1] >= 0.6)
            contact_seq_mask = np.zeros(full_length, dtype=float)
            if contacts[0].size > 0:
                global_i = contacts[0] + start_i
                global_j = contacts[1] + start_j
                contact_positions = np.unique(np.concatenate([global_i, global_j]))
                contact_seq_mask[contact_positions] = 1.0
                actif_binary = predicted_tm_score_modified(
                    logits, breaks, residue_weights=contact_seq_mask, asym_id=asym_id, use_jnp=True
                )
                metrics["pairwise_actifptm_binary"][key] = round(float(actif_binary.max()), 3)
                full_binary_weights += contact_seq_mask[None, :] * contact_seq_mask[:, None]
            else:
                metrics["pairwise_actifptm_binary"][key] = 0.0

            iptm = predicted_tm_score_modified(
                logits, breaks, residue_weights=seq_mask, asym_id=asym_id, use_jnp=True
            )
            metrics["pairwise_iptm"][key] = round(float(iptm.max()), 3)

    if len(labels) > 1:
        residue_weights = np.ones(full_length, dtype=float)
        actif_prob = predicted_tm_score_modified(
            logits, breaks, residue_weights=residue_weights, asym_id=asym_id,
            pair_residue_weights=full_prob_weights, use_jnp=True
        )
        pair_mask = asym_id[:, None] != asym_id[None, :]
        full_binary_weights *= pair_mask
        actif_binary = predicted_tm_score_modified(
            logits, breaks, residue_weights=residue_weights, asym_id=asym_id,
            pair_residue_weights=full_binary_weights, use_jnp=True
        )
        metrics["actifptm"] = round(float(actif_prob.max()), 3)
        metrics["actifptm_prob"] = metrics["actifptm"]
        metrics["actifptm_binary"] = round(float(actif_binary.max()), 3)
    else:
        metrics["actifptm"] = 0.0
        metrics["actifptm_prob"] = 0.0
        metrics["actifptm_binary"] = 0.0

    return metrics


def validate_tm_executable(tm_executable):
    if not tm_executable:
        return None
    path = Path(tm_executable).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"TMscore executable not found: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"TMscore executable is not executable: {path}")
    return str(path)


def tmscore(xyz_a, xyz_b, tm_executable):
    tm_executable = validate_tm_executable(tm_executable)
    if tm_executable is None:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_a = os.path.join(tmpdir, "a.pdb")
        pdb_b = os.path.join(tmpdir, "b.pdb")
        write_ca_pdb(xyz_a, pdb_a)
        write_ca_pdb(xyz_b, pdb_b)
        proc = subprocess.run(
            [tm_executable, pdb_a, pdb_b],
            check=False,
            capture_output=True,
            text=True,
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"TMscore exited with status {proc.returncode}: {proc.stderr.strip()[:200]}"
        )

    for line in proc.stdout.splitlines():
        if line.startswith("TM-score"):
            try:
                return float(line.split()[2])
            except (IndexError, ValueError):
                break
    raise ValueError("TMscore completed without a parseable TM-score value")


def write_ca_pdb(xyz, filename):
    with open(filename, "w", encoding="utf-8") as handle:
        for i, (x, y, z) in enumerate(np.array(xyz), start=1):
            handle.write(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        handle.write("TER\nEND\n")


class AF2RankRunner:
    def __init__(self, pdb, chain, model_name, data_dir):
        from colabdesign import mk_af_model
        from colabdesign.shared.utils import copy_dict

        self.copy_dict = copy_dict
        self.args = {
            "pdb": pdb,
            "chain": chain,
            "model_name": model_name,
            "use_multimer": "multimer" in model_name,
        }
        self.model = mk_af_model(
            protocol="fixbb",
            use_templates=True,
            use_multimer=self.args["use_multimer"],
            debug=True,
            model_names=[model_name],
            data_dir=data_dir,
        )
        self.model.prep_inputs(pdb, chain=chain)
        self.model.set_seq(mode="wildtype")
        self.wt_batch = copy_dict(self.model._inputs["batch"])
        self.wt = self.model._wt_aatype

    def set_pdb(self, pdb, chain):
        self.model.prep_inputs(pdb, chain=chain)
        self.model.set_seq(mode="wildtype")
        self.wt_batch = self.copy_dict(self.model._inputs["batch"])
        self.wt = self.model._wt_aatype

    def predict(self, pdb, chain, recycles, iterations, rm_seq, rm_sc, rm_ic, output_pdb=None,
                no_template=None):
        """Score one structure; `no_template` marks positions (e.g. a linker) with no template."""
        if iterations < 1:
            raise ValueError("iterations must be at least 1")
        if pdb is not None:
            self.set_pdb(pdb, chain)

        self.model._inputs["batch"]["aatype"] = self.wt
        self.model.set_opt(template=dict(rm_ic=rm_ic), num_recycles=recycles)
        self.model._inputs["rm_template"][:] = False
        self.model._inputs["rm_template_sc"][:] = rm_sc
        self.model._inputs["rm_template_seq"][:] = rm_seq
        if no_template is not None:
            no_template = np.asarray(no_template, dtype=bool)
            if no_template.shape[0] != len(self.wt):
                raise ValueError(
                    f"Linked input has {len(self.wt)} prepared residues, expected {no_template.shape[0]}")
            self.model._inputs["rm_template"][no_template] = True
            self.model._inputs["batch"]["all_atom_mask"][no_template] = 0

        initial_atoms = self.model._inputs["batch"]["all_atom_positions"].copy()
        for i in range(iterations):
            self.model.predict(models=self.args["model_name"], verbose=False)
            if i < iterations - 1:
                self.model._inputs["batch"]["all_atom_positions"] = self.model.aux["atom_positions"]
            else:
                self.model._inputs["batch"]["all_atom_positions"] = initial_atoms

        if output_pdb is not None:
            self.model.save_pdb(output_pdb)

        return self.model.aux


def flatten_metrics(row):
    out = {}
    for key, value in row.items():
        if key in {"pairwise_actifptm", "pairwise_actifptm_binary", "pairwise_iptm", "per_chain_ptm"}:
            continue
        out[key] = value

    for key, value in row.get("pairwise_actifptm", {}).items():
        out[f"actifptm_{key}"] = value
    for key, value in row.get("pairwise_actifptm_binary", {}).items():
        out[f"actifptm_binary_{key}"] = value
    for key, value in row.get("pairwise_iptm", {}).items():
        out[f"iptm_{key}"] = value
    for key, value in row.get("per_chain_ptm", {}).items():
        out[f"chain_ptm_{key}"] = value
    return out


def write_csv(rows, filename):
    flattened = [flatten_metrics(row) for row in rows]
    fieldnames = sorted({key for row in flattened for key in row.keys()})
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="N/A")
        writer.writeheader()
        for row in flattened:
            writer.writerow(row)


def maybe_save_raw_npz(out_dir, pdb_id, aux, asym_id):
    outputs = aux["debug"]["outputs"]
    raw_dir = Path(out_dir) / "raw_npz"
    raw_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        raw_dir / f"{pdb_id}_raw_outputs.npz",
        pae=np.array(aux.get("pae")),
        pae_logits=np.array(outputs["predicted_aligned_error"]["logits"]),
        pae_breaks=np.array(outputs["predicted_aligned_error"]["breaks"]),
        distogram_logits=np.array(outputs["distogram"]["logits"]),
        distogram_bins=np.array(outputs["distogram"]["bin_edges"]),
        asym_id=np.array(asym_id),
    )


def save_pae_json(pae_json_dir, pdb_id, aux, gzip_enabled=True):
    pae_json_dir = Path(pae_json_dir)
    pae_json_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".json.gz" if gzip_enabled else ".json"
    path = pae_json_dir / f"{pdb_id}{suffix}"
    tmp_path = path.with_name(f".{path.name}.tmp")
    payload = {"model_id": pdb_id, "pae": np.asarray(aux.get("pae"), dtype=float).tolist()}
    if gzip_enabled:
        with gzip.open(tmp_path, "wt") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.write("\n")
    else:
        with tmp_path.open("w") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.write("\n")
    tmp_path.replace(path)
    return path


def main():
    args = parse_args()
    args.tm = validate_tm_executable(args.tm)
    setup_imports(args.colabdesign)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pae_json_dir = Path(args.pae_json_dir) if args.pae_json_dir else out_dir / "pae_json"
    pdb_out_dir = out_dir / "scored_pdbs"
    pdb_out_dir.mkdir(exist_ok=True)
    linked_dir = out_dir / "linked_inputs"
    if args.linker_length:
        linked_dir.mkdir(exist_ok=True)

    version = args.multimer_version if args.model_mode == "alphafold-multimer" else "ptm"
    model_name = model_name_from_args(args.model_mode, args.model_num, args.multimer_version)
    validate_params(args.params, args.model_mode, args.model_num, args.multimer_version)
    data_dir = get_data_dir(args.params)
    chain_labels = parse_chain_labels(args.chain)

    pdb_files = sorted(Path(args.decoy_dir).glob("*.pdb"))
    if not pdb_files:
        raise SystemExit(f"No .pdb files found in {args.decoy_dir}")

    rows = []
    # One model per weight, reused for every structure: building one per structure reloads the
    # weights each time and, in long ensembles, keeps growing host memory until the process is killed.
    runner = None
    for pdb_path in pdb_files:
        pdb_id = pdb_path.stem
        print(f"[{args.target}] scoring {pdb_path.name} with {model_name}", flush=True)

        output_pdb = str(pdb_out_dir / f"{pdb_id}_scored.pdb")
        if args.linker_length:
            # Score one linked chain so the complex is a single intra-chain template.
            model_pdb = linked_dir / pdb_path.name
            position_map = write_linked_pdb(pdb_path, chain_labels, args.linker_length, model_pdb)
            model_chain = LINKED_CHAIN
            no_template = np.array([entry is None for entry in position_map])
            keep = np.flatnonzero(~no_template)
            raw_output_pdb = str(linked_dir / f"{pdb_id}_scored_linked.pdb")
        else:
            model_pdb, model_chain, no_template, keep = pdb_path, args.chain, None, None
            raw_output_pdb = output_pdb
        if runner is None:
            runner = AF2RankRunner(str(model_pdb), model_chain, model_name, data_dir)
        aux = runner.predict(
            str(model_pdb),
            model_chain,
            recycles=args.recycle,
            iterations=args.iterations,
            rm_seq=args.mask_sequence,
            rm_sc=args.mask_sidechains,
            rm_ic=args.mask_interchain,
            output_pdb=raw_output_pdb,
            no_template=no_template,
        )

        log = deepcopy(aux["log"])
        input_ca = runner.model._inputs["batch"]["all_atom_positions"][:, 1]
        output_ca = np.array(aux["atom_positions"][:, 1])
        native_ca = runner.wt_batch["all_atom_positions"][:, 1]
        if args.linker_length:
            unlink_pdb(raw_output_pdb, position_map, output_pdb)
            # Record which residues of the linked PDBs are linker glycines, for viewers.
            linker_residues = [resseq for (_, resseq, _), entry in zip(read_pdb_residues(model_pdb), position_map)
                               if entry is None]
            (linked_dir / f"{pdb_id}_linker.json").write_text(json.dumps(
                {"chain": LINKED_CHAIN, "linker_length": args.linker_length,
                 "linker_residues": linker_residues}) + "\n")
            asym_id = linked_asym_id(position_map)
            outputs = drop_linker_outputs(aux["debug"]["outputs"], keep)
            confidence = linked_confidence(aux, outputs, keep, asym_id)
            plddt, ptm, iptm = confidence["plddt"], confidence["ptm"], confidence["i_ptm"]
            pae_mean = confidence["pae_mean"]
            saved_aux = {"pae": confidence["pae"], "debug": {"outputs": outputs}}
            input_ca, output_ca, native_ca = input_ca[keep], output_ca[keep], native_ca[keep]
        else:
            asym_id = get_asym_id(runner, str(pdb_path), args.chain)
            outputs = aux["debug"]["outputs"]
            plddt = float(log.get("plddt", np.nan))
            ptm = float(log.get("ptm", np.nan))
            iptm = float(log.get("i_ptm", np.nan)) if "i_ptm" in log else np.nan
            pae_mean = float(31.0 * log["pae"]) if "pae" in log else np.nan
            saved_aux = aux
        actif_metrics = compute_pair_metrics(outputs, asym_id, chain_labels)
        tm_i = tmscore(native_ca, input_ca, args.tm)
        tm_o = tmscore(native_ca, output_ca, args.tm)
        tm_io = tmscore(input_ca, output_ca, args.tm)
        composite = ptm * plddt * tm_io if tm_io is not None and not np.isnan(ptm) else np.nan

        row = {
            "target": args.target,
            "id": pdb_path.name,
            "model_name": model_name,
            "model_mode": args.model_mode,
            "version": version,
            "model_num": args.model_num,
            "recycle": args.recycle,
            "iterations": args.iterations,
            "chain": args.chain or "",
            "rm_seq": bool(args.mask_sequence),
            "rm_sc": bool(args.mask_sidechains),
            "rm_ic": bool(args.mask_interchain),
            "linker_length": args.linker_length,
            "plddt": round(plddt, 4) if not np.isnan(plddt) else "N/A",
            "ptm": round(ptm, 4) if not np.isnan(ptm) else "N/A",
            "i_ptm": round(iptm, 4) if not np.isnan(iptm) else "N/A",
            "tm_i": round(tm_i, 4) if tm_i is not None else "N/A",
            "tm_o": round(tm_o, 4) if tm_o is not None else "N/A",
            "tm_io": round(tm_io, 4) if tm_io is not None else "N/A",
            "composite": round(float(composite), 4) if not np.isnan(composite) else "N/A",
            "pae": round(pae_mean, 4) if not np.isnan(pae_mean) else "N/A",
        }
        row.update(actif_metrics)
        rows.append(row)

        if args.keep_raw:
            maybe_save_raw_npz(out_dir, pdb_id, saved_aux, asym_id)
        save_pae_json(pae_json_dir, pdb_id, saved_aux, gzip_enabled=args.pae_json_gzip)

        # Release large arrays between structures.
        del aux, outputs

    csv_path = out_dir / f"{args.target}_af2rank_actifptm.csv"
    write_csv(rows, csv_path)

    metadata = {
        "target": args.target,
        "decoy_dir": args.decoy_dir,
        "chain": args.chain,
        "model_name": model_name,
        "params": args.params,
        "colabdesign": args.colabdesign,
        "num_pdbs": len(pdb_files),
        "csv": str(csv_path),
        "raw_npz_saved": bool(args.keep_raw),
        "iterations": args.iterations,
        "pae_json_dir": str(pae_json_dir),
        "pae_json_gzip": bool(args.pae_json_gzip),
        "scored_pdbs_saved": True,
        "actifptm_variant": ACTIFPTM_VARIANT,
        "linker_length": args.linker_length,
    }
    with open(out_dir / f"{args.target}_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Wrote {csv_path}", flush=True)


if __name__ == "__main__":
    main()
