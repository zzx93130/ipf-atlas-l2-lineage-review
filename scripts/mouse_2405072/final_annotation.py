"""Final annotated h5ad for mouse 2405072 (L1 lineages + user L2 annotation).

Usage:  python final_annotation.py
Reads, per lineage, the clustered h5ad of the round the user annotated (L2/<lineage>/ or L2/<lineage>/round<k>/) and
the cluster ids at the resolution in that round's cluster_table.csv; maps clusters -> labels below.

Writes  final/mouse_2405072_L2_final.h5ad          kept cells; X = raw counts, layers counts + lognorm;
                                                    obs = metadata, QC, L1 + L2 columns; obsm X_umap_L1 (whole-lung
                                                    L1 UMAP), X_umap_lineage (UMAP of the annotated round, per lineage)
        final/mouse_2405072_L2_labels_all_cells.parquet   every cell incl. QC / doublet / L1 / L2 removals + reason
        final/L2_decisions_final.csv                one row per lineage x cluster
        final/L2_cell_type_counts_by_sample.csv, final/figures/*.png
"""
import os
import sys
from datetime import date
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_mouse_final")
sys.modules["tensorflow"] = None

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

ROOT = Path(os.environ.get("MOUSE_ROOT", "/Users/joe/SSD2_1T/IPF_ATLAS/data/processed/mouse_2405072"))
OUT = ROOT / "final"
STATES = {"", "Proliferating", "IFN response", "Stress response", "Metallothionein-high", "Mitochondrial-high"}
# L2 Scrublet filter applied in l2_scanpy_lineage.py before clustering (keeps scrublet_score < threshold);
# must match the thresholds used there. Cells it dropped are recorded as removed_L2 with their score.
SCRUBLET_MAX = {"endothelial": 0.07}
SCRUBLET_MAX_DEFAULT = 0.2

# ------------------------------------------------------------------ user annotation (2026-11-02)
# cluster: (as typed by the user, cell_type, L2_identity, L2_state)   or   "remove"
# cell_type = user's wording completed to full names; L2_identity = cell_type without the state;
# L2_state = closed state list (heat-shock / stressed -> Stress response).
R = "remove"
ANNOTATION = {
    "myeloid": {"round": 1, "res": None, "clusters": {
        "18": ("MK", "Megakaryocyte", "Megakaryocyte", ""),
        "14": ("pDC", "pDC", "pDC", ""),
        "9": ("IFN response Neutrophil", "IFN-response neutrophil", "Neutrophil", "IFN response"),
        "17": ("Neutrophil", "Neutrophil", "Neutrophil", ""),
        "8": ("Neutrophil", "Neutrophil", "Neutrophil", ""),
        "16": ("AM", "Alveolar macrophage", "Alveolar macrophage", ""),
        "5": ("AM", "Alveolar macrophage", "Alveolar macrophage", ""),
        "10": ("proliferating myeloid", "Proliferating myeloid", "Myeloid", "Proliferating"),
        "1": ("Spp1 macro", "Spp1+ macrophage", "Spp1+ macrophage", ""),
        "4": ("activated macro", "Activated macrophage", "Activated macrophage", ""),
        "3": ("mregDC", "mregDC", "mregDC", ""),
        "12": ("mregDC", "mregDC", "mregDC", ""),
        "0": ("IM", "Interstitial macrophage", "Interstitial macrophage", ""),
        "2": ("IM", "Interstitial macrophage", "Interstitial macrophage", ""),
        "11": ("cDC1", "cDC1", "cDC1", ""),
        "6": ("cDC1", "cDC1", "cDC1", ""),
        "13": ("Non-classic Mono", "Non-classical monocyte", "Non-classical monocyte", ""),
        "15": ("classic Mono", "Classical monocyte", "Classical monocyte", ""),
        "7": ("cDC2", "cDC2", "cDC2", ""),
    }},
    "epithelial": {"round": 1, "res": None, "clusters": {
        "13": ("proliferating epi", "Proliferating epithelial", "Epithelial", "Proliferating"),
        "0": ("AT1", "AT1", "AT1", ""),
        "4": ("AT1", "AT1", "AT1", ""),
        "5": ("AT1", "AT1", "AT1", ""),
        "12": ("VIM+AT1", "Vim+ AT1", "Vim+ AT1", ""),
        "6": ("AT2", "AT2", "AT2", ""),
        "7": ("AT2", "AT2", "AT2", ""),
        "1": ("AT2", "AT2", "AT2", ""),
        "14": R,
        "8": ("Goblet", "Goblet", "Goblet", ""),
        "11": ("Basal", "Basal", "Basal", ""),
        "10": ("Multiciliated", "Multiciliated", "Multiciliated", ""),
        "9": ("Multiciliated", "Multiciliated", "Multiciliated", ""),
        "2": ("Aberrant basaloid", "Aberrant basaloid", "Aberrant basaloid", ""),
        "3": ("Club", "Club", "Club", ""),
    }},
    "endothelial": {"round": 1, "res": None, "clusters": {
        "8": ("System Venous", "Systemic venous EC", "Systemic venous EC", ""),
        "10": ("Arterial", "Arterial EC", "Arterial EC", ""),
        "9": ("Venous", "Venous EC", "Venous EC", ""),
        "5": ("proliferating Endo", "Proliferating endothelial", "Endothelial", "Proliferating"),
        "11": R,
        "2": ("aCap", "aCap", "aCap", ""),
        "3": ("aCap", "aCap", "aCap", ""),
        "7": ("Sox_high Endo", "Sox-high endothelial", "Sox-high endothelial", ""),
        "1": ("gCap", "gCap", "gCap", ""),
        "6": ("gCap", "gCap", "gCap", ""),
        "13": R,
        "4": ("Stressed Endo", "Stressed endothelial", "Endothelial", "Stress response"),
        "12": ("Heat-shock Endo", "Heat-shock endothelial", "Endothelial", "Stress response"),
        "0": ("Low-feature Endo", "Low-feature endothelial", "Low-feature endothelial", ""),
    }},
    "stromal": {"round": 2, "res": None, "clusters": {
        "15": ("IFN-respones stromal", "IFN-response stromal", "Stromal", "IFN response"),
        "5": ("proliferating stromal", "Proliferating stromal", "Stromal", "Proliferating"),
        "14": ("mesothelial", "Mesothelial", "Mesothelial", ""),
        "3": ("Peribronchial Fibro", "Peribronchial fibroblast", "Peribronchial fibroblast", ""),
        "2": ("SMC", "SMC", "SMC", ""),
        "11": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "9": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "13": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "6": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "12": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "0": ("Alveolar Fibro", "Alveolar fibroblast", "Alveolar fibroblast", ""),
        "1": ("Pericyte", "Pericyte", "Pericyte", ""),
        "10": ("CD36+Fibro", "Cd36+ fibroblast", "Cd36+ fibroblast", ""),
        "7": ("Adventitial Fibro", "Adventitial fibroblast", "Adventitial fibroblast", ""),
        "4": ("Myofibroblast", "Myofibroblast", "Myofibroblast", ""),
        "8": ("Myofibroblast", "Myofibroblast", "Myofibroblast", ""),
    }},
    "tnk": {"round": 2, "res": None, "clusters": {
        "14": R,
        "9": ("IFN-response T/NK", "IFN-response T/NK", "T/NK", "IFN response"),
        "0": ("proliferating", "Proliferating T/NK", "T/NK", "Proliferating"),
        "13": ("ILC2", "ILC2", "ILC2", ""),
        "10": ("CD4", "CD4 T", "CD4 T", ""),
        "2": ("activated CD8", "Activated CD8 T", "Activated CD8 T", ""),
        "3": ("memory CD4", "Memory CD4 T", "Memory CD4 T", ""),
        "1": ("heat-shock CD4", "Heat-shock CD4 T", "CD4 T", "Stress response"),
        "8": ("DP", "DP T (CD4+CD8+)", "DP T (CD4+CD8+)", ""),
        "12": ("ILC3", "ILC3", "ILC3", ""),
        "6": ("Treg", "Treg", "Treg", ""),
        "11": ("NK", "NK", "NK", ""),
        "4": ("activated CD4", "Activated CD4 T", "Activated CD4 T", ""),
        "7": ("Naive CD8", "Naive CD8 T", "Naive CD8 T", ""),
        "5": ("Naive CD4", "Naive CD4 T", "Naive CD4 T", ""),
    }},
    "bplasma": {"round": 2, "res": None, "clusters": {
        "11": ("Plasma", "Plasma cell", "Plasma cell", ""),
        **{c: ("stressed B", "Stressed B cell", "B cell", "Stress response") for c in ["3", "4", "9", "0"]},
        **{c: ("B cell", "B cell", "B cell", "") for c in ["1", "2", "5", "6", "7", "8", "10", "12"]},
    }},
}
LINEAGE_NAME = {"epithelial": "Epithelial", "endothelial": "Endothelial", "stromal": "Stromal",
                "myeloid": "Myeloid", "tnk": "T/NK", "bplasma": "B/plasma"}


# ------------------------------------------------------------------ helpers
def round_dir(lin, k):
    return ROOT / "L2" / lin if k == 1 else ROOT / "L2" / lin / f"round{k}"


def round_file(lin, k):
    return round_dir(lin, k) / (f"mouse_2405072_{lin}_L2_clustered.h5ad" if k == 1
                                else f"mouse_2405072_{lin}_L2_round{k}_clustered.h5ad")


def reviewed_res(lin, k):
    """resolution the user reviewed = index name of that round's cluster table (leiden_L2_<res>)."""
    tag = lin if k == 1 else f"{lin}_round{k}"
    f = round_dir(lin, k) / f"{tag}_cluster_table.csv"
    key = pd.read_csv(f, nrows=0).columns[0]
    assert key.startswith("leiden_L2_"), f"unexpected cluster table header {key} in {f}"
    return key, pd.read_csv(f, dtype={key: str}).set_index(key)["n_cells"]


def say(msg):
    print(msg, flush=True)


# ------------------------------------------------------------------ build
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(exist_ok=True)
sc.settings.figdir = OUT / "figures"
today = str(date.today())

L1 = pd.read_parquet(ROOT / "L1" / "annotation" / "L1_labels.parquet").set_index("cell_id")
labels = L1[["sample_id", "condition", "L1_status", "L1_lineage", "L1_call", "L1_state", "L1_removal_reason"]].copy()
for c in ["L2_round", "L2_cluster_key", "L2_cluster", "L2_cluster_uid", "user_label_raw", "cell_type", "L2_identity",
          "L2_state", "L2_removal_reason"]:
    labels[c] = ""
labels["final_status"] = np.where(labels["L1_status"] == "kept", "",
                                  labels["L1_status"].replace({"removed": "removed_L1"}))

parts, decisions = [], []
for lin, cfg in ANNOTATION.items():
    k = cfg["round"]
    key, sizes = reviewed_res(lin, k)
    if cfg["res"] is not None and key != f"leiden_L2_{cfg['res']}":
        key = f"leiden_L2_{cfg['res']}"
        sizes = None
    a = sc.read_h5ad(round_file(lin, k))
    cl = a.obs[key].astype(str)
    if sizes is not None:  # the h5ad must be the one the cluster table was made from
        vc = cl.value_counts()
        assert set(vc.index) == set(sizes.index) and (vc.reindex(sizes.index).values == sizes.astype(int).values).all(), \
            f"{lin}: cluster sizes in {round_file(lin, k).name} differ from the round's cluster_table.csv"
    have, given = set(cl), set(cfg["clusters"])
    if have != given:
        raise SystemExit(f"{lin} round {k} {key}: clusters in data {sorted(have, key=int)}; "
                         f"missing from ANNOTATION {sorted(have - given, key=int)}; "
                         f"not in data {sorted(given - have, key=int)}")
    for c, v in cfg["clusters"].items():
        assert v == R or (len(v) == 4 and v[3] in STATES), f"{lin} cluster {c}: bad entry {v}"
    say(f"== {lin}: round {k}, {key}, {a.n_obs:,} cells, {len(have)} clusters")

    # cells of this lineage removed in earlier L2 rounds
    l1_cells = labels.index[(labels["L1_status"] == "kept") & (labels["L1_lineage"] == lin)]
    for j in range(2, k + 1):
        f = round_dir(lin, j) / f"{lin}_round{j}_removed_cells.csv"
        rm = pd.read_csv(f, dtype=str)
        labels.loc[rm["cell_id"], "final_status"] = "removed_L2"
        labels.loc[rm["cell_id"], "L2_round"] = str(j - 1)
        labels.loc[rm["cell_id"], "L2_removal_reason"] = rm["reason"].values
    # cells dropped by the L2 Scrublet filter (round 1, before clustering; not written to a removed-cells file)
    thr = SCRUBLET_MAX.get(lin, SCRUBLET_MAX_DEFAULT)
    if "scrublet_score" in a.obs and (a.obs["scrublet_score"].astype(float) >= thr).any():
        raise SystemExit(f"{lin}: round {k} still has cells with scrublet_score >= {thr}; set SCRUBLET_MAX to the "
                         "threshold used in l2_scanpy_lineage.py")
    missing = pd.Index(sorted(set(l1_cells) - set(a.obs_names) - set(labels.index[labels["final_status"] == "removed_L2"])))
    if len(missing):
        src = sc.read_h5ad(ROOT / "L1" / "lineages" / f"mouse_2405072_{lin}_L1.h5ad", backed="r")
        score = src.obs["scrublet_score"].astype(float).reindex(missing)
        src.file.close()
        dbl = missing[~(score < thr).values]  # NaN counts as filtered, as in the L2 script
        labels.loc[dbl, "final_status"] = "removed_L2"
        labels.loc[dbl, "L2_round"] = "1"
        labels.loc[dbl, "L2_removal_reason"] = [f"round1 L2 Scrublet filter: scrublet_score {x:.3f} >= {thr}"
                                                for x in score.loc[dbl]]
        say(f"   {len(dbl):,} cells removed by the L2 Scrublet filter (score >= {thr})")
        missing = missing.difference(dbl)
    if len(missing):
        sc_ = score.loc[missing]
        raise SystemExit(f"{lin}: {len(missing)} L1 cells are neither in round {k}, nor in a removed-cells file, nor "
                         f"above the Scrublet threshold {thr} (their scores: min {sc_.min():.3f}, max {sc_.max():.3f})")
    assert set(a.obs_names) <= set(l1_cells), f"{lin}: round {k} has cells outside the L1 lineage"

    m = cl.map(lambda c: cfg["clusters"][c])
    rem = (m == R).values
    labels.loc[a.obs_names, "L2_round"] = str(k)
    labels.loc[a.obs_names, "L2_cluster_key"] = key
    labels.loc[a.obs_names, "L2_cluster"] = cl.values
    labels.loc[a.obs_names, "L2_cluster_uid"] = f"{lin}|round{k}|{key}|" + cl.values
    kept_ids = a.obs_names[~rem]
    rm_ids = a.obs_names[rem]
    labels.loc[rm_ids, "final_status"] = "removed_L2"
    labels.loc[rm_ids, "L2_removal_reason"] = [f"round{k}_{key}_cluster_{c}: removed by user ({today})"
                                               for c in cl.values[rem]]
    vals = np.array([v for v in m.values[~rem]], dtype=object)
    labels.loc[kept_ids, "final_status"] = "kept"
    labels.loc[kept_ids, "user_label_raw"] = [v[0] for v in vals]
    labels.loc[kept_ids, "cell_type"] = [v[1] for v in vals]
    labels.loc[kept_ids, "L2_identity"] = [v[2] for v in vals]
    labels.loc[kept_ids, "L2_state"] = [v[3] for v in vals]
    for c, v in sorted(cfg["clusters"].items(), key=lambda x: int(x[0])):
        decisions.append({"lineage": lin, "round": k, "cluster_key": key, "cluster": c, "n_cells": int((cl == c).sum()),
                          "action": "remove" if v == R else "keep", "user_label_raw": "" if v == R else v[0],
                          "cell_type": "" if v == R else v[1], "L2_identity": "" if v == R else v[2],
                          "L2_state": "" if v == R else v[3], "decided_on": today})

    b = a[~rem].copy()
    counts = b.layers["counts"] if "counts" in b.layers else b.X
    obs = b.obs[[c for c in b.obs.columns if not c.startswith("leiden_L2_") and not c.startswith("round")]].copy()
    obs = obs.rename(columns={c: f"L1_{c}" for c in obs.columns if c.startswith("leiden_")})
    obsm = {"X_umap_lineage": np.asarray(b.obsm["X_umap"], dtype=np.float32)}
    if "X_umap_L1" in b.obsm:
        obsm["X_umap_L1"] = np.asarray(b.obsm["X_umap_L1"], dtype=np.float32)
    parts.append(ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32), obs=obs,
                            var=b.var[[c for c in ["gene_ids", "feature_types", "gene_symbols"] if c in b.var]].copy(),
                            obsm=obsm))
    del a, b

# ------------------------------------------------------------------ assemble
names = [p.var_names for p in parts]
assert all(n.equals(names[0]) for n in names), "gene order differs between lineages"
A = ad.concat(parts, join="inner", merge="first")
A.var = parts[0].var.copy()
assert A.obs_names.is_unique
lab_cols = ["L1_lineage", "L1_call", "L1_state", "L2_round", "L2_cluster_key", "L2_cluster", "L2_cluster_uid",
            "user_label_raw", "cell_type", "L2_identity", "L2_state"]
for c in lab_cols:
    A.obs[c] = labels.loc[A.obs_names, c].values
A.obs["L1_lineage_name"] = A.obs["L1_lineage"].map(LINEAGE_NAME)
assert (labels.loc[A.obs_names, "final_status"] == "kept").all()
assert (A.obs["cell_type"] != "").all(), "kept cell without a cell_type"
A.layers["counts"] = A.X.copy()
L_ = A.copy()
sc.pp.normalize_total(L_)
sc.pp.log1p(L_)
A.layers["lognorm"] = L_.X
del L_
for c in A.obs.columns:
    if A.obs[c].dtype == object or str(A.obs[c].dtype) in ("string", "str"):
        A.obs[c] = A.obs[c].astype(str).astype("category")
A.uns["annotation_info"] = {"built": today, "script": "final_annotation.py",
                            "X": "raw counts; layers counts (same) and lognorm (normalize_total + log1p)",
                            "X_umap_L1": "whole-lung L1 UMAP", "X_umap_lineage": "UMAP of the annotated L2 round, "
                            "computed per lineage (not comparable across lineages)",
                            "rounds": {lin: cfg["round"] for lin, cfg in ANNOTATION.items()}}

# ------------------------------------------------------------------ accounting + outputs
labels["final_status"] = labels["final_status"].replace("", "removed_unaccounted")
assert not (labels["final_status"] == "removed_unaccounted").any(), "cells without a final status"
assert (labels["final_status"] == "kept").sum() == A.n_obs
labels.reset_index().to_parquet(OUT / "mouse_2405072_L2_labels_all_cells.parquet", index=False)
pd.DataFrame(decisions).to_csv(OUT / "L2_decisions_final.csv", index=False)
comp = pd.crosstab([A.obs["L1_lineage_name"], A.obs["cell_type"]], A.obs["sample_id"], margins=True)
comp.to_csv(OUT / "L2_cell_type_counts_by_sample.csv")

# random cells must carry the same counts as the L1 lineage files
rng = np.random.default_rng(0)
for lin in ANNOTATION:
    ids = A.obs_names[A.obs["L1_lineage"] == lin]
    if not len(ids):
        continue
    pick = rng.choice(ids, min(50, len(ids)), replace=False)
    src = sc.read_h5ad(ROOT / "L1" / "lineages" / f"mouse_2405072_{lin}_L1.h5ad", backed="r")
    x_src = src[pick].to_memory().X
    x_out = A[pick].X
    assert abs(sp.csr_matrix(x_src) - sp.csr_matrix(x_out)).sum() == 0, f"{lin}: counts differ from the L1 file"
    src.file.close()

out_h5 = OUT / "mouse_2405072_L2_final.h5ad"
A.write_h5ad(out_h5, compression="gzip")

if "X_umap_L1" in A.obsm:
    A.obsm["X_umap"] = A.obsm["X_umap_L1"]
    sc.pl.umap(A, color=["L1_lineage_name", "sample_id"], wspace=0.35, save="_final_L1_lineage.png", show=False)
    sc.pl.umap(A, color="cell_type", legend_fontsize=6, save="_final_cell_type.png", show=False)
for lin in ANNOTATION:
    s_ = A[A.obs["L1_lineage"] == lin].copy()
    if s_.n_obs:
        s_.obsm["X_umap"] = s_.obsm["X_umap_lineage"]
        sc.pl.umap(s_, color=["cell_type", "L2_state"], wspace=0.5, legend_fontsize=7,
                   save=f"_final_{lin}.png", show=False)

say("\nfinal status of every cell:\n" + labels["final_status"].value_counts().to_string())
say("\ncell types (kept):\n" + comp.to_string())
say(f"\nwrote {out_h5} ({A.n_obs:,} cells x {A.n_vars:,} genes)")
