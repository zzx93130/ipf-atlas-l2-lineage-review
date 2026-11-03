"""Standard scanpy analysis of one L1 lineage h5ad (mouse 2405072).

Usage:  python l2_scanpy_lineage.py stromal                       round 1 (from the L1 lineage h5ad)
        python l2_scanpy_lineage.py stromal --remove 3,7         drop clusters 3 and 7 of the latest round
                                                                 (ids at MAIN_RES) and recluster -> round 2, 3, ...
        python l2_scanpy_lineage.py stromal --remove 3,7 --from-round 1 --from-res 0.8 --reason "doublets"
        Lineages: epithelial | endothelial | stromal | myeloid | tnk | bplasma.
        In Jupyter: set LINEAGE / REMOVE_CLUSTERS / FROM_ROUND / REMOVE_FROM_RES / REMOVE_REASON below.
Each later round restarts from raw counts of the kept cells (new HVGs, PCA, Harmony, UMAP, Leiden); outputs go to
L2/<lineage>/round<k>/ (round 1 stays in L2/<lineage>/); an existing round<k>/ is renamed, never deleted.
Steps:  counts -> Scrublet filter -> normalize_total + log1p -> HVGs (seurat_v3_paper, batch = sample_id;
        no mt/ribo/Hb genes; Ig/TCR V/J genes kept)
        -> regress_out(total_counts, pct_counts_mt) -> scale -> PCA 50 -> Harmony(sample_id) -> kNN 15 -> UMAP
        -> Leiden 0.5/0.8/1.0/1.5/2.0 -> Wilcoxon DE (scores, RP/MT removed) -> figures + tables + h5ad.
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_mouse_L2")
sys.modules["tensorflow"] = None  # umap imports TF at load; segfaults on macOS arm64

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ settings
LINEAGE = "stromal"
REMOVE_CLUSTERS = []         # e.g. ["3", "7"]: drop these clusters of the previous round, then recluster
FROM_ROUND = None            # round the ids refer to (None = latest round on disk)
REMOVE_FROM_RES = None       # resolution the ids refer to (None = MAIN_RES)
REMOVE_REASON = "removed by user"
if "ipykernel" not in sys.modules:  # command line overrides the values above
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("lineage")
    ap.add_argument("--remove", default="", help="comma-separated cluster ids to drop before reclustering")
    ap.add_argument("--from-round", type=int, default=None)
    ap.add_argument("--from-res", type=float, default=None)
    ap.add_argument("--reason", default=REMOVE_REASON)
    _a = ap.parse_args()
    LINEAGE, FROM_ROUND, REMOVE_FROM_RES, REMOVE_REASON = _a.lineage, _a.from_round, _a.from_res, _a.reason
    REMOVE_CLUSTERS = [x.strip() for x in _a.remove.split(",") if x.strip()]
ROOT = Path(os.environ.get("MOUSE_ROOT", "/Users/joe/SSD2_1T/IPF_ATLAS/data/processed/mouse_2405072"))
IN_FILE = ROOT / "L1" / "lineages" / f"mouse_2405072_{LINEAGE}_L1.h5ad"
BASE = ROOT / "L2" / LINEAGE
BATCH = "sample_id"          # one sample per condition here: Harmony on it can also remove BLM signal
USE_HARMONY = True           # False -> cluster on the regressed PCA
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
RESOLUTIONS = [0.5, 0.8, 1.0]
MAIN_RES = 0.8 if LINEAGE == "endothelial" else 1.0
SCRUBLET_MAX = 0.07 if LINEAGE == "endothelial" else 0.2   # keep cells with scrublet_score < this (round 1 onwards);
                                                          # final_annotation.py SCRUBLET_MAX must match
SEED = 0

# L2 markers: references/marker_panel_mouse_v2.md (zzx93130/ipf-atlas-l2-lineage-review; gene lists as of ede2050), ## L2 lists
# incl. nested lists (ILC2/3/1, gdT17, MZB, B1) and LEC from ## L1; "?" genes and free text skipped,
# "+" genes kept. STATES / PROGRAMS from ## L3. Genes absent from the data are skipped and printed.
MARKERS = {
    "epithelial": {
        "AT1": ['Rtkn2', 'Hopx', 'Aqp5', 'Pdpn', 'Emp2'],
        "AT2": ['Sftpc', 'Sftpb', 'Sftpa1', 'Sftpd', 'Lamp3', 'Abca3', 'Slc34a2'],
        "Secretory": ['Krt19', 'Cxcl17'],
        "Club": ['Scgb3a2', 'Cyp2f2', 'Gabrp', 'Retnla'],
        "Aberrant basaloid": ['Krt8', 'Cldn4', 'Krt19', 'Sfn', 'Lgals3', 'Cdkn1a', 'Krt18'],
        "Basal": ['Krt5', 'Trp63', 'Krt15', 'Krt17', 'Aqp3'],
        "Goblet": ['Muc5b', 'Ltf', 'Bpifb1', 'Scgb1a1'],
        "Multiciliated": ['Foxj1', 'Rsph1', 'Ccdc153'],
        "Deuterosomal": ['Deup1', 'Foxn4', 'Cdc20b', 'Ccno', 'Mcidas'],
        "PNEC": ['Ascl1', 'Calca', 'Chga', 'Chgb', 'Scg5', 'Syp', 'Resp18'],
        "Ionocyte": ['Foxi1', 'Ascl3', 'Atp6v1b1', 'Atp6v0a4'],
        "Tuft": ['Pou2f3', 'Trpm5', 'Dclk1', 'Avil', 'Gnat3', 'Il25', 'Alox5ap'],
    },
    "endothelial": {
        "aCap": ['Car4', 'Ednrb', 'Apln', 'Tbx2'],
        "gCap": ['Ntrk2', 'Cdk19'],
        "Arterial": ['Gja5', 'Gja4', 'Bmx'],
        "Venous": ['Vcam1', 'Selp', 'Vwf', 'Slc6a2'],
        "Systemic venous": ['Col15a1', 'Vwa1', 'Plvap'],
        "Tip / angiogenic": ['Esm1', 'Kcne3', 'Angpt2', 'Dll4', 'Pgf', 'Cxcr4', 'Mcam'],
        "Lymphatic EC": ['Prox1', 'Mmrn1', 'Ccl21a', 'Pdpn', 'Lyve1'],
    },
    "stromal": {
        "Alveolar fibroblast": ['Npnt'],
        "Adventitial fibroblast": ['Pi16'],
        "Peribronchial fibroblast": ['Hhip', 'Aspn', 'Fgf18', 'Actg2', 'Lgr6'],
        "Myofibroblast": ['Cthrc1', 'Postn', 'Tnc', 'Acta2', 'Col1a1', 'Ltbp2'],
        "Lipofibroblast": ['Plin2', 'Lpl', 'Tcf21', 'Fgf10', 'Ces1d'],
        "SMC": ['Myh11', 'Acta2', 'Tagln', 'Cnn1'],
        "Pericyte": ['Rgs5', 'Higd1b', 'Cox4i2'],
        "Mesothelial": ['Msln', 'Upk3b', 'Wt1'],
    },
    "myeloid": {
        "Alveolar macrophage": ['Siglecf', 'Ear2', 'Chil3', 'Krt79', 'Plet1'],
        "Interstitial macrophage": ['C1qa', 'C1qb', 'C1qc', 'Mrc1', 'Pf4', 'Cx3cr1', 'Selenop', 'Rnase1'],
        "Monocyte-derived macrophage": ['Lyz2', 'Cd68', 'Mertk', 'Fcgr1', 'Ccr2', 'Cx3cr1'],
        "Monocyte": ['Plac8', 'Ifitm3'],
        "Classical monocyte": ['Plac8', 'Vcan', 'Cd14'],
        "Non-classical monocyte": ['Cx3cr1', 'Ace', 'Treml4'],
        "cDC1": ['Clec9a', 'Xcr1', 'Irf8', 'Wdfy4', 'Cadm1', 'Itgae'],
        "cDC2": ['Cd209a', 'Ccl17', 'Mgl2'],
        "Mature / migratory DC": ['Ccr7', 'Cacnb3', 'Fscn1', 'Mreg'],
        "pDC": ['Siglech', 'Bst2', 'Ly6d', 'Tcf4', 'TFF3'],
        "Neutrophil": ['Csf3r', 'Cxcr2', 'S100a8', 'S100a9', 'G0s2'],
        "Mast": ['Tpsab1', 'Tpsb2', 'Cpa3', 'Ms4a2', 'Mcpt4', 'Cma1', 'Kit', 'Gata2'],
        "Basophil": ['Mcpt8', 'Cd200r3', 'Ms4a2', 'Fcer1a'],
        "Megakaryocyte": ['Alox12', 'Tubb1', 'Itga2b', 'Pf4', 'Ppbp', 'Gp1ba', 'Gp1bb', 'Gp9'],
    },
    "tnk": {
        "CD4 T": ['Cd3d', 'Trac', 'Cd4'],
        "CD8 T": ['Cd3d', 'Trac', 'Cd8a', 'Cd8b1'],
        "Treg": ['Foxp3', 'Il2ra', 'Ctla4', 'Ikzf2', 'Tigit'],
        "NK": ['Klrd1', 'Ncr1', 'Nkg7', 'Prf1', 'Ctsw', 'Gzma', 'Klrb1c', 'Eomes', 'Klrk1'],
        "Mature NK": ['Itgam', 'Cx3cr1', 'Klrg1', 'Zeb2', 'S1pr5', 'Fgfbp2', 'Fcgr3'],
        "Immature / tissue NK": ['Cd27', 'Xcl1', 'Sell', 'Itga1', 'Cxcr6'],
        "MAIT": ['Trav1', 'Slc4a10', 'Zbtb16', 'Klrb1c', 'Il18rap', 'Rorc'],
        "Gamma-delta T": ['Trdc', 'Trgc1', 'Trgc2', 'Trgc4'],
        "gdT17": ['Il17a', 'Il17f', 'Rorc', 'Cd163l1', 'Blk', 'Sox13', 'Maf', 'Ccr6'],
        "ILC2": ['Il1rl1', 'Il5', 'Il13', 'Gata3', 'Klrg1', 'Areg', 'Arg1', 'Icos'],
        "ILC3": ['Rorc', 'Il22', 'Il23r', 'Kit', 'Ccr6'],
        "ILC1": ['Tbx21', 'Ncr1', 'Il7r'],
    },
    "bplasma": {
        "Naive B": ['Fcer2a', 'Ighd', 'Ighm', 'Cd19', 'Ms4a1', 'Sell'],
        "Memory B": ['Cd80', 'Pdcd1lg2', 'Nt5e', 'Cd38', 'Ccr6', 'Zbtb32'],
        "MZB": ['Cr2', 'Cd1d1', 'Cd9'],
        "B1": ['Cd5', 'Spn', 'Cd9', 'Itgam'],
        "GC B": ['Rgs13', 'Lmo2', 'Bcl6', 'Aicda', 'Hmces', 'S1pr2', 'Gcsam', 'Fas'],
        "Plasmablast": ['Prdm1', 'Xbp1', 'Mzb1', 'Jchain', 'Tnfrsf17', 'Mki67', 'Slamf7'],
        "Plasma": ['Mzb1', 'Jchain', 'Xbp1', 'Sdc1', 'Tnfrsf17', 'Igkc', 'Igha', 'Ighg1'],
    },
}
STATES = {  # marker_panel_mouse_v2.md ## L3, general states (all lineages)
    "Proliferation": ['Mki67', 'Top2a', 'Ube2c', 'Stmn1'],
    "IFN response": ['Isg15', 'Ifit1', 'Ifit3', 'Mx1', 'Oas1a', 'Rsad2', 'Irf7', 'Usp18', 'Ifi27l2a'],
    "Stress response": ['Fos', 'Jun', 'Junb', 'Egr1', 'Atf3', 'Dusp1'],
    "Heat shock": ['Hspa1a', 'Hspa1b', 'Hspb1', 'Hsp90aa1', 'Dnajb1'],
    "Metallothionein-high": ['Mt1', 'Mt2'],
}
PROGRAMS = {  # marker_panel_mouse_v2.md ## L3, lineage programs
    "myeloid": {"AM IFN response": ['Isg15', 'Rsad2', 'Oasl2'], "Macrophage SPP1-high": ['Spp1', 'Gpnmb', 'Trem2', 'Lgals3', 'Cd9', 'Fabp5', 'Mmp12', 'Chil1'], "Macrophage inflammatory 1": ['Ccl2', 'Il1rn', 'Tnf'], "Macrophage inflammatory 2": ['Il1b', 'Cxcl1', 'Cxcl2', 'Cxcl3']},
    "endothelial": {"EC IFN response": ['Isg20', 'Lap3', 'Cx3cl1', 'Irf1']},
    "tnk": {"T naive": ['Ccr7', 'Sell', 'Lef1', 'Tcf7', 'Foxp1'], "T central memory": ['Ccr7', 'Sell', 'Il7r', 'Tcf7', 'Cd27'], "TRM": ['Cd69', 'Itgae', 'Itga1', 'Cxcr6', 'Zfp683', 'Cd101'], "CD8 Gzmk-high": ['Gzmk', 'Cd27', 'Il7r', 'Ccl5'], "CD8 cytotoxic": ['Gzmb', 'Gzma', 'Prf1', 'Ifng', 'Cx3cr1', 'Tbx21'], "Terminal differentiation": ['Klrg1', 'Klrc1', 'Cx3cr1', 'Zeb2', 'Klrk1']},
    "stromal": {"Fibroblast inflammatory": ['Ccl2', 'Cxcl12', 'Cxcl14', 'Il6', 'Il33'], "Resting": ['Apod', 'Lepr'], "Immune-recruiting": ['Cxcl13', 'Ccl19', 'Ccl21a', 'Il7']},
}
EXCLUDE_HVG = r"^mt-|^(?:Rps|Rpl|Mrps|Mrpl)\d|^Hb[ab]-"   # mito, ribosomal, haemoglobin (Hba-a1, Hbb-bs, ...);
                                                              # Ig/TCR V/J genes are NOT excluded
RP_MT_DE = r"^(?:Rps|Rpl)\d|^Rplp|^Rpsa$|^mt-"

# ------------------------------------------------------------------ run
def round_dir(k):
    return BASE if k == 1 else BASE / f"round{k}"


def round_file(k):
    return round_dir(k) / (f"mouse_2405072_{LINEAGE}_L2_clustered.h5ad" if k == 1
                           else f"mouse_2405072_{LINEAGE}_L2_round{k}_clustered.h5ad")


def latest_round():
    ks = [1] if round_file(1).exists() else []
    ks += [int(p.name[5:]) for p in BASE.glob("round*") if p.name[5:].isdigit() and round_file(int(p.name[5:])).exists()]
    if not ks:
        raise SystemExit(f"no clustered round found under {BASE}: run round 1 first (no --remove)")
    return max(ks)


if REMOVE_CLUSTERS:
    prev_k = FROM_ROUND or latest_round()
    ROUND = prev_k + 1
    prev = sc.read_h5ad(round_file(prev_k))
    ckey = f"leiden_L2_{REMOVE_FROM_RES or MAIN_RES}"
    if ckey not in prev.obs:
        raise SystemExit(f"{ckey} not in round {prev_k}; available: {[c for c in prev.obs if c.startswith('leiden_L2_')]}")
    pc = prev.obs[ckey].astype(str)
    bad = [c for c in REMOVE_CLUSTERS if c not in set(pc)]
    if bad:
        raise SystemExit(f"clusters {bad} not in round {prev_k} {ckey} (has {sorted(set(pc), key=float)})")
    rm = pc.isin(REMOVE_CLUSTERS).values
    OUT = round_dir(ROUND)
    if OUT.exists():  # never delete an earlier attempt
        OUT.rename(OUT.with_name(f"{OUT.name}_prev_{pd.Timestamp.now():%Y%m%d_%H%M%S}"))
    OUT.mkdir(parents=True)
    removed = pd.DataFrame({"cell_id": prev.obs_names[rm], "sample_id": prev.obs["sample_id"].values[rm],
                            "L1_call": prev.obs["L1_call"].astype(str).values[rm], "from_round": prev_k,
                            "cluster_key": ckey, "cluster": pc.values[rm],
                            "reason": f"round{prev_k}_{ckey}_cluster_" + pc.values[rm] + f": {REMOVE_REASON}"})
    removed.to_csv(OUT / f"{LINEAGE}_round{ROUND}_removed_cells.csv", index=False)
    print(f"== {LINEAGE} round {ROUND}: removing clusters {REMOVE_CLUSTERS} of round {prev_k} {ckey} "
          f"({rm.sum():,} cells: {removed['cluster'].value_counts().to_dict()}); {(~rm).sum():,} cells kept")
    obs = prev.obs.loc[~rm, [c for c in prev.obs.columns if not c.startswith("leiden_L2_")]].copy()
    obs[f"round{prev_k}_{ckey}"] = pc.values[~rm]
    adata = sc.AnnData(X=prev.layers["counts"][~rm].copy() if "counts" in prev.layers else prev.X[~rm].copy(),
                       obs=obs, var=prev.var[["gene_ids", "feature_types", "gene_symbols"]].copy()
                       if "gene_ids" in prev.var else None)
    if "X_umap_L1" in prev.obsm:
        adata.obsm["X_umap_L1"] = prev.obsm["X_umap_L1"][~rm]
    adata.layers["counts"] = adata.X.copy()
    PREV_KEY = f"round{prev_k}_{ckey}"
    del prev
else:
    ROUND, OUT, PREV_KEY = 1, BASE, None
    adata = sc.read_h5ad(IN_FILE)
    adata.X = adata.layers["counts"].copy()
TAG = LINEAGE if ROUND == 1 else f"{LINEAGE}_round{ROUND}"
(OUT / "figures").mkdir(parents=True, exist_ok=True)
sc.settings.figdir = OUT / "figures"
sc.settings.set_figure_params(dpi=80, dpi_save=200, frameon=False)
print(f"== {TAG}: {adata.n_obs:,} cells -> {OUT}")
print(adata)
print(adata.obs.groupby(["sample_id"], observed=True).size().to_string())

# 1. normalise (QC metrics come from L1; recomputed here on all genes for the regression covariates)
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
n0 = adata.n_obs
adata = adata[adata.obs["scrublet_score"] < SCRUBLET_MAX].copy()
print(f"Scrublet filter (scrublet_score < {SCRUBLET_MAX}): removed {n0 - adata.n_obs:,} cells, {adata.n_obs:,} kept")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# 2. HVGs within the lineage (excluded gene classes never enter the HVG set)
# seurat_v3's loess fit can fail ("reciprocal condition number") when a small sample has many near-zero genes;
# fall back to a stricter gene filter + wider span, then to the dispersion-based 'seurat' flavor (lognorm).
excl = adata.var_names.str.contains(EXCLUDE_HVG, regex=True)
n_det = np.asarray((adata.layers["counts"] > 0).sum(0)).ravel()
HVG_ATTEMPTS = [("seurat_v3_paper", 3, 0.3), ("seurat_v3_paper", max(10, int(0.005 * adata.n_obs)), 0.5),
                ("seurat", max(10, int(0.005 * adata.n_obs)), None)]
for flavor, min_cells, span in HVG_ATTEMPTS:
    tmp = adata[:, ~excl & (n_det >= min_cells)].copy()
    try:
        if flavor.startswith("seurat_v3"):
            sc.pp.highly_variable_genes(tmp, n_top_genes=N_HVG, flavor=flavor, batch_key=BATCH, layer="counts",
                                        span=span)
        else:
            sc.pp.highly_variable_genes(tmp, n_top_genes=N_HVG, flavor=flavor, batch_key=BATCH)
    except ValueError as e:
        print(f"HVG {flavor} (genes in >={min_cells} cells, span {span}) failed: {str(e).strip()[:80]}")
        continue
    break
else:
    raise SystemExit("all HVG attempts failed")
adata.var["highly_variable"] = adata.var_names.isin(tmp.var_names[tmp.var["highly_variable"]])
adata.uns["hvg_method"] = {"flavor": flavor, "min_cells": int(min_cells), "span": span, "batch_key": BATCH,
                           "n_top_genes": N_HVG}
del tmp
print(f"HVGs: {adata.var['highly_variable'].sum()} ({flavor}, genes in >={min_cells} cells, span {span})")

# 3. regress out depth + mito, scale, PCA (on HVGs only; adata.X keeps the full lognorm matrix)
hv = adata[:, adata.var["highly_variable"]].copy()
sc.pp.regress_out(hv, ["total_counts", "pct_counts_mt"])
sc.pp.scale(hv, max_value=10)
sc.pp.pca(hv, n_comps=50, svd_solver="arpack", random_state=SEED)
adata.obsm["X_pca"] = hv.obsm["X_pca"]
adata.uns["pca"] = hv.uns["pca"]
del hv
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f"_{TAG}.png", show=False)
for k in ["total_counts", "pct_counts_mt"]:
    v = np.log1p(adata.obs[k]) if k == "total_counts" else adata.obs[k]
    r = max(abs(np.corrcoef(adata.obsm["X_pca"][:, i], v)[0, 1]) for i in range(10))
    print(f"max |r| top-10 PCs vs {k}: {r:.2f}")

# 4. batch correction
rep = "X_pca"
if USE_HARMONY and adata.obs[BATCH].nunique() > 1:
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"][:, :N_PCS], adata.obs[[BATCH]], [BATCH], random_state=SEED)
    Z = np.asarray(ho.Z_corr)
    adata.obsm["X_pca_harmony"] = (Z.T if Z.shape[0] != adata.n_obs else Z).astype(np.float32)
    rep = "X_pca_harmony"

# 5. graph, UMAP, Leiden ladder
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep=rep, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
for r in RESOLUTIONS:
    sc.tl.leiden(adata, resolution=r, key_added=f"leiden_L2_{r}", flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED)
    print(f"Leiden {r}: {adata.obs[f'leiden_L2_{r}'].nunique()} clusters")
key = f"leiden_L2_{MAIN_RES}"

# 6. UMAPs
sc.pl.umap(adata, color=[key], legend_loc="on data", legend_fontsize=8, save=f"_{TAG}_clusters.png", show=False)
sc.pl.umap(adata, color=[f"leiden_L2_{r}" for r in RESOLUTIONS], legend_loc="on data", legend_fontsize=6, ncols=3,
           save=f"_{TAG}_leiden_ladder.png", show=False)
meta = [c for c in ["sample_id", "L1_call", "L1_state", "L1_cluster", PREV_KEY] if c and c in adata.obs]
adata.obs["L1_state"] = adata.obs["L1_state"].astype(str).replace({"": "none", "nan": "none"})
sc.pl.umap(adata, color=meta, ncols=2, wspace=0.4, save=f"_{TAG}_metadata.png", show=False)
sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt", "scrublet_score"], ncols=2,
           cmap="viridis", save=f"_{TAG}_qc.png", show=False)

# 7. cluster table: size, sample composition, QC
g = adata.obs.groupby(key, observed=True)
tab = pd.DataFrame({"n_cells": g.size(), "median_genes": g["n_genes_by_counts"].median(),
                    "median_counts": g["total_counts"].median(), "median_pct_mt": g["pct_counts_mt"].median().round(2),
                    "scrublet_median": g["scrublet_score"].median().round(3)})
tab = tab.join(pd.crosstab(adata.obs[key], adata.obs["sample_id"], normalize="index").round(3).add_prefix("frac_"))
tab = tab.join(pd.crosstab(adata.obs[key], adata.obs["L1_call"]).add_prefix("L1_"))
tab = tab.sort_values("n_cells", ascending=False)
tab["ratio_to_prev"] = (tab["n_cells"] / tab["n_cells"].shift(1)).round(2)

# 8. DE: Wilcoxon cluster-vs-rest on lognorm, full gene space, ranked by score; RP/MT genes removed
sc.tl.rank_genes_groups(adata, key, method="wilcoxon", layer="lognorm", use_raw=False, pts=True, key_added="de")
de = sc.get.rank_genes_groups_df(adata, group=None, key="de")
de = de.rename(columns={"names": "gene", "logfoldchanges": "logFC", "pvals_adj": "padj",
                        "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
de = de[(de["padj"] < 0.05) & (de["logFC"] > 0) & ~de["gene"].str.contains(RP_MT_DE, regex=True)]
de = de.sort_values(["group", "scores"], ascending=[True, False])
top30 = de.groupby("group", observed=True).head(30)
top30.to_csv(OUT / f"{TAG}_DE_top30.csv", index=False)
tab["top10_DE"] = [", ".join(f"{r.gene} {r.pct_in:.0%}/{r.pct_out:.0%}" for r in
                             top30[top30["group"] == c].head(10).itertuples()) for c in tab.index]
tab.index.name = key
tab.to_csv(OUT / f"{TAG}_cluster_table.csv")
print(tab[["n_cells", "ratio_to_prev", "median_genes", "top10_DE"]].to_string(max_colwidth=90))
if PREV_KEY:  # where did the previous round's clusters go?
    ct = pd.crosstab(adata.obs[PREV_KEY], adata.obs[key])
    ct.to_csv(OUT / f"{TAG}_prev_vs_new_clusters.csv")
    print(f"previous round clusters ({PREV_KEY}, rows) x new clusters ({key}, columns):\n{ct.to_string()}")

# 9. dot plots: panel markers + state genes; top-5 DEGs per cluster
def present(groups):
    seen, out = set(), {}
    for k, gs in groups.items():
        keep = [x for x in gs if x in adata.var_names and x not in seen]
        seen.update(keep)
        if keep:
            out[k] = keep
    return out

groups = {**MARKERS[LINEAGE], **STATES, **PROGRAMS.get(LINEAGE, {})}
panel = present(groups)
missing = sorted({x for gs in groups.values() for x in gs} - set(adata.var_names))
print(f"panel genes not in the data: {missing}")
order = [str(c) for c in tab.index]
adata.obs[key] = adata.obs[key].cat.reorder_categories(order)
sc.pl.dotplot(adata, panel, groupby=key, layer="lognorm", standard_scale="var", dendrogram=False,
              save=f"{TAG}_marker_dotplot.png", show=False)
top5 = present({c: list(top30.loc[top30["group"] == c, "gene"].head(5)) for c in order})
thin = [c for c in order if c not in top5]
if thin:
    print(f"clusters with no passing DE gene (padj<0.05, logFC>0, RP/MT removed): {thin}")
if sum(len(v) for v in top5.values()) >= 2:
    sc.pl.dotplot(adata, top5, groupby=key, layer="lognorm", standard_scale="var", dendrogram=False,
                  save=f"{TAG}_top5_DEG_dotplot.png", show=False)
else:
    print("top-5 DEG dot plot skipped: fewer than 2 DE genes in total")

# 10. save
adata.X = adata.layers["counts"]  # X back to raw counts; lognorm kept as a layer
del adata.uns["de"]
for c in adata.obs.columns:
    if adata.obs[c].dtype == object:
        adata.obs[c] = adata.obs[c].astype(str)
adata.write_h5ad(round_file(ROUND), compression="gzip")
print(f"== done: {OUT}")
