# Mouse 2405072 - L1 lineage analysis

Three samples: `2405072_BLM`, `2405072_BLM_870`, `2405072_Saline` (10x filtered matrices, zipped).
Output units match the human atlas: epithelial, endothelial, stromal, myeloid, tnk, bplasma.

## Files
- `l1_mouse.py` - all stages (`--stage ...`)
- `doublets_r.R` - scDblFinder + DoubletFinder, called by the doublets stage if `Rscript` is on PATH
- marker panel: `../../references/marker_panel_mouse_v2.md` (read by `--stage markers`; a copy next to the script takes precedence)
- `user_decisions_L1.csv` - your L1 decisions for Leiden 0.5 clusters 0-27 (read by `--stage apply`)
- `l2_scanpy_lineage.py` - standard scanpy analysis of one lineage h5ad (`python l2_scanpy_lineage.py <lineage>`), outputs in `L2/<lineage>/`; `--remove 3,7` drops clusters of the latest round and reclusters the rest from raw counts into `L2/<lineage>/round<k>/` (`--from-round`, `--from-res`, `--reason`). Scrublet filter `SCRUBLET_MAX` (0.07 endothelial, 0.2 others); HVGs exclude mt/ribo/Hb genes, Ig/TCR V/J kept; dot plots use the L2 + L3 lists of `marker_panel_mouse_v2.md` (commit 973545f).
- `final_annotation.py` - user L2 annotation -> `final/mouse_2405072_L2_final.h5ad` + per-cell label table for every cell (QC / doublet / L1 / L2 / Scrublet removals with reasons). `SCRUBLET_MAX` must match the L2 script.
- `sample_sheet.csv` - sample_id, zip path (relative to the data dir), condition, batch. **Check `condition`** (what is 870?).

## Run order (fresh process per stage)
```bash
cd <repo>/scripts/mouse_2405072
PY=~/.claude-science/conda/envs/ipf-atlas/bin/python
export NUMBA_CACHE_DIR=/tmp/numba_mouse_L1
$PY l1_mouse.py --stage markers          # reads ../../references/marker_panel_mouse_v2.md
for s in qc doublets embed evidence figures; do $PY l1_mouse.py --stage $s || break; done
# optional L1 round 2, 3, ...: drop L1 clusters and re-embed the rest, then review the new round
$PY l1_mouse.py --stage recluster --remove 3,7 --reason "low quality"
$PY l1_mouse.py --stage evidence && $PY l1_mouse.py --stage figures   # latest round -> L1/review/round<k>/
# after review: decisions in user_decisions_L1.csv (round 1) or user_decisions_L1_round<k>.csv (round k)
$PY l1_mouse.py --stage apply   # L1 labels for every cell (latest round; --round k to choose)
$PY l1_mouse.py --stage split   # one h5ad per lineage in L1/lineages/
for L in epithelial endothelial stromal myeloid tnk bplasma; do $PY l2_scanpy_lineage.py $L; done
$PY l2_scanpy_lineage.py stromal --remove 3,7 --reason "..."   # optional further rounds
$PY final_annotation.py
```
Do not re-run the L1 stages after `apply`: re-clustering renumbers the L1 clusters that `user_decisions_L1.csv` refers to.
Defaults: data `/Users/joe/SSD2_1T/IPF_ATLAS/data/mouse_data`, output
`/Users/joe/SSD2_1T/IPF_ATLAS/data/processed/mouse_2405072` (override with `--data-dir`, `--out-dir`).
Needs `scikit-misc` (seurat_v3 HVGs), `harmonypy`, `leidenalg`/`igraph`, `openpyxl`, `scikit-image` (Scrublet).

## What each stage does
| stage | rule |
|---|---|
| markers | Parses `marker_panel_mouse_v2.md`: the `## L1` block routed to units (Blood + Lymphatic endothelial -> endothelial; Fibroblast + Mural + Mesothelial -> stromal; Myeloid + Mast -> myeloid; Lymphoid T_lineage + NK_lineage -> tnk, B_lineage -> bplasma; pan-immune Ptprc not scored). `+` markers kept, `?` markers not scored. Genes claimed by >1 unit anywhere in L1+L2 are dropped (`multi_unit_claims.csv`). A unit with <3 L1 genes (epithelial: only Epcam) gets the first unique anchor of each of its L2 types. States from `## L3`: Proliferation, IFN response, Stress response (Fos/Jun/Egr1...), Metallothionein-high. `--human-panel` / `--mouse-panel` CSV still work. |
| qc | Unzip; integer counts, `mt-` genes, unique barcodes, identical features across samples; per sample 5 MAD counts/genes, mito 3 MAD capped at 20%, min 200 genes; stops if a sample loses >40%. |
| doublets | Scrublet per sample (rate 0.8%/1k cells) + scDblFinder + DoubletFinder; remove when >=2 of 3 agree (2 of 2 if one R method fails; flag only if Scrublet alone). |
| embed | Median-normalise + log1p; 3,000 seurat_v3_paper HVGs by sample (no mt/ribo/Hb/IG-TR V/J); regress total_counts + pct_mt, scale, PCA 50; Harmony on sample_id; gates G1 iLISI and G4 kNN separation (B-T, myeloid-T, endo-stromal, epi-stromal) for PCA vs Harmony; k=15, UMAP, Leiden 0.2/0.5/1.0. |
| evidence | At Leiden 0.5: AUCell per unit, cluster z + margin, support floor 20%, detection >=40% or margin >=0.75; Wilcoxon cluster-vs-rest (scores, padj<0.05, RP/MT removed); size drop (n/n_prev < 0.5); per-cell dual-lineage test; ambient vs doublet diagnosis; states relative to lineage median; keep/remove/review proposals with reasons. |
| figures | UMAP by cluster, by proposed L1 (removals grey), by sample; panel marker dot plot; top-10 DEG dot plot (own-row outlined, sharing counts); state-gene dot plot; QC violins. |
| recluster | `--remove 3,7` drops Leiden 0.5 clusters of the latest L1 round (`--from-round`, `--from-res`, `--reason`); removed cells go to `work/L1_round<k>_removed_cells.csv`; the rest are re-embedded from raw counts (HVGs, regress/scale/PCA, Harmony, gates, UMAP, Leiden ladder) into `work/03_embedded_round<k>.h5ad`; `L1/review/round<k>/L1_prev_vs_new_clusters.csv`. An existing round<k> is renamed `*_prev_<timestamp>`, never deleted. evidence / figures / apply / split then use the latest round (`--round k` to choose). |
| apply | Reads `user_decisions_L1.csv` (round 1) or `user_decisions_L1_round<k>.csv` (round k; or `--decisions`) (columns cluster, user_label = unit, user_call = finer L1 name, user_state, user_action keep/remove, route, user_notes; every cluster must be listed). `route=per_cell` (clusters 17 Proliferating, 10 IFN response) assigns each cell the unit with the highest L1-panel AUCell z if its margin over the runner-up is >= 1.0 (`ROUTE_MARGIN`); cells below are removed with the reason. Writes `L1/annotation/L1_labels.parquet` (every cell: QC/doublet/user removals with reasons), `L1_state_cluster_routing.csv`, `L1_counts_by_sample.csv`, `L1/units/<unit>_cells.parquet`, `user_decisions_L1_0.5.csv`, UMAPs `umap_L1_final*.png` (lineage, call, state). |
| split | `L1/lineages/mouse_2405072_<unit>_L1.h5ad` for epithelial, endothelial, stromal, myeloid, tnk, bplasma: kept cells, X = raw counts (float32) = layers['counts'], all measured genes; obs = sample/condition, QC, doublet scores, Leiden ladder, AUCell, L1_cluster/lineage/call/state/route; obsm X_pca, X_pca_harmony, X_umap_L1. L2 re-derives HVGs and latents from these counts. |

## Caveats
- One sample per condition: `sample_id` (Harmony key) == condition, and there are no replicates. Harmony can remove BLM-specific signal; the gates are logged for both latents. Set `LATENT = "pca"` to skip Harmony.
- No ambient correction (filtered matrices only). If `raw_feature_bc_matrix` exists, add SoupX/CellBender before QC.
- The embed stage prints panel genes missing from the matrix (the panel asks for symbols to be checked against var_names).
- The panel notes its study attributions were written from memory (PMIDs checked); spot-check key genes.
- No reference engine (CellTypist / scANVI) is run; L1 evidence is the marker panel plus DE. Add a mouse lung reference if you want a third vote.
