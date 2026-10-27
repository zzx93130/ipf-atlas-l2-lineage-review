---
name: ipf-atlas-l2-lineage-review
description: "Run the L2 (minor cell type) clustering and annotation review loop for one major lineage of the IPF / GSv2 / scIPF lung scRNA-seq atlas, scanpy style: lineage HVGs, Harmony vs scVI with separation gates, Leiden ladder, full-data Wilcoxon DE, the user's size-drop + DEG filtering rule, preliminary labels with keep/merge/remove recommendations, the fixed review package (cluster and label UMAPs, canonical-marker dot plot by subtype, top-DEG dot plot, review table with blank user columns, README), and after the user agrees: record decisions, per-cell labels, recluster kept cells, re-run DE. Use for any lineage track (epithelial, endothelial, stromal, myeloid, T/NK, B/plasma) of this atlas, and whenever the user asks for L2, minor cell types, subtype clustering, a review package, DEGs and markers to review, reclustering after removal, or agent-proposed labels to compare against. Not for L1 calling, QC or single-dataset annotation."
metadata:
  authored_via: "claude-science"
  published_by: "agent"
  published_at: "2026-10-26T13:28:17.459Z"
  authored_by: "agent"
  authored_at: "2026-10-26T13:21:16.208Z"
  authoring_session: "95d0e9a7-4fa3-4566-a47d-990864a90457"
  last_modified_by: "agent"
  last_modified_at: "2026-10-26T13:28:17.444Z"
---

# L2 lineage review loop (IPF atlas)

One major lineage at a time. The agent clusters, computes evidence, applies the
user's filtering rule, proposes labels, and delivers a fixed package. The user
compares their judgement with the proposals and iterates; nothing becomes a
label until the user agrees. Project paths, env rules and lineage-specific
pitfalls are in `references/ipf_atlas_scipf.md` - read it first.

## Why the loop looks like this

- The user judges clusters from **three things together**: cluster size and
  where the size distribution drops, the score-ranked DEGs, and whether those
  DEGs name a known cell state or rare subtype. Give all three for every
  cluster, in one table and two dot plots, so the comparison is fast.
- Labels are the user's. Agent proposals exist to be compared against, so each
  carries a confidence and a one-line reason naming the deciding genes.
- Established rare types (PNEC, ionocyte, tuft, deuterosomal, aberrant
  basaloid, mesothelial, pDC, mast, GC B, LEC, etc.) are small and the size
  rule would remove them. The user has kept such clusters when the signature
  was shown and removed them when it was not; so **flag the signature with its
  genes explicitly** and let the user decide.
- Dataset-dominated clusters (one dataset >= ~80-95%) are common in a
  multi-study atlas. A distinct program carried by one dataset is not a cell
  type until another dataset shows it: say so in the reason.

## Stage A - integrate and cluster

1. Cells = the lineage's kept L1 unit list. Counts via the project reader
   (raw counts, full gene space; unmeasured genes are NOT zero - use the
   measured-by-dataset table for detection percentages).
2. Lineage HVGs: scanpy `seurat_v3`, `batch_key=SampleID`, ranked
   nbatches-first (seurat_v3_paper ordering), IG/TR V/J clonotype genes
   excluded. normalize_total + log1p, PCA 50.
3. Two latents: PCA+Harmony(SampleID) (call `harmonypy.run_harmony`
   directly) and unit scVI (batch SampleID, ~20 epochs on 400k cells). Gates:
   G1 iLISI(SampleID, DatasetID) above unintegrated PCA; G3 dataset-dominated
   cluster fraction <= 0.10 at Leiden 1.0; G4 known-pair kNN separation
   (lineage-specific pairs, worst pair <= max(0.05, 1.5x PCA)). Prefer Harmony
   when it passes (user prefers scanpy style); otherwise scVI and say why.
   Report both latents' gate values.
4. k=15 neighbours, UMAP, Leiden ladder 0.5/0.8/1.0/1.5/2.0. Primary =
   coarsest resolution that resolves every engine-consensus population with
   stability >= 0.85; name a finer one (stability >= 0.80). Keep the ladder as
   obs columns.

## Stage B - evidence at the primary resolution

- Full-data Wilcoxon cluster-vs-rest (`sc.tl.rank_genes_groups`, ranked by
  `scores`) over the full gene space, in gene blocks normalised with all-gene
  totals, own BH; top-30 per cluster with logFC, padj, measured-aware
  pct_in/pct_out. No subsampling: the user confirms clusters from these lists.
- Per cluster: n_cells, size rank, ratio to the previous (larger) cluster,
  n_donors, n_datasets, top dataset and donor share, frac_IPF (or condition
  share), median genes/counts/mt, doublet scores, ambient score, other-lineage
  genes in top 10/30 (contamination), closed-list state genes in top 10/30,
  CellTypist and HLCA majority + fraction.

## Stage C - the user's filtering rule and the proposals

Sort clusters by size, descending. Find the first position where
`n / n_previous < 0.5` (the user's rule has no number; 0.5 is the working
threshold - state it). Clusters from that drop onward are recommended for
removal **unless** their top DEGs show a distinct state or an established
rare subtype. Above the drop, recommend removal for contamination (another
lineage's genes in the top DE), doublets (doublet score + two programs, no own
program), ambient/low-quality (few genes, high ambient score, other-lineage
genes), ribosomal-only "low information", and near-single-donor clusters.
Recommend `merge_with:<cl>` when a cluster shares the target's identity genes
and has no distinct DEGs (dataset splits, canonical-only duplicates, NF-kB /
immediate-early or lipogenic variants of the same type).

Per cluster give: `proposed_label` (HLCA-style nomenclature for the lineage),
`proposed_state` (closed list only: Proliferating, IFN response, Stress
response, Metallothionein-high, Mitochondrial-high; blank otherwise - states
never enter a label), `recommendation` (keep / remove / merge_with:<cl>),
`confidence` (high/medium/low), `reason` (deciding genes with pct_in/pct_out,
size and size-drop position, dataset/donor dominance). Do not hard-code
markers from memory: take them from the project marker tables, adding
literature markers with a source column only when a proposed subtype has none.

## Stage D - deliverables (figure-style skill; save each as an artifact)

1. `umap_<lin>_clusters.png` - primary Leiden, cluster numbers on the map.
2. `umap_<lin>_proposed_labels.png` - proposed labels; clusters recommended
   for removal in grey, legend with cell counts.
3. `<lin>_marker_dotplot.png` - canonical markers grouped and labelled by the
   subtype they mark; rows = clusters ordered by proposed label, label
   `cluster | n | proposed label`; dot size = % detected over datasets that
   measure the gene; colour = z of the mean across the rows shown.
4. `<lin>_top_DEG_dotplot.png` - top-5 (or top-10) score-ranked DEGs per
   cluster, padj<0.05, logFC>0, own-row genes outlined; same rows.
   Plus `<lin>_DE_top30.csv`.
5. `<lin>_review_table.csv/.xlsx` - one row per cluster with all Stage B
   evidence, Stage C proposals, and blank `user_label`, `user_state`,
   `user_action`, `user_notes`.
6. `README_<lineage>.md` - latent + gates, resolution choice, counts of
   keep/merge/remove, proposed label set, flagged rare types.

In the chat summary: a short table of the clusters that most need the user's
judgement (rare types, single-dataset programs, low-confidence calls), the
counts, and the two things the user must know about the evidence (DE basis,
how the size rule was operationalised).

## Stage E - after the user decides

1. Write `review/user_decisions_<lineage>_<res>.csv`: one row per cluster with
   the proposal columns and the user's `user_label / user_state / user_action
   / user_notes` (when the user agrees wholesale, copy the proposals and note
   the date).
2. Per-cell labels `<lineage>_L2_labels.parquet`: cell_id, the primary Leiden
   column (string), label (merged clusters take the target's label), state,
   status kept/removed, removal reason. Removed cells keep their old label in
   a `_before_removal` column.
3. Recluster the kept cells with the unchanged Stage A pipeline (new HVGs,
   retrained latents, same gates and resolution rule). Move the previous
   outputs to `prev_<date>/` - never delete. Deliver the new UMAP coloured by
   the agreed labels, the new UMAP by new clusters, and a new-cluster x label
   crosstab with majority label and purity; name clusters with purity < 0.8.
4. On request, re-run the full-data DE on the new clusters and deliver a
   top-10 DEG dot plot (rows `cluster | n | majority label (purity)`). Read
   it critically: name clusters whose top 10 lack identity genes, and clusters
   that are ~95% one dataset.

## Hygiene

- Never `rm`, never `delete_host_files`: list stray files in the report.
- Read cluster ids as strings (sub-cluster ids such as `16.10` collapse as
  floats).
- Run scripts as fresh processes with the env recipe in the reference; log
  per stage with peak RSS; keep inside the RAM budget stated in the brief.
- Report deviations from the brief explicitly (subsampling, substituted
  marker sources, gates dropped, thresholds chosen).
