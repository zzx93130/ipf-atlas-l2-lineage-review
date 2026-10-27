# IPF atlas (GSv2 / scIPF) - project facts for lineage tracks

## Paths
- ROOT = <path to IPF_ATLAS/data/processed on your machine>
- Build spec: ROOT/code/atlas_gsv2_scIPF/spec.md (scIPF) ; ROOT/code/atlas_gsv2/atlas_gsv2_spec.md (GSv2, all modalities)
- Counts: ROOT/atlas_gsv2/core/shards/<DatasetID>.h5ad (18,120 genes = 17,934 coding + 186 IG/TR; raw counts).
  Reader: `sys.path.insert(0, ROOT+"/code/atlas_gsv2/01_integration"); import gsv2_lib as G; G.read_cells(cell_ids=..., genes=...)`
- Measured-aware detection: ROOT/atlas_gsv2/core/gene_measured_by_dataset.csv (102 files hold gene subsets; unmeasured != 0)
- L1 units (scIPF): ROOT/atlas_gsv2_scIPF/L1/units/<unit>_cells.parquet (epithelial, endothelial, stromal, myeloid, tnk, bplasma); labels ROOT/atlas_gsv2_scIPF/L1/annotation/L1_labels.parquet
- Engine labels: ROOT/atlas/engines/celltypist_labels.parquet, hlca_labels.parquet (key cell_id = `<GSM>|<barcode>`). HLCA must not vote on plasma or B/T routing.
- Marker sources: ROOT/atlas/final/marker_evidence/marker_panel_genes_with_sources.csv; ROOT/atlas/L2/<lineage>/marker_dictionary.csv; ROOT/atlas/L1/tables/L1_marker_panel.csv; panel evaluation artifact panel_eval_gene_level.csv (894e37c9).
- Reusable code: ROOT/code/atlas_gsv2_scIPF/03_epithelial/ (eplib.py, 01_preprocess, 02_integrate_cluster, 03_evidence, 04_review, 05_readme, 07_v2_figures, 08_v2_de) is the reference implementation of this skill; GSv2 originals per lineage in ROOT/code/atlas_gsv2/03_<lineage>/.
- Writable: scripts ROOT/code/atlas_gsv2_scIPF/03_<lineage>/ ; outputs ROOT/atlas_gsv2_scIPF/L2/<lineage>/ . Everything else read-only. Never touch the Claude working folder.

## Environment
- env `ipf-atlas` (scanpy 1.11.5, anndata 0.12, scvi-tools 1.4, harmonypy 2.1, torch 2.14). Fresh processes:
  `NUMBA_CACHE_DIR=/tmp/numba_<track> PYTHONPATH=<script dir>:ROOT/code:ROOT/code/atlas_gsv2/01_integration <ipf-atlas env prefix>/bin/python script.py`
- `sys.modules['tensorflow'] = None` before umap/scvi (gsv2_lib.setup_env does it). scanpy wilcoxon resets numba threads: call setup_env(n) in DE scripts.
- harmonypy 2.1: call `harmonypy.run_harmony` directly (sc.external wrapper transposes wrongly); its compute_lisi is O(n^2) - use kNN-based LISI.
- `np.isin` on object strings is very slow: use `pd.Index.isin`. Long jobs: bash background:true with a log file (ps/top blocked).
- Machine 28 cores / 96 GB. Observed peaks: preprocess 429k cells 22 GB; scVI 14 GB; full-data DE 406k cells x 18,120 genes in 4 gene blocks 18-20 GB / ~5 min; UMAP 400k cells ~11 min.

## Annotation principles (user, binding)
1. Elevated stress/IFN/other state genes = a functional state of the compartment's type, never a distinct type. Closed state list: Proliferating, IFN response, Stress response, Metallothionein-high, Mitochondrial-high (two different states).
2. Small clusters are kept only if top-10 (then 20, 30) score-ranked DEGs overlap established markers of a program or rare type.
3. DE = scanpy rank_genes_groups ranked by scores.
4. A cluster whose top DE genes are another lineage's genes is contamination and may be removed.

## User filtering history (what they actually decided)
- L1 (2026-10-26): removed small tails (proliferating AT2, cycling secretory, single-donor fragments, 16-81 cell clusters) and all low-complexity unresolved cells; KEPT mesothelial (946 cells, CALB2/UPK3B/ITLN1/MSLN) and aberrant-basaloid-like (LAMC2/ITGA2/LAMB3/SFN) after the agent flagged them.
- Epithelial L2 (2026-10-26): agreed with all 35 agent proposals at Leiden 1.0 (8 AT2 clusters merged, 2 ciliated merged, 10 removed = contamination, doublets, ribosomal-only, ambient, single-donor, below-drop fragments). Then asked to recluster kept cells and re-run DE to confirm.

## Lineage-specific pitfalls
- Epithelial: Harmony fails G4 (AT1 vs AT2 0.15; AT2 vs pre-TB 0.06) -> scVI. SFTPC/SCGB1A1/SCGB3A1 are ambient (37-96% in every label). Transitional AT2 (CLDN4/SFN), secretory (KRT19/WFDC2) and aberrant basaloid (KRT17/LAMB3/CDH2) lack a +15-point discriminating gene vs suprabasal. GSE157996 (deep FACS-sorted) forms its own LAMC2/ITGA2/PPP1R14B clusters, control and IPF alike. No tuft cells in this atlas.
- Endothelial: APLN marks aCap, not tip cells (ESM1 <= 8%); ACKR1 does not separate pulmonary from systemic venous EC; PDPN is higher in mesothelium than LEC. G4 pair: LEC vs blood EC.
- Stromal: curated other-compartment dictionary lacks DCN/COL1A2/LUM/BGN/C1S/COL6A2/TPM2 -> supplement contamination genes from L1 DE and report separately. Mesothelial is Stromal (KRT8/18/19 high; contamination flag only via HP/other). G4 pairs: pericyte vs SMC, alveolar vs adventitial FB. Marker panel v2 gaps: myofibroblast (ACTA2 -47 vs SMC), lipofibroblast, fibromyocyte lack own genes; SMC MYH11 works.
- Myeloid: pDC = CLEC4C/LILRA4 (GZMB+JCHAIN definition unsafe: GZMB 91% in CD16+ NK, JCHAIN 73% in plasma); mast KIT/HPGDS/TPSAB1; IM SELENOP, classical monocyte VCAN; SPP1 is bimodal (a subtype-like program); IL1B/CXCL8 higher in AM/monocytes than IM. Low-complexity L1 cluster 0.5-7 (18,746, GSE227136-rich) was kept at L1 - check it. G4 pair: mast vs other myeloid.
- T/NK: the L1 Harmony latent mixes NK with CD8 (kNN 0.31 vs scVI 0.07) - check cluster-level NK/CD8 co-clustering before accepting a latent. IL7R is 66-82% in naive/CM/MAIT/TRM, so it cannot mark CM/CD56bright NK/ILC. GSE326573 is CITE-seq sorted T cells (122k) - report its share per cluster. TR_C genes are in the gene space (TRAC, TRBC1/2, TRDC, TRGC1/2); TRAV1-2 is not.
- B/plasma: weakest L1 lineage (CellTypist 0.86). Plasma needs MZB1/JCHAIN/XBP1 (CD79A near-absent in plasma); GC B = RGS13/IRAG2 (LRMP). IG_C genes (IGKC, IGHG1-4, IGHA1/2, IGHM, IGHD) are in the gene space; B->T doublets (CD3E/IL32) and myeloid (LYZ/TYROBP) contamination clusters recur. G4 pair: plasma vs B.
