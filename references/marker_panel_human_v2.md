# Human Lung Single-Cell Annotation Marker Panel (v2)

Recorded verbatim from the user (2026-11-02). Typos in the source are kept and listed under "Source issues" at the end; the mouse panel (`marker_panel_mouse_v2.md`) corrects them.

**Scope:** Human lung scRNA-seq/snRNA-seq, including healthy lung and disease contexts such as interstitial lung disease (ILD) and idiopathic pulmonary fibrosis (IPF).  
**Status:** Curated candidate panel for annotation and review; not a dataset-validated classifier.  
**Principle:** Distinguish lineage (L1), cell identity/subtype (L2), and potentially overlapping transcriptional states (L3).

The markers below are **sets of supporting evidence, not mandatory co-expression lists**. Single-cell dropout, differences between scRNA-seq and snRNA-seq, disease-associated remodeling, and anatomical sampling must be considered. The `Notes` beneath a class are part of its interpretation.

## L1 — Compartment and lineage markers

```yaml
Epithelial:
  markers: [EPCAM]

Blood endothelial:
  markers: [PECAM1, CDH5, CLDN5, VWF, RAMP2]

Lymphatic endothelial:
  markers: [PROX1, FLT4, MMRN1, CCL21, PDPN]

Fibroblast:
  markers: [PDGFRA, DCN, LUM, COL1A2, COL3A1]

Mural:
  markers: [PDGFRB, RGS5, ACTA2, MYH11, TAGLN]

Myeloid:
  markers: [LYZ, TYROBP, FCER1G, AIF1, SPI1]
  neutrophil_support: [S100A8, S100A9, FCGR3B]

Lymphoid:
  T_lineage: [CD3D, TRAC, CD2]
  B_lineage: [CD79A, MS4A1]
  NK_lineage: [KLRD1, NKG7, GNLY, NCR1]

Mast:
  markers: [TPSAB1, TPSB2, CPA3, MS4A2]

Mesothelial:
  markers: [MSLN, UPK3B, CALB2, WT1]
```

**L1 interpretation:** Lymphoid is the union of T-, B-, and NK-lineage evidence, not a single co-expression signature. The mural lineage likewise combines smooth muscle and pericyte markers and must be resolved at L2. Blood and lymphatic endothelial cells may be handled as sibling L1 categories or nested under an overarching endothelial parent, provided the hierarchy remains consistent.

## L2 — Cell identity and subtype markers

### Epithelial

```yaml
AT1:
  markers: [AGER, RTKN2, HOPX]

AT2:
  markers: [SFTPA1, SFTPA2]

AT0:
  markers: [SCGB3A2, SFTPB, SFTPC]
  expected_low: [SCGB1A1]

pre-TB secretory:
  markers: [SCGB3A2, SFTPB, SCGB1A1]
  expected_low: [SFTPC]

Terminal/respiratory bronchiole (TRB) secretory:
  markers: [SCGB3A2, SFTPB]
  expected_low: [SFTPC, SCGB1A1]
  note: Conditional call; establish separation from AT0 and pre-TB cells.

Secretory:
  markers: [KRT19, WFDC2, CXCL17]

Club:
  markers: [SCGB1A1, SCGB3A1, CYP2F1, BPIFB1]

Aberrant basaloid:
  markers: [KRT17, LAMC2]
  required_context: Epithelial identity (EPCAM/keratins); KRT5 low or absent.
  note: Disease-associated; exclude epithelial-stromal doublets or contamination.

Basal:
  markers: [KRT5, TP63, KRT15, KRT17]

Suprabasal:
  markers: [KRT4, KRT13, KRT16, SCEL]

Goblet:
  markers: [MUC5AC, MUC5B, SPDEF, TFF3, AGR2]

Multiciliated:
  markers: [TPPP3, RSPH1, CAP5]

Deuterosomal:
  markers: [DEUP1, FOXN4, CDC20B, CCNO]

PNEC:
  markers: [ASCL1, CALCA, CHGA, CHGB, SCGN]

Ionocyte:
  markers: [FOXI1, CFTR, ATP6V1A, ASCL3]

Tuft:
  markers: [POU2F3, TRPM5, LRMP, HES6]

Submucosal gland (SMG) serous:
  markers: [LTF, LPO, AZGP1, PIP, ZG16B]

Submucosal gland (SMG) mucous:
  markers: [MUC5B, TFF3, PIGR, BPIFB1]
  note: Conditional call; confirm gland-associated identity.
```

### Myeloid and mast

```yaml
Alveolar macrophage (AM):
  markers: [FABP4, MARCO, MRC1, C1QA, C1QB]

Interstitial macrophage (IM):
  markers: [FOLR2, LYVE1, CD163, SELENOP, RNASE1]
  note: Confirm macrophage lineage; exclude LYVE1-positive lymphatic ECs.

Monocyte-derived macrophage (MoMac):
  markers: [LYZ, CD68, CSF1R, MERTK]
  note: Requires evidence of a monocyte-to-macrophage transition. These genes alone do not prove developmental origin.

Classical monocyte:
  markers: [FCN1, VCAN, CD14]

Non-classical monocyte:
  markers: [FCN1, FCGR3A, CX3CR1]
  note: Exclude FCGR3A-positive NK cells.

cDC1:
  markers: [CLEC9A, XCR1, IRF8, WDFY4, CADM1]

cDC2:
  markers: [CD1C, FCER1A, CLEC10A, CD1E, CLEC4A]

Mature / migratory DC:
  markers: [LAMP3, CCR7, FSCN1, LAD1, CD83, CCL19]

Plasmacytoid DC (pDC):
  markers: [CLEC4C, IL3RA, TCF4, PLD4, LILRA4]

Neutrophil:
  markers: [FCGR3B, CSF3R, CXCR2, S100A8, S100A9, G0S2]

Mast:
  markers: [TPSAB1, TPSB2, CPA3, MS4A2]

Eosinophil:
  markers: [CLC, PRG2, RNASE2, SIGLEC8]
  note: Conditional call; require coherent, cluster-level evidence.
```

### T, NK, and ILC

```yaml
CD4 T:
  markers: [CD3D, TRAC, CD4]

CD8 T:
  markers: [CD3D, TRAC, CD8A, CD8B]

Regulatory T (Treg):
  markers: [FOXP3, IL2RA, CTLA4, IKZF2, TIGIT]

NK:
  markers: [KLRD1, NCR1, NKG7, GNLY, PRF1, CTSW]
  note: Confirm that a convincing T-lineage program is absent.

CD16-positive NK:
  markers: [FCGR3A, FGFBP2, PRF1, KLRD1, CX3CR1]

CD56bright NK:
  markers: [NCAM1, XCL1, XCL2, KLRD1]
  note: NCAM1 transcript detection can be sparse.

MAIT:
  markers: [TRAV1-2, KLRB1, SLC4A10, IL18RAP]

Gamma-delta T:
  markers: [TRDC, TRGC1, TRGC2]
  note: Require T-lineage and preferably V(D)J evidence; TCR constant transcripts alone are insufficient.

ILC:
  markers: [IL7R, KIT, RORA, GATA3, KLRB1]
  note: Conditional; establish separation from T and NK lineages.
```

### B and plasma

```yaml
Naive B:
  markers: [TCL1A, FCER2, IGHD, IGHM]

Memory B:
  markers: [CD27, TNFRSF13B, CD44]
  note: Review IGHD/IGHM expression and class-switched immunoglobulin genes.

Germinal center (GC) B:
  markers: [RGS13, LMO2, BCL6, AICDA, HMCES]
  note: AICDA can be sparsely detected.

Plasmablast:
  markers: [PRDM1, XBP1, MZB1, JCHAIN, TNFRSF17]
  note: Often proliferative; may retain some B-lineage markers.

Plasma:
  markers: [MZB1, JCHAIN, XBP1, SDC1, TNFRSF17]

pDC:
  makers: [GZMB, JCHAIN]
  note: From myeloid lineage, but usually cluster with B lineage and emerge as a split cluster in B lineage subclustering.
```

### Endothelial

```yaml
Aerocyte capillary (aCap):
  markers: [CA4, HPGD, EDNRB, EMCN]

General capillary (gCap):
  markers: [RGCC, FCN3, EPAS1]

Arterial EC:
  markers: [GJA5, IGFBP3, CXCL12, DKK2]

venous EC:
  markers: [ACKR1]

Pulmonary venous EC:
  markers: [HDAC9, CPE, CLU1]
  expected_low: [COL15A1]

Systemic venous EC:
  markers: [COL15A1, VWA1, PLVAP, ADAMTS9]

Tip / angiogenic EC:
  markers: [ESM1, APLN, ANGPT2, DLL4]
  note: Conditional angiogenic state; do not call by default.
```

**Venous EC caveat:** Pulmonary versus systemic venous identity can be difficult to separate with a small expression panel. Validate using relative subtype scores, anatomical context when available, and reference mapping. Fibrotic remodeling can alter their spatial distribution.

### Stromal

```yaml
Alveolar fibroblast:
  markers: [NPNT, GPC3, INMT, COL13A1, FMO2]

Adventitial fibroblast:
  markers: [PI16, MFAP5, CD34, IGFBP6]

Peribronchial fibroblast:
  markers: [FGF18, LGR6, WIF1]

Myofibroblast:
  markers: [CTHRC1, POSTN, TNC]

Smooth muscle cell (SMC):
  markers: [MYH11, ACTA2, TAGLN, CNN1, LMOD1, MYL9]

Pericyte:
  markers: [RGS5, CSPG4, PDGFRB, HIGD1B, COX4I2]

Mesothelial:
  markers: [MSLN, UPK3B, CALB2, WT1, KRT8, KRT18]
```

## L3 — cell-state programs

L3 labels are parallel assigned **with L2 identity**. A cluster could be assigned to a cell-state program and keep its lineage identity.

### General programs

```yaml
Stress response:
  markers: [FOS, JUN, JUNB, EGR1, ATF3, DUSP1]

Interferon (IFN) response:
  markers: [ISG15, IFIT1, IFIT3, MX1, OAS1, IFI6]

Proliferation:
  markers: [MKI67, TOP2A, UBE2C, STMN1]

Hypoxia response:
  markers: [HSPA1A, HSPB1, HSP90AA1]
  
Metallothionein-high:
  markers: [MT1E, MT1X, MT2A, MT1G]

```

### Macrophage programs

```yaml
Macrophage SPP1-high:
  markers: [SPP1, MMP9, CHI3L1, CHIT1, LILRB4]

Macrophage inflammatory type1:
  markers: [CCL2, IL1RN, TNF]

Macrophage inflammatory type2:
  markers: [IL1B, CXCL8, CXCL2]
```

### Endothelial programs

```yaml
IFN-response:
  makers: [ISG20, LAP3, CX3CL1, IRF1]
```

### T-cell programs

```yaml
T naive:
  markers: [CCR7, SELL, LEF1, TCF7, FOXP1]

T central memory:
  markers: [CCR7, SELL, IL7R, TCF7, CD27]

T tissue-resident memory (TRM):
  markers: [CD69, ITGAE, ITGA1, CXCR6, ZNF683]

CD8 GZMK-high:
  markers: [GZMK, CD27, IL7R, CCL5]
  note: Require CD8 T identity.

CD8 cytotoxic:
  markers: [GZMB, PRF1, GNLY, CX3CR1, TBX21]
  note: Require CD8 T identity; exclude NK assignment.

Terminally differentiation:
  markers: [KLRG1, KLRC2, KLRC3]
  note: Require CD8 T identity; exclude NK assignment.
```

### Fibroblast programs

```yaml
Fibroblast inflammatory:
  markers: [CCL2, CXCL12, CXCL14, IL6, IL33]

Resting:
  markers: [APOD, LEPR]

Fibroblast immune-recruiting:
  markers: [CXCL13, CCL19, CCL21]
  note: Exclude lymphatic ECs; spatial validation is especially useful for TLS-related interpretation.
```

## Annotation decision rules

| Rule | Required behavior |
| --- | --- |
| **1. Confirm lineage first** | Assign an L2 identity and L3 program only after sufficient L1 evidence supports its parent lineage. |
| **2. Score sibling specificity** | Require markers that discriminate the proposed subtype from competing subtypes within the same compartment, not merely a high absolute signature score. |
| **3. Check differential expression** | Compare the proposed signature with cluster-versus-rest DE genes, detection rates, and effect sizes. |
| **4. Evaluate negative evidence** | Reject or downgrade labels when a competing lineage/subtype explains the transcriptome better. |
| **5. Keep states parallel** | Take L3 programs as parallel to identity. |
| **6. Fail upward** | Return to a validated parent type if subtype evidence is inadequate. |
| **7. Audit small clusters** | Inspect the top 10–30 DE genes, potential doublets, ambient RNA, sample distribution, and biological coherence. |
| **8. Validate externally** | Compare final labels with reference mapping (e.g., HLCA) and manually inspect DotPlots/feature plots for discrepant clusters. |

## Selected primary references

1. Sikkema L, et al. **An integrated cell atlas of the lung in health and disease.** *Nature Medicine* (2023). https://doi.org/10.1038/s41591-023-02327-2
2. Travaglini KJ, et al. **A molecular cell atlas of the human lung from single-cell RNA sequencing.** *Nature* (2020). https://doi.org/10.1038/s41586-020-2922-4
3. Kadur Lakshminarasimha Murthy P, et al. **Human distal lung maps and lineage hierarchies reveal a bipotent progenitor.** *Nature* (2022). https://pmc.ncbi.nlm.nih.gov/articles/PMC9169066/
4. Habermann AC, et al. **Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis.** *Science Advances* (2020). https://pmc.ncbi.nlm.nih.gov/articles/PMC7439444/
5. Schupp JC, et al. **Integrated Single-Cell Atlas of Endothelial Cells of the Human Lung.** (2021). https://pmc.ncbi.nlm.nih.gov/articles/PMC8300155/
6. **Human Lung v2 (HLCA) Reference Explorer.** Azimuth/HuBMAP. https://azimuth.hubmapconsortium.org/references/human_lung_v2/

## Source issues (not changed above)

- `makers` key in pDC (B and plasma) and IFN-response (endothelial programs) should be `markers`.
- `CAP5` (multiciliated) is likely `CAPS`; `CLU1` (pulmonary venous EC) is likely `CLU`.
- "Hypoxia response" lists heat-shock genes (HSPA1A, HSPB1, HSP90AA1), which mainly report dissociation/heat-shock stress, not hypoxia.
- `venous EC` heading is lower-case; "Terminally differentiation" should be "Terminal differentiation"; pDC appears under both myeloid and B/plasma.
