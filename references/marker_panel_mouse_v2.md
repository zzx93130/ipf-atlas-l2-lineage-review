# Mouse Lung Single-Cell Annotation Marker Panel (v2-mouse, draft)

Derived from `marker_panel_human_v2.md` by (1) ortholog mapping, (2) review against public mouse lung studies.
**Status:** Curated candidate panel. Before use, confirm every symbol against the dataset's var_names (and MGI `HOM_MouseHumanSequence.rpt`). 
**Scope:** Mouse lung scRNA/snRNA-seq (homeostasis, bleomycin/other fibrosis, injury). 

## Species-level cautions

- Mice have no respiratory bronchioles and submucosal glands exist only in the upper trachea, so AT0, pre-TB, TRB and SMG calls usually do not apply to distal-lung datasets.
- Basal cells are restricted to trachea/main bronchi, except after injury (ectopic Krt5+ cells).
- Many human immune markers are absent in mouse: GNLY, FCGR3B, CD1C, CD1E, LILRA4, XCL2, FCN3, IFI6, CXCL8, KLRC3, CD27-on-B-memory logic.

## L1

```yaml
Epithelial:         {markers: [Epcam]}
Blood endothelial:  {markers: [Pecam1, Cdh5, Cldn5, Vwf, Ramp2]}
Lymphatic endothelial:
  markers: [Prox1, Mmrn1, Ccl21a, Pdpn, Lyve1]
Fibroblast:         {markers: [Pdgfra, Dcn, Lum, Col1a2, Col3a1]}
Mural:              {markers: [Pdgfrb, Rgs5, Acta2, Myh11, Tagln]}
Myeloid:
  markers: [Lyz2, Tyrobp, Fcer1g, Aif1, Spi1]
  neutrophil_support: [S100a8, S100a9, +Ly6g, +Retnlg]
Lymphoid:
  markers: [Ptprc]
  T_lineage: [Cd3d, Trac, Cd2]
  B_lineage: [Cd79a, Ms4a1, +Cd19]
  NK_lineage: [Klrd1, Nkg7, Ncr1, +Klrb1c, +Gzma] 
Mast:               {markers: [Tpsab1, Tpsb2, Cpa3, Ms4a2, +Mcpt4, +Kit]}
Mesothelial:        {markers: [Msln, Upk3b, Calb2, Wt1, +Lrrn4]}
```

## L2

### Epithelial

```yaml
AT1:
  markers: [Ager, Rtkn2, Hopx, +Aqp5, +Pdpn, +Emp2]
AT2:
  markers: [Sftpc, Sftpb, Sftpa1, Sftpd, +Lamp3, +Abca3, +Slc34a2]
  note: SFTPA1/SFTPA2 → single Sftpa1 (no Sftpa2). Sftpc is the standard mouse AT2 anchor; Lyz1/Lyz2 are also high in mouse AT2.
Secretory:
  markers: [Krt19, Wfdc2, Cxcl17, +Bpifa1]
Club:
  markers: [Scgb1a1, Scgb3a2, Cyp2f2, Bpifb1, +Cbr2]
Aberrant basaloid → "Krt8-high transitional AT2→AT1 state" (PATS / DATP / ADI):
  markers: [Krt8, Cldn4, Krt19, Sfn, Lgals3, Cdkn1a, +Krt18, +Ndrg1, +Il1r1]
Basal:
  markers: [Krt5, Trp63, Krt15, Krt17, +Krt14, +Ngfr]
Suprabasal / hillock:
  markers: [Krt4, Krt13, +Ecm1, +Lgals3, +Anxa1]
Goblet:
  markers: [Muc5ac, Muc5b, Spdef, Agr2, +Gcnt3, +Bpifb1]
Multiciliated:
  markers: [Foxj1, Tppp3, Rsph1, +Dynlrb2, +Cdhr3]
Deuterosomal:
  markers: [Deup1, Foxn4, Cdc20b, Ccno, +Mcidas]
PNEC:
  markers: [Ascl1, Calca, Chga, Chgb, Scg5, +Syp, +Resp18]
Ionocyte:
  markers: [Foxi1, Cftr, Ascl3, +Atp6v1b1, +Atp6v0a4]  # Atp6v1a is ubiquitous, drop as sole marker
Tuft:
  markers: [Pou2f3, Trpm5, Dclk1, Avil, Gnat3, +Il25, +Alox5ap]   # HES6, LRMP are weak in mouse

```

### Myeloid and mast

```yaml
Alveolar macrophage (AM):
  markers: [Siglecf, Itgax, Ear2, Chil3, Krt79, Pparg, Marco, Fabp1, Fabp4]
  note: ⚠ C1qa/C1qb are higher in interstitial macrophages in mouse; do not use as AM markers. Mrc1 is moderate. Siglecf alone is shared with eosinophils (add Ear2/Pparg, exclude Prg2/Epx).
Interstitial macrophage (IM):
  markers: [C1qa, C1qb, Lyve1, Folr2, Mrc1, Pf4, Cx3cr1, Selenop, Rnase1]   # CD163 weaker in mouse
  subsets: Lyve1-hi/Folr2-hi perivascular; MHCII-hi (H2-Aa, Cd74) Cx3cr1-hi
  note: exclude Lyve1+ lymphatic ECs.
Monocyte-derived macrophage (MoMac):
  markers: [Lyz2, Cd68, Csf1r, Mertk, Fcgr1, Ccr2, +Cx3cr1]
  note: In bleomycin, recruited Siglecf-low Itgax+ monocyte-derived AM (Misharin 2017) replace resident AM; do not infer origin from expression alone.
Classical monocyte:
  markers: [Ly6c2, Ccr2, Plac8, Vcan, Cd14, Sell]   # FCN1 not discriminating in mouse
Non-classical monocyte:
  markers: [Cx3cr1, Nr4a1, Pparg?, Cd36, Ace, Itgax-low, Treml4, Fcgr4]   # Fcgr3a → Fcgr3/Fcgr4
  note: exclude Fcgr3+ NK.
cDC1:
  markers: [Clec9a, Xcr1, Irf8, Wdfy4, Cadm1, +Batf3, +Itgae]
cDC2:
  markers: [Sirpa, Cd209a, Clec10a, Clec4a4, Fcer1a?, Irf4, Itgam]   # CD1C/CD1E no ortholog; Fcer1a is weak on mouse cDC2
Mature / migratory DC (mregDC):
  markers: [Ccr7, Fscn1, Lamp3, Cd83, Ccl22, Cd274, Pdcd1lg2, Cd200, Traf1]   # CCL19 weak in mouse DC
pDC:
  markers: [Siglech, Bst2, Ly6d, Tcf4, Pld4, Irf7, Ccr9, Gzmb]   # CLEC4C, LILRA4 no ortholog; IL3RA low
Neutrophil:
  markers: [Csf3r, Cxcr2, S100a8, S100a9, G0s2, Ly6g, Retnlg, Mmp9, Il1b]
  note: immature: Ltf, Camp, Ngp, Lcn2, Chil3.
Mast:
  markers: [Tpsab1, Tpsb2, Cpa3, Ms4a2, Mcpt4, Cma1, Kit, Gata2]
Basophil (mouse add):
  markers: [Mcpt8, Cd200r3, Ms4a2, Fcer1a]
Eosinophil:
  markers: [Siglecf, Prg2, Epx, Ear6, Ccr3, Il5ra, Clc]   # RNASE2 → Ear-family; SIGLEC8 → Siglecf (shared with AM)
```

### T, NK, ILC

```yaml
CD4 T:   {markers: [Cd3d, Trac, Cd4]}
CD8 T:   {markers: [Cd3d, Trac, Cd8a, Cd8b1]}     # CD8B → Cd8b1
Treg:    {markers: [Foxp3, Il2ra, Ctla4, Ikzf2, Tigit]}
NK:
  markers: [Klrd1, Ncr1, Nkg7, Prf1, Ctsw, Gzma, Klrb1c, Eomes, Klrk1]   # GNLY dropped
  note: Klrb1c (NK1.1) is strain-dependent (absent in BALB/c). Confirm no Cd3 program.
Mature NK (replaces CD16+ NK):
  markers: [Itgam, Cx3cr1, Klrg1, Zeb2, S1pr5, Fgfbp2, Fcgr3]
Immature / tissue NK (replaces CD56bright):
  markers: [Cd27, Xcl1, Sell, Kit?, Itga1, Cxcr6]   # NCAM1, XCL2 absent
MAIT:
  markers: [Trav1, Slc4a10, Zbtb16, Klrb1c, Il18rap, Rorc]   # TRAV1-2 → Trav1 (Vα19-Jα33)
  note: Rare in SPF C57BL/6 lung; conditional.
Gamma-delta T:
  markers: [Trdc, Trgc1, Trgc2, Trgc4, Tcrg-V gene usage]
  gdT17: [Il17a, Il17f, Rorc, Cd163l1 (Scart2), Blk, Sox13, Maf, Ccr6]
  note: lung γδT are mostly IL-17 type; require T-lineage.
ILC:
  ILC2: [Il1rl1, Il5, Il13, Gata3, Klrg1, Areg, Arg1, Icos]
  ILC3: [Rorc, Il22, Il23r, Kit, Ccr6]
  ILC1: [Tbx21, Ncr1, Il7r, Eomes-low]
  note: Il7r/Kit/Rora/Klrb1 shared with T and NK; separate by Cd3-negative and Il1rl1/Gata3.
```

### B and plasma

```yaml
Naive B:
  markers: [Fcer2a, Ighd, Ighm, Cd19, Ms4a1, Sell, Cd55?]   # TCL1A → Tcl1
Memory B:
  markers: [Cd80, Pdcd1lg2, Nt5e, Cd38, Ccr6, Zbtb32, class-switched Ighg/Igha]   # CD27 is not a mouse memory-B marker
Marginal-zone / B1 (mouse add):
  MZB: [Cr2, Cd1d1, Cd9, Cd1d1]   B1: [Cd5, Spn, Cd9, Itgam, high Ighm]
GC B:
  markers: [Rgs13, Lmo2, Bcl6, Aicda, Hmces, S1pr2, Gcsam, Fas]
Plasmablast:
  markers: [Prdm1, Xbp1, Mzb1, Jchain, Tnfrsf17, Mki67, Slamf7]
Plasma:
  markers: [Mzb1, Jchain, Xbp1, Sdc1, Tnfrsf17, Igkc, Iglc1/2, Igha, Ighg1]
pDC:  see myeloid pDC (Siglech, Bst2, Ly6d); do not use Jchain as the pDC marker in mouse.
```

### Endothelial

```yaml
aCap:
  markers: [Car4, Ednrb, Apln, Tbx2, Hpgd, Emcn, Sema3c]   # Car4/Apln/Tbx2 are the Gillich 2020 anchors
gCap:
  markers: [Plvap, Gpihbp1, Kit, Rgcc, Epas1, Aplnr]       # FCN3 no ortholog; Kit marks the progenitor-like gCap subset
Arterial:
  markers: [Gja5, Gja4, Bmx, Sox17, Efnb2, Hey1, Fbln5, Dkk2, Cxcl12, Igfbp3]
Venous:
  markers: [Ackr1, Vcam1, Selp, Nr2f2, Vwf, Slc6a2?]
Pulmonary venous:
  markers: [Hdac9?, Cpe?, Clu, Vcam1]  expected_low: [Col15a1]   # [verify] weaker mouse evidence
Systemic (bronchial) venous / systemic vessel:
  markers: [Col15a1, Vwa1, Plvap, Adamts9]   # Col15a1+ systemic vascular EC around airways in mouse (Gillich 2020)
Tip / angiogenic:
  markers: [Esm1, Kcne3, Angpt2, Dll4, Pgf, Cxcr4, Mcam]
  note: ⚠ Apln is an aCap marker in mouse; remove from the tip panel.
```

### Stromal

```yaml
Alveolar fibroblast:
  markers: [Npnt, Gpc3, Inmt, Col13a1, Fmo2?, Scube2, Adh1b?]   # Tsukui 2020; Npnt is the cleanest
Adventitial fibroblast:
  markers: [Pi16, Dpt, Mfap5, Cd34, Igfbp6, Sfrp2]
Peribronchial fibroblast:
  markers: [Hhip, Aspn, Fgf18, Actg2, Lgr6]    # Hhip, Aspn per Tsukui 2020; WIF1 less specific
Myofibroblast:
  markers: [Cthrc1, Postn, Tnc, Acta2, Col1a1, Ltbp2]   # Cthrc1+ is the bleomycin fibrotic population
Lipofibroblast (mouse add):
  markers: [Plin2, Lpl, Fabp4?, Tcf21, Fgf10, Ces1d]
SMC:
  markers: [Myh11, Acta2, Tagln, Cnn1, Lmod1, Myl9, Actg2 (airway SMC), Des]
Pericyte:
  markers: [Rgs5, Cspg4, Pdgfrb, Higd1b, Cox4i2, Kcnj8, Abcc9]
Mesothelial: {markers: [Msln, Upk3b, Calb2, Wt1, Lrrn4, Krt8, Krt18, Krt19]}
```

## L3

```yaml
Stress response:          [Fos, Jun, Junb, Egr1, Atf3, Dusp1]
IFN response:             [Isg15, Ifit1, Ifit3, Mx1, Oas1a, Rsad2, Irf7, Usp18, Ifi27l2a]   # OAS1 → Oas1a/g; IFI6 no ortholog
Proliferation:            [Mki67, Top2a, Ube2c, Stmn1]
Heat shock / dissociation stress (rename from "Hypoxia"): [Hspa1a, Hspa1b, Hspb1, Hsp90aa1, Dnajb1]
Hypoxia (true):           [Vegfa, Bnip3, Ndrg1, Slc2a1, Ldha, Pgk1, Egln3]
Metallothionein-high:     [Mt1, Mt2]      # MT1E/X/G, MT2A collapse to Mt1, Mt2

Macrophage SPP1-high:     [Spp1, Gpnmb, Trem2, Lgals3, Cd9, Fabp5, Mmp12, Chil1, Arg1?]   # CHI3L1 → Chil1; CHIT1, LILRB4 weak
Macrophage inflammatory 1: [Ccl2, Il1rn, Tnf]
Macrophage inflammatory 2: [Il1b, Cxcl1, Cxcl2, Cxcl3]    # CXCL8 no ortholog
EC IFN response:          [Isg20, Lap3, Cx3cl1, Irf1]

T naive:   [Ccr7, Sell, Lef1, Tcf7, Foxp1]
T central memory: [Ccr7, Sell, Il7r, Tcf7, Cd27]
TRM:       [Cd69, Itgae, Itga1, Cxcr6, Zfp683, Cd101]     # ZNF683 → Zfp683
CD8 Gzmk-high: [Gzmk, Cd27, Il7r, Ccl5]
CD8 cytotoxic: [Gzmb, Gzma, Prf1, Ifng, Cx3cr1, Tbx21]    # GNLY dropped
Terminal differentiation: [Klrg1, Klrc1, Cx3cr1, Zeb2, Klrk1]   # KLRC2/KLRC3 not usable in mouse

Fibroblast inflammatory: [Ccl2, Cxcl12, Cxcl14, Il6, Il33]
Resting:                 [Apod, Lepr]
Immune-recruiting:       [Cxcl13, Ccl19, Ccl21a, Il7]   # exclude Ccl21a+ LECs
```

## Key changes vs a naive symbol map

1. **Dropped, no ortholog:** GNLY, FCGR3B, CD1C, CD1E, CLEC4C, LILRA4, XCL2, FCN3, IFI6, CXCL8, KLRC3, NCAM1 (as NK discriminator), SFTPA2, TRAV1-2 (renamed Trav1).
2. **Renamed or collapsed:** LYZ→Lyz2, CCL21→Ccl21a, CYP2F1→Cyp2f2, CD8B→Cd8b1, CHI3L1→Chil1, OAS1→Oas1a, MT1*/MT2A→Mt1/Mt2, ZNF683→Zfp683, CD16 (FCGR3A)→Fcgr3/Fcgr4, TCL1A→Tcl1, FCER2→Fcer2a, CAP5→Caps, CLU1→Clu.
3. **Markers that behave differently:** C1qa/C1qb (IM-high, not AM), Lyve1 (macrophage and LEC), Apln (aCap, not tip), Siglecf (AM and eosinophil), Calca (PNEC and ILC2), Klrb1c (strain-dependent), CD27 (not memory B).
4. **Structures with no mouse equivalent:** AT0, TRB secretory, SMG (distal lung), human-style aberrant basaloid (replaced by the Krt8-high transitional state).
5. **Added from mouse literature:** Ly6g/Retnlg (neutrophil), Ly6c2/Ccr2 (monocyte), Siglech/Bst2/Ly6d (pDC), Il1rl1/Gata3 (ILC2), Car4/Tbx2/Gpihbp1/Kit (capillary), Hhip/Aspn/Dpt (fibroblast), Cthrc1 plus Krt8/Cldn4/Sfn state markers, lipofibroblast and basophil classes.

## Public references used

- Gillich A, et al. Capillary cell-type specialization in the alveolus. *Nature* 2020. PMID 33057196
- Strunz M, et al. Alveolar regeneration through a Krt8+ transitional stem cell state. *Nat Commun* 2020. PMID 32678092
- Choi J, et al. Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors. *Cell Stem Cell* 2020. PMID 32750316
- Kobayashi Y, et al. Persistence of a regeneration-associated, transitional alveolar epithelial cell state in pulmonary fibrosis. *Nat Cell Biol* 2020. PMID 32661339
- Tsukui T, et al. Collagen-producing lung cell atlas identifies multiple subsets with distinct localization and relevance to fibrosis. *Nat Commun* 2020. PMID 32317643
- Montoro DT, et al. A revised airway epithelial hierarchy includes CFTR-expressing ionocytes. *Nature* 2018. PMID 30069044
- Angelidis I, et al. An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics. *Nat Commun* 2019. PMID 30814501
- Misharin AV, et al. Monocyte-derived alveolar macrophages drive lung fibrosis and persist in the lung over the life span. *J Exp Med* 2017. PMID 28694385
- Zilionis R, et al. Single-cell transcriptomics of human and mouse lung cancers reveals conserved myeloid populations across individuals and species. *Immunity* 2019. PMID 30979687
- Schyns J, et al. Non-classical tissue monocytes and two functionally distinct populations of interstitial macrophages populate the mouse lung. *Nat Commun* 2019. PMID 31481690

Titles and the specific marker-to-study attributions were written from memory of these papers; only the PMIDs were checked against PubMed. Spot-check key genes in each paper's supplementary tables before treating them as final.
