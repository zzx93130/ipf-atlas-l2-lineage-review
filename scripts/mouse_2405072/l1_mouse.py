#!/usr/bin/env python
"""Mouse lung scRNA-seq (2405072: BLM, BLM_870, Saline) - L1 lineage analysis.

Stages (run in order, each as a fresh process):
  markers   build/validate the mouse L1 marker panel (from the lab's human L1 panel or a mouse CSV)
  qc        unzip 10x filtered matrices, input contract, MAD QC per sample
  doublets  Scrublet (+ scDblFinder, DoubletFinder via doublets_r.R if Rscript is available)
  embed     remove consensus doublets, lognorm, HVGs, regress/scale/PCA, Harmony, gates, UMAP, Leiden ladder
  evidence  per-cluster evidence at the L1 resolution, DE, proposals, review table
  figures   review package (cluster UMAP, proposed-L1 UMAP, sample UMAP, panel / top-DEG / state dot plots)
  recluster drop L1 clusters (--remove 3,7) and re-embed the rest from raw counts -> L1 round 2, 3, ...;
            evidence / figures / apply / split then use the latest round (or --round k)
  apply     apply the user's decisions (user_decisions_L1.csv): L1 labels for every cell, per-unit cell lists,
            state clusters (Proliferating, IFN response) routed to a lineage per cell
  split     one h5ad per lineage (raw counts, L1 metadata) under L1/lineages/ for L2

Example:
  python l1_mouse.py --stage markers --human-panel ROOT/atlas/L1/tables/L1_marker_panel.csv
  for s in qc doublets embed evidence figures; do python l1_mouse.py --stage $s || break; done
  # fill user_* columns in review/L1_review_table.xlsx (or .csv), then
  python l1_mouse.py --stage apply
"""
import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import date
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_mouse_L1")
sys.modules["tensorflow"] = None  # umap imports TF at load; on macOS arm64 that segfaults

import numpy as np
import pandas as pd
import scipy.sparse as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

# ----------------------------------------------------------------------------- config
DATA_DIR = Path("/Users/joe/SSD2_1T/IPF_ATLAS/data/mouse_data")
OUT_DIR = Path("/Users/joe/SSD2_1T/IPF_ATLAS/data/processed/mouse_2405072")
SCRIPT_DIR = Path(__file__).resolve().parent
SAMPLE_SHEET = SCRIPT_DIR / "sample_sheet.csv"
RSCRIPT = "Rscript"

UNITS = ["epithelial", "endothelial", "stromal", "myeloid", "tnk", "bplasma"]  # same units as the human atlas

# QC (sc-annotation-pipeline Stage 2: permissive MAD rule)
MAD_COUNTS, MAD_GENES, MAD_MT = 5, 5, 3
MT_HARD_CAP = 20.0        # % mito, hard cap on top of the 3-MAD rule (working value; log it)
MIN_GENES = 200           # absolute floor (working value; log it)
MAX_GENES = 8000          # absolute cap (working value; log it)
MAX_SAMPLE_LOSS = 0.40    # halt if a sample loses more than this at QC

# doublets: expected rate ~0.8% per 1,000 recovered cells (10x)
DBL_RATE_PER_1K = 0.008

# embedding
N_HVG = 3000
N_PCS = 50
N_NEIGHBORS = 15
LEIDEN_LADDER = [0.2, 0.5, 1.0]
L1_RES = 0.5              # human atlas L1 was reviewed at Leiden 0.5 (cluster ids "0.5-k")
LATENT = "auto"           # auto | harmony | pca

# marker scoring / calling (sc-annotation-pipeline Stage 8b)
AUCELL_TOP_FRAC = 0.05
SUPPORT_FLOOR = 0.20
DETECT_STRONG = 0.40
MARGIN_MIN = 0.75
DUAL_RATIO, DUAL_FRAC = 0.80, 0.60

# user's size-drop rule (operationalised; the user's rule has no number)
SIZE_DROP_RATIO = 0.5

# patterns (mouse)
MT_PAT = r"^mt-"
RIBO_PAT = r"^(Rps|Rpl|Mrps|Mrpl)\d"
RP_DE_PAT = r"^(Rps|Rpl)\d|^(Rps|Rpl)[0-9]+[a-z0-9-]*$|^Rplp|^Rpsa$|^mt-"
HB_PAT = r"^Hb[ab]-"
IGTR_VJ_PAT = r"^(Igh|Igk|Igl)[vj]|^Tr[abdg][vj]\d"

# state sets (human symbols + source; title-cased or ortholog-mapped at the markers stage)
STATE_DEFAULT = {
    "Proliferating": (["MKI67", "TYMS", "STMN1"], "lung_ipf_marker_panels.md: Proliferating myeloid"),
    "IFN response": (["ISG15", "MX1", "IFI6"], "lung_ipf_marker_panels.md: IFN-response AT2"),
    "Stress response": (["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1"], "sc-annotation-pipeline Stage 8c inducible chaperones"),
    "Metallothionein-high": (["MT1E", "MT1X"], "lung_ipf_marker_panels.md: AM metallothionein-high"),
}
STATE_FOLD = {"Proliferating": 2.5, "IFN response": 2.5, "Stress response": 2.0, "Metallothionein-high": 2.5}
STATE_FLOOR = 0.30

# keyword routing of the human panel's compartment names onto the units (checked at the markers stage)
UNIT_KEYWORDS = [
    ("epithelial", ["epithel", "at1", "at2", "alveolar type", "club", "ciliated", "basal", "secretory", "pnec", "ionocyte"]),
    ("endothelial", ["endothel", "lymphatic", "lec", "capillary", "vein", "arter"]),
    ("stromal", ["stroma", "fibro", "smooth", "smc", "pericyte", "mesothel", "mesenchym"]),
    ("myeloid", ["myeloid", "macrophage", "monocyte", "dendritic", "dc", "mast", "neutrophil", "granulocyte"]),
    ("tnk", ["t cell", "t/nk", "tnk", "nk", "ilc", "cd4", "cd8", "t_nk"]),
    ("bplasma", ["b cell", "b/plasma", "bplasma", "plasma", "b_plasma"]),
]


# ----------------------------------------------------------------------------- utilities
ROUND = 1  # L1 round the review / apply / split stages work on (set in main)


def emb_file(k):
    w = OUT_DIR / "work"
    return w / "03_embedded.h5ad" if k == 1 else w / f"03_embedded_round{k}.h5ad"


def latest_l1_round():
    ks = [k for k in range(1, 100) if emb_file(k).exists()]
    if not ks:
        raise SystemExit(f"no L1 embedding in {OUT_DIR / 'work'}: run --stage embed first")
    return max(ks)


def review_dir(k):
    return OUT_DIR / "L1" / "review" if k == 1 else OUT_DIR / "L1" / "review" / f"round{k}"


def paths():
    p = {
        "work": OUT_DIR / "work", "logs": OUT_DIR / "logs", "markers": OUT_DIR / "markers",
        "review": review_dir(ROUND), "figs": review_dir(ROUND) / "figures",
        "annot": OUT_DIR / "L1" / "annotation", "units": OUT_DIR / "L1" / "units", "final": OUT_DIR / "final",
    }
    for v in p.values():
        v.mkdir(parents=True, exist_ok=True)
    return p


def peak_rss_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e9 if sys.platform == "darwin" else r / 1e6


T0 = time.time()


def say(msg):
    print(f"[{time.strftime('%H:%M:%S')} +{time.time() - T0:6.0f}s rss {peak_rss_gb():5.1f}GB] {msg}", flush=True)


def _ver(pkg):
    from importlib.metadata import version, PackageNotFoundError
    try:
        return version(pkg)
    except PackageNotFoundError:
        return None


def log_event(stage, **kw):
    f = paths()["logs"] / "pipeline_log.json"
    log = json.loads(f.read_text()) if f.exists() else []
    log.append({"stage": stage, "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "versions": {k: _ver(k) for k in ("scanpy", "anndata", "harmonypy", "numpy")}, **kw})
    f.write_text(json.dumps(log, indent=1, default=str))


def read_samples():
    s = pd.read_csv(SAMPLE_SHEET, dtype=str)
    need = {"sample_id", "zip_path", "condition", "batch", "species"}
    miss = need - set(s.columns)
    if miss:
        raise SystemExit(f"sample sheet missing columns: {miss}")
    if s["sample_id"].duplicated().any():
        raise SystemExit("duplicate sample_id in sample sheet")
    if set(s["species"].str.lower()) != {"mouse"}:
        raise SystemExit("this script is for mouse only; do not mix species")
    return s


def markers_table():
    f = paths()["markers"] / "mouse_L1_markers.csv"
    if not f.exists():
        raise SystemExit(f"{f} missing: run --stage markers first")
    m = pd.read_csv(f, dtype=str).fillna("")
    return m


def marker_sets(m, var_names):
    """identity sets per unit and state sets, restricted to measured genes (genes missing from var reported)."""
    vn = pd.Index(var_names)
    ident, states, req, missing = {}, {}, {}, []
    for s, g in m.groupby("set"):
        genes = list(dict.fromkeys(g["gene"]))
        present = [x for x in genes if x in vn]
        missing += [(s, x) for x in genes if x not in vn]
        if s.startswith("state:"):
            if len(present) >= 2:
                states[s.split(":", 1)[1]] = present
        else:
            ident[s] = present
            r = [x for x in g.loc[g["required"].str.lower().isin(["1", "true", "yes"]), "gene"] if x in vn]
            if r:
                req[s] = r
    bad = [u for u in UNITS if len(ident.get(u, [])) < 2]
    if bad:
        raise SystemExit(f"units with <2 measured marker genes: {bad}; fix markers/mouse_L1_markers.csv")
    for s, r in req.items():  # invariant: a required gene must be a member of the set it gates
        assert set(r) <= set(ident[s]), f"required genes of {s} not in its set"
    return ident, states, req, missing


def select_hvgs(H):
    """seurat_v3_paper by sample; on a loess failure retry with a stricter gene filter + wider span, then fall back
    to the dispersion-based 'seurat' flavor on lognorm. Returns (hvg names, method dict)."""
    n_det = np.asarray((H.layers["counts"] > 0).sum(0)).ravel()
    for flavor, min_cells, span in [("seurat_v3_paper", 1, 0.3), ("seurat_v3_paper", max(10, int(0.005 * H.n_obs)), 0.5),
                                    ("seurat", max(10, int(0.005 * H.n_obs)), None)]:
        t = H[:, n_det >= min_cells].copy()
        try:
            if flavor == "seurat":
                t.X = t.layers["lognorm"]
                sc.pp.highly_variable_genes(t, n_top_genes=N_HVG, flavor=flavor, batch_key="sample_id")
            else:
                sc.pp.highly_variable_genes(t, n_top_genes=N_HVG, flavor=flavor, batch_key="sample_id",
                                            layer="counts", span=span)
        except ValueError as e:
            say(f"HVG {flavor} (genes in >={min_cells} cells, span {span}) failed: {str(e).strip()[:80]}")
            continue
        return t.var_names[t.var["highly_variable"]], {"flavor": flavor, "min_cells": int(min_cells), "span": span}
    raise SystemExit("all HVG attempts failed")


def regress_scale_pca(a, keys=("total_counts", "pct_counts_mt"), n_pcs=N_PCS, block=500, seed=0, n_jobs=8):
    """sc-annotation-pipeline Stage 4b: blocked regress_out -> scale(10) -> arpack PCA (kernel.regress_scale_pca)."""
    keys = list(keys)
    X = a.X
    n, g = X.shape
    cov = pd.DataFrame(a.obs[keys].to_numpy(np.float64), columns=keys, index=np.arange(n).astype(str))
    Xc = X.tocsc() if sp.issparse(X) else np.asarray(X)
    Y = np.empty((n, g), dtype=np.float32)
    for j0 in range(0, g, block):
        j1 = min(g, j0 + block)
        blk = Xc[:, j0:j1]
        blk = blk.toarray() if sp.issparse(blk) else np.array(blk)
        if keys:
            b = ad.AnnData(X=blk.astype(np.float32), obs=cov)
            sc.pp.regress_out(b, keys, n_jobs=n_jobs)
            blk = b.X
        Y[:, j0:j1] = blk
    t = ad.AnnData(X=Y)
    sc.pp.scale(t, max_value=10)
    nc = min(n_pcs, min(t.shape) - 1)
    sc.pp.pca(t, n_comps=nc, svd_solver="arpack", random_state=seed)
    return np.asarray(t.obsm["X_pca"], dtype=np.float32)


def qc_pc_correlation(Z, obs, n_top=10):
    out = {}
    for k in ["total_counts", "pct_counts_mt"]:
        v = obs[k].to_numpy(float)
        v = np.log1p(v) if k == "total_counts" else v
        out[k] = float(np.nanmax([abs(np.corrcoef(Z[:, i], v)[0, 1]) for i in range(min(n_top, Z.shape[1]))]))
    return out


def aucell(X, sets, top_frac=AUCELL_TOP_FRAC, chunk=2000):
    """AUCell-style rank score per cell (own gene ranking; independent of composition). X: cells x genes lognorm CSR.
    sets: dict name -> list of column indices. Undetected genes get no rank."""
    X = X.tocsr()
    n, g = X.shape
    n_up = max(1, int(round(g * top_frac)))
    names = list(sets)
    out = np.zeros((n, len(names)), dtype=np.float32)
    for i0 in range(0, n, chunk):
        blk = X[i0:i0 + chunk].toarray()
        m = blk.shape[0]
        top = np.argpartition(-blk, n_up - 1, axis=1)[:, :n_up]
        vals = np.take_along_axis(blk, top, 1)
        top = np.take_along_axis(top, np.argsort(-vals, axis=1, kind="stable"), 1)
        rank = np.full(blk.shape, n_up, dtype=np.int32)
        rank[np.arange(m)[:, None], top] = np.arange(n_up, dtype=np.int32)[None, :]
        rank[blk <= 0] = n_up
        for k, nm in enumerate(names):
            idx = sets[nm]
            out[i0:i0 + m, k] = np.clip(n_up - rank[:, idx], 0, None).sum(1) / (len(idx) * n_up)
    return pd.DataFrame(out, columns=names)


def knn_indices(Z, k):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Z)
    return nn.kneighbors(Z, return_distance=False)[:, 1:]


def ilisi(idx, labels):
    codes, uniq = pd.factorize(labels)
    nb = codes[idx]
    p = np.stack([(nb == c).mean(1) for c in range(len(uniq))], 1)
    return float(np.median(1.0 / (p ** 2).sum(1)))


def detect_frac(X, rows, cols):
    sub = X[rows][:, cols]
    return np.asarray((sub > 0).mean(0)).ravel()


def cluster_detection(X, groups, cols):
    """fraction of cells detecting each gene, per group -> DataFrame groups x genes."""
    g = pd.Series(groups).astype(str)
    out = {}
    Xs = X[:, cols].tocsc()
    for k in g.unique():
        rows = np.where(g.values == k)[0]
        out[k] = np.asarray((Xs[rows] > 0).mean(0)).ravel()
    return pd.DataFrame(out).T


def find_10x_dir(root):
    for d, _, files in os.walk(root):
        if "__MACOSX" in d:
            continue
        if any(f in files for f in ("matrix.mtx.gz", "matrix.mtx")):
            return Path(d)
    raise SystemExit(f"no matrix.mtx(.gz) found under {root}")


# ----------------------------------------------------------------------------- stage: markers

# --- mouse panel in the lab's markdown format (marker_panel_mouse_v2.md)
MD_L1_ROUTE = {  # L1 compartment (or Lymphoid sub-list) -> unit; None = pan-lineage, not scored at L1
    "Epithelial": "epithelial", "Blood endothelial": "endothelial", "Lymphatic endothelial": "endothelial",
    "Fibroblast": "stromal", "Mural": "stromal", "Mesothelial": "stromal",
    "Myeloid": "myeloid", "Mast": "myeloid",
    "Lymphoid/markers": None, "Lymphoid/T_lineage": "tnk", "Lymphoid/NK_lineage": "tnk",
    "Lymphoid/B_lineage": "bplasma",
}
MD_L2_ROUTE = {"Epithelial": "epithelial", "Myeloid and mast": "myeloid", "T, NK, ILC": "tnk",
               "B and plasma": "bplasma", "Endothelial": "endothelial", "Stromal": "stromal"}
MD_STATE_ROUTE = {"Proliferation": "Proliferating", "IFN response": "IFN response",
                  "Stress response": "Stress response", "Metallothionein-high": "Metallothionein-high"}
MIN_L1_GENES = 3


def _md_tokens(txt):
    """bracket lists on a line -> [(gene, uncertain)]; '+' = supporting marker (kept), '?' = uncertain (not scored),
    '-low'/'-hi' qualifiers and free text are skipped."""
    out = []
    for lst in re.findall(r"\[([^\]]*)\]", txt.split("#")[0]):
        for tok in lst.split(","):
            tok = tok.strip()
            mm = re.match(r"^\+?([A-Za-z][A-Za-z0-9.\-]*)(\??)(.*)$", tok)
            if not mm or re.search(r"-(low|hi|high)$", mm.group(1)):
                continue
            rest = mm.group(3).strip()
            if rest and not rest.startswith("("):
                continue
            out.append((mm.group(1), mm.group(2) == "?"))
    return out


def _md_blocks(text):
    """{(level, section): [yaml lines]} for each ```yaml block under ## L1 / ## L2 ### X / ## L3."""
    blocks, level, section, inside, buf = {}, None, None, False, []
    for line in text.splitlines():
        if line.startswith("## "):
            level, section = line[3:].strip(), None
        elif line.startswith("### "):
            section = line[4:].strip()
        elif line.strip().startswith("```"):
            if inside:
                blocks[(level, section)] = buf
                buf = []
            inside = not inside
        elif inside:
            buf.append(line)
    return blocks


def parse_marker_md(path):
    """-> rows for the L1 sets (+ supplements), state rows, and every (gene, unit) claim in L1+L2 (for the level-0 rule)."""
    text = Path(path).read_text()
    blocks = _md_blocks(text)
    l1 = blocks.get(("L1", None))
    if not l1:
        raise SystemExit(f"no ```yaml block under '## L1' in {path}")
    rows, claims, uncertain, top = [], [], [], None
    for line in l1:
        if not line.strip():
            continue
        if not line.startswith(" "):
            top = line.split(":")[0].strip()
            sub = "markers"
        else:
            sub = line.strip().split(":")[0].strip()
        key = top if top != "Lymphoid" else f"Lymphoid/{sub}"
        if key not in MD_L1_ROUTE:
            raise SystemExit(f"L1 compartment '{key}' in {path} has no unit route; add it to MD_L1_ROUTE")
        unit = MD_L1_ROUTE[key]
        for g, unc in _md_tokens(line):
            if unc:
                uncertain.append(g)
                continue
            if unit is None:
                rows.append({"gene": g, "set": "dropped:pan-lineage", "source": f"{Path(path).name} L1 {key}"})
                continue
            claims.append((g, unit))
            rows.append({"gene": g, "set": unit, "source": f"{Path(path).name} L1 {key}"})
    l2_first = {}  # unit -> [(type, [anchor genes in order])]
    for (lev, sec), lines in blocks.items():
        if lev != "L2":
            continue
        unit = MD_L2_ROUTE.get(sec)
        if unit is None:
            raise SystemExit(f"L2 section '{sec}' has no unit route; add it to MD_L2_ROUTE")
        typ = None
        for line in lines:
            if line and not line.startswith(" "):
                typ = line.split(":")[0].strip()
            toks = _md_tokens(line)
            claims += [(g, unit) for g, _ in toks]
            if line.strip().startswith("markers:") or (line and not line.startswith(" ") and "[" in line):
                anchors = [g for g, unc in toks if not unc]
                plus = set(re.findall(r"\+([A-Za-z][A-Za-z0-9.\-]*)", line))
                l2_first.setdefault(unit, []).append((typ, [g for g in anchors if g not in plus]))
    states = []
    for line in blocks.get(("L3", None), []):
        if ":" not in line:
            continue
        nm = line.split(":")[0].strip()
        if nm in MD_STATE_ROUTE:
            genes = [g for g, unc in _md_tokens(line) if not unc]
            states += [{"gene": g, "set": f"state:{MD_STATE_ROUTE[nm]}", "source": f"{Path(path).name} L3 {nm}"}
                       for g in genes]
    return rows, states, claims, l2_first, uncertain

def default_marker_md():
    """marker_panel_mouse_v2.md next to this script, else the repo copy (<repo>/references/)."""
    for f in [SCRIPT_DIR / "marker_panel_mouse_v2.md", SCRIPT_DIR.parent.parent / "references" / "marker_panel_mouse_v2.md"]:
        if f.exists():
            return str(f)
    raise SystemExit("marker_panel_mouse_v2.md not found next to the script or in ../../references; use --marker-md")


def stage_markers(args):
    p = paths()
    rows = []
    claims = []
    md = args.marker_md or (None if (args.mouse_panel or args.human_panel) else default_marker_md())
    if md:
        r1, st_rows, claims, l2_first, uncertain = parse_marker_md(md)
        units_of = pd.DataFrame(claims, columns=["gene", "unit"]).drop_duplicates()
        n_units = units_of.groupby("gene")["unit"].nunique()
        for r in r1:
            r.update(required="", ortholog_method="mouse_panel")
        rows = [r for r in r1 if r["set"] in UNITS]
        dropped_pan = [r["gene"] for r in r1 if r["set"] == "dropped:pan-lineage"]
        supp = {}
        for u in UNITS:  # L1 set too small -> add the first unique non-'+' anchor of each L2 type of that unit
            have = [r["gene"] for r in rows if r["set"] == u and n_units.get(r["gene"], 1) == 1]
            if len(have) >= MIN_L1_GENES:
                continue
            for typ, anchors in l2_first.get(u, []):
                g = next((x for x in anchors if n_units.get(x, 1) == 1 and x not in have), None)
                if g:
                    have.append(g)
                    supp.setdefault(u, []).append(f"{g} ({typ.split(' →')[0]})")
                    rows.append({"gene": g, "set": u, "source": f"{Path(md).name} L2 {typ.split(' →')[0]} (L1 supplement)",
                                 "required": "", "ortholog_method": "mouse_panel"})
        for r in st_rows:
            r.update(required="", ortholog_method="mouse_panel")
        rows += st_rows
        shutil.copy(md, p["markers"] / Path(md).name)
        src_note = f"markdown panel {md}; pan-lineage dropped {dropped_pan}; uncertain ('?') not scored {uncertain}; " \
                   f"L1 supplements {supp}"
        say(f"parsed {md}: pan-lineage dropped {dropped_pan}; L1 supplements from L2 anchors {supp}")
    elif args.mouse_panel:
        m = pd.read_csv(args.mouse_panel, dtype=str).fillna("")
        need = {"gene", "set", "source"}
        if need - set(m.columns):
            raise SystemExit(f"mouse panel needs columns {need}")
        if "required" not in m:
            m["required"] = ""
        if "ortholog_method" not in m:
            m["ortholog_method"] = "mouse_panel"
        rows = m[["gene", "set", "source", "required", "ortholog_method"]].to_dict("records")
        src_note = f"mouse panel {args.mouse_panel}"
    else:
        if not args.human_panel:
            raise SystemExit("give --human-panel (lab human L1 panel CSV) or --mouse-panel")
        h = pd.read_csv(args.human_panel, dtype=str).fillna("")
        gcol = next((c for c in h.columns if c.lower() in ("gene", "symbol", "gene_symbol", "marker", "genes")), None)
        scol = next((c for c in h.columns if c.lower() in ("compartment", "lineage", "l1", "l1_lineage", "unit",
                                                             "set", "cell_type", "celltype", "group")), None)
        if gcol is None or scol is None:
            raise SystemExit(f"cannot find gene / compartment columns in {list(h.columns)}; "
                             "pass --mouse-panel with columns gene,set,source instead")
        orth = None
        if args.orthologs:
            o = pd.read_csv(args.orthologs, dtype=str)
            hc = next(c for c in o.columns if "human" in c.lower())
            mc = next(c for c in o.columns if "mouse" in c.lower())
            o = o[[hc, mc]].dropna()
            one2one = o.groupby(hc)[mc].nunique() == 1
            orth = o[o[hc].isin(one2one[one2one].index)].drop_duplicates(hc).set_index(hc)[mc]

        def unit_of(name):
            n = name.lower()
            hits = [u for u, kws in UNIT_KEYWORDS if any(k in n for k in kws)]
            return hits[0] if len(set(hits)) == 1 else ("AMBIGUOUS:" + "|".join(sorted(set(hits))) if hits else "")

        comp = pd.DataFrame({"compartment": sorted(h[scol].unique())})
        comp["unit"] = comp["compartment"].map(unit_of)
        comp.to_csv(p["markers"] / "human_compartment_to_unit.csv", index=False)
        bad = comp[(comp["unit"] == "") | comp["unit"].str.startswith("AMBIGUOUS")]
        if len(bad) and not args.force:
            raise SystemExit("compartments not routed to a unit (edit UNIT_KEYWORDS or use --force to drop them):\n"
                             + bad.to_string(index=False))
        cmap = comp.set_index("compartment")["unit"]
        for _, r in h.iterrows():
            u = cmap.get(r[scol], "")
            if not u or u.startswith("AMBIGUOUS"):
                continue
            hg = r[gcol].strip()
            if orth is not None:
                mg, how = orth.get(hg, ""), "ortholog_table_1to1"
            else:
                mg, how = hg.capitalize(), "titlecase"
            if mg:
                req = r.get("required", "") if "required" in h.columns else ""
                rows.append({"gene": mg, "set": u, "source": f"human L1 panel ({Path(args.human_panel).name}): "
                             f"{r[scol]} / {hg}", "required": req, "ortholog_method": how})
        src_note = f"human panel {args.human_panel}"
    for st, (genes, src) in (STATE_DEFAULT.items() if not md else []):
        for hg in genes:
            rows.append({"gene": hg.capitalize(), "set": f"state:{st}", "source": src, "required": "",
                         "ortholog_method": "titlecase"})
    m = pd.DataFrame(rows).drop_duplicates(["gene", "set"])
    # a gene claimed by two units is uninformative at L1 (level-0 rule): drop and report
    ident = m[~m["set"].str.startswith("state:")]
    cl_ = pd.concat([ident[["gene", "set"]], pd.DataFrame(claims, columns=["gene", "set"])]).drop_duplicates()
    multi = cl_.groupby("gene")["set"].nunique()
    multi = multi[multi > 1].index
    if len(multi):
        cl_[cl_["gene"].isin(multi)].sort_values("gene").to_csv(p["markers"] / "multi_unit_claims.csv", index=False)
    if len(multi):
        ident[ident["gene"].isin(multi)].to_csv(p["markers"] / "dropped_multi_unit_genes.csv", index=False)
    m = m[~(m["gene"].isin(multi) & ~m["set"].str.startswith("state:"))]
    m.to_csv(p["markers"] / "mouse_L1_markers.csv", index=False)
    summ = m.groupby("set")["gene"].apply(lambda x: ", ".join(x)).to_dict()
    say(f"marker panel written ({len(m)} rows). Dropped multi-unit genes: {list(multi)}")
    for k, v in summ.items():
        print(f"  {k}: {v}")
    log_event("markers", source=src_note, n_rows=len(m), dropped_multi_unit=list(multi), sets=summ)


# ----------------------------------------------------------------------------- stage: qc
def stage_qc(args):
    p = paths()
    s = read_samples()
    adatas, ret, gate_fail = [], [], []
    for _, r in s.iterrows():
        z = DATA_DIR / r["zip_path"]
        dest = p["work"] / "raw10x" / r["sample_id"]
        if not dest.exists():
            say(f"unzip {z}")
            with zipfile.ZipFile(z) as zf:
                zf.extractall(dest)
        d = find_10x_dir(dest)
        a = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True, gex_only=True)
        a.var["gene_symbols"] = a.var_names
        # Stage 0 contract
        assert not a.obs_names.duplicated().any(), f"duplicate barcodes in {r['sample_id']}"
        x = a.X.data
        assert np.allclose(x, np.round(x)) and (x >= 0).all(), f"{r['sample_id']}: X is not raw integer counts"
        n_mt_lower = a.var_names.str.match(MT_PAT).sum()
        n_mt_upper = a.var_names.str.match(r"^MT-").sum()
        if n_mt_lower == 0 or n_mt_upper > 0:
            raise SystemExit(f"{r['sample_id']}: species check failed (mt- genes {n_mt_lower}, MT- genes {n_mt_upper})")
        a.obs["barcode"] = a.obs_names
        a.obs_names = r["sample_id"] + "|" + a.obs_names
        for c in ["sample_id", "condition", "batch"]:
            a.obs[c] = r[c]
        a.var["mt"] = a.var_names.str.match(MT_PAT)
        a.var["ribo"] = a.var_names.str.match(RIBO_PAT)
        a.var["hb"] = a.var_names.str.match(HB_PAT)
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True)
        o = a.obs

        def mad_out(v, nm, lo=True):
            med = np.median(v)
            mad = np.median(np.abs(v - med)) * 1.4826
            return ((v < med - nm * mad) if lo else False) | (v > med + nm * mad), med - nm * mad, med + nm * mad

        f_c, c_lo, c_hi = mad_out(o["log1p_total_counts"].values, MAD_COUNTS)
        f_g, g_lo, g_hi = mad_out(o["log1p_n_genes_by_counts"].values, MAD_GENES)
        med = np.median(o["pct_counts_mt"])
        mt_hi = min(MT_HARD_CAP, med + MAD_MT * np.median(np.abs(o["pct_counts_mt"] - med)) * 1.4826)
        f_m = o["pct_counts_mt"].values > mt_hi
        f_n = (o["n_genes_by_counts"].values < MIN_GENES) | (o["n_genes_by_counts"].values > MAX_GENES)
        o["qc_fail_counts"], o["qc_fail_genes"], o["qc_fail_mt"], o["qc_fail_mingenes"] = f_c, f_g, f_m, f_n
        o["qc_pass"] = ~(f_c | f_g | f_m | f_n)
        loss = 1 - o["qc_pass"].mean()
        ret.append({"sample_id": r["sample_id"], "condition": r["condition"], "n_input": a.n_obs,
                    "n_pass": int(o["qc_pass"].sum()), "frac_lost": round(loss, 4),
                    "fail_counts": int(f_c.sum()), "fail_genes": int(f_g.sum()), "fail_mt": int(f_m.sum()),
                    "fail_min_max_genes": int(f_n.sum()),
                    "counts_range": f"{np.expm1(c_lo):.0f}-{np.expm1(c_hi):.0f}",
                    "genes_range": f"{max(np.expm1(g_lo), MIN_GENES):.0f}-{min(np.expm1(g_hi), MAX_GENES):.0f}",
                    "mt_max_pct": round(float(mt_hi), 2),
                    "median_counts": float(o["total_counts"].median()),
                    "median_genes": float(o["n_genes_by_counts"].median()),
                    "median_pct_mt": float(med)})
        if loss > MAX_SAMPLE_LOSS:
            gate_fail.append(r["sample_id"])
        adatas.append(a)
        say(f"{r['sample_id']}: {a.n_obs} cells, {a.n_vars} genes, pass {o['qc_pass'].sum()} ({loss:.1%} lost)")
    ids = [tuple(x.var["gene_ids"]) for x in adatas]
    if len(set(ids)) != 1:
        raise SystemExit("feature lists differ between samples (different references?) - stop and check")
    A = ad.concat(adatas, join="outer", merge="same")
    A.var = adatas[0].var[["gene_ids", "feature_types", "gene_symbols", "mt", "ribo", "hb"]].copy()
    A.layers["counts"] = A.X.copy()
    rt = pd.DataFrame(ret)
    rt.to_csv(p["review"] / "QC_retention.csv", index=False)
    print(rt.to_string(index=False))
    # before/after violins
    fig, axs = plt.subplots(2, 3, figsize=(13, 7))
    for j, (k, lab) in enumerate([("total_counts", "UMIs"), ("n_genes_by_counts", "genes"),
                                  ("pct_counts_mt", "% mito")]):
        for i, sel in enumerate([np.ones(A.n_obs, bool), A.obs["qc_pass"].values]):
            ax = axs[i, j]
            data = [A.obs.loc[sel & (A.obs["sample_id"] == sid).values, k].values for sid in s["sample_id"]]
            ax.violinplot(data, showmedians=True)
            ax.set_xticks(range(1, len(s) + 1), s["sample_id"], rotation=20, fontsize=8)
            if k != "pct_counts_mt":
                ax.set_yscale("log")
            ax.set_title(f"{lab} - {'before' if i == 0 else 'after'} QC", fontsize=10)
    fig.tight_layout()
    fig.savefig(p["figs"] / "QC_before_after.png", dpi=150)
    plt.close(fig)
    A.write_h5ad(p["work"] / "01_qc.h5ad", compression="gzip")
    log_event("qc", retention=ret, rules={"mad_counts": MAD_COUNTS, "mad_genes": MAD_GENES, "mad_mt": MAD_MT,
                                          "mt_hard_cap": MT_HARD_CAP, "min_genes": MIN_GENES, "max_genes": MAX_GENES},
              ambient="SKIPPED: only filtered matrices supplied; SoupX/CellBender need raw_feature_bc_matrix",
              gate_failed=gate_fail)
    if gate_fail and not args.force:
        raise SystemExit(f"QC gate: samples losing >{MAX_SAMPLE_LOSS:.0%}: {gate_fail}. Inspect QC_before_after.png; "
                         "rerun with --force to continue")


# ----------------------------------------------------------------------------- stage: doublets
def stage_doublets(args):
    p = paths()
    A = sc.read_h5ad(p["work"] / "01_qc.h5ad")
    A = A[A.obs["qc_pass"].values].copy()
    A.X = A.layers["counts"]
    res = []
    rscript = shutil.which(RSCRIPT)
    for sid in A.obs["sample_id"].unique():
        a = A[A.obs["sample_id"] == sid].copy()
        rate = DBL_RATE_PER_1K * a.n_obs / 1000
        sc.pp.scrublet(a, expected_doublet_rate=rate, random_state=0)
        df = pd.DataFrame({"cell_id": a.obs_names, "scrublet_score": a.obs["doublet_score"].values,
                           "scrublet_call": a.obs["predicted_doublet"].astype(bool).values})
        df["scdbl_score"], df["scdbl_call"], df["df_pann"], df["df_call"] = np.nan, np.nan, np.nan, np.nan
        if rscript:
            bc = p["work"] / "doublets" / f"{sid}_barcodes.txt"
            bc.parent.mkdir(exist_ok=True)
            pd.Series(a.obs["barcode"].values).to_csv(bc, index=False, header=False)
            out = p["work"] / "doublets" / f"{sid}_R.csv"
            tenx = find_10x_dir(p["work"] / "raw10x" / sid)
            cmd = [rscript, str(SCRIPT_DIR / "doublets_r.R"), str(tenx), str(bc), f"{rate:.5f}", str(out)]
            say(" ".join(cmd))
            rr = subprocess.run(cmd, capture_output=True, text=True)
            (p["logs"] / f"doublets_r_{sid}.log").write_text(rr.stdout + "\n" + rr.stderr)
            if rr.returncode == 0 and out.exists():
                r = pd.read_csv(out)
                r["cell_id"] = sid + "|" + r["barcode"]
                df = df.drop(columns=["scdbl_score", "scdbl_call", "df_pann", "df_call"]).merge(
                    r.drop(columns="barcode"), on="cell_id", how="left")
            else:
                say(f"R doublet methods failed for {sid}; see logs/doublets_r_{sid}.log")
        df["sample_id"], df["expected_rate"] = sid, rate
        res.append(df)
        say(f"{sid}: scrublet calls {int(df['scrublet_call'].sum())} / {len(df)} (expected rate {rate:.3f})")
    D = pd.concat(res, ignore_index=True)
    meth = [c for c in ["scrublet_call", "scdbl_call", "df_call"] if D[c].notna().any()]
    for c in meth:
        D[c] = D[c].astype("boolean").fillna(False).astype(bool)
    D["n_doublet_methods"] = len(meth)
    D["n_doublet_calls"] = D[meth].sum(axis=1).astype(int)
    need = 2 if len(meth) >= 2 else 99   # >=2 of 3 (or 2 of 2); one method -> flag only, no removal
    D["doublet_consensus"] = D["n_doublet_calls"] >= need
    D["doublet_candidate"] = (D["n_doublet_calls"] >= 1) & ~D["doublet_consensus"]
    D.to_parquet(p["work"] / "02_doublets.parquet", index=False)
    summ = D.groupby("sample_id").agg(n=("cell_id", "size"), consensus=("doublet_consensus", "sum"),
                                      candidate=("doublet_candidate", "sum")).reset_index()
    print(summ.to_string(index=False))
    log_event("doublets", methods=meth, consensus_rule=f">={need} of {len(meth)}", per_sample=summ.to_dict("records"),
              deviation=None if len(meth) == 3 else f"only {meth} ran; 3-method consensus not available")


# ----------------------------------------------------------------------------- stage: embed
def stage_embed(args):
    p = paths()
    m = markers_table()
    A = sc.read_h5ad(p["work"] / "01_qc.h5ad")
    A = A[A.obs["qc_pass"].values].copy()
    D = pd.read_parquet(p["work"] / "02_doublets.parquet").set_index("cell_id")
    A.obs = A.obs.join(D.drop(columns=["sample_id"]), how="left")
    n0 = A.n_obs
    A = A[~A.obs["doublet_consensus"].astype(bool).values].copy()
    say(f"removed {n0 - A.n_obs} consensus doublets; {A.n_obs} cells")
    keep_g = np.asarray((A.layers["counts"] > 0).sum(0)).ravel() > 0
    A = A[:, keep_g].copy()
    embed_cells(A, m, 1)


def embed_cells(A, m, rnd, extra_log=None):
    """lognorm -> HVGs -> regress/scale/PCA -> Harmony -> AUCell -> gates -> kNN/UMAP/Leiden ladder; writes emb_file(rnd).
    A: cells x all genes with layers['counts'] (raw) and the QC / doublet obs columns."""
    A.X = A.layers["counts"].copy()
    sc.pp.normalize_total(A)
    sc.pp.log1p(A)
    A.layers["lognorm"] = A.X.copy()
    # HVGs on counts, excluding mt/ribo/hb and IG/TR V/J genes
    excl = (A.var_names.str.match(MT_PAT) | A.var_names.str.match(RIBO_PAT) | A.var_names.str.match(HB_PAT)
            | A.var_names.str.match(IGTR_VJ_PAT))
    hvg, hvg_method = select_hvgs(A[:, ~excl].copy())
    A.var["highly_variable"] = A.var_names.isin(hvg)
    Ah = A[:, A.var["highly_variable"].values].copy()
    Ah.X = Ah.layers["lognorm"]
    say(f"HVGs: {len(hvg)} ({hvg_method}); regress/scale/PCA")
    Z0 = regress_scale_pca(Ah, keys=[])
    Z = regress_scale_pca(Ah)
    corr = {"unregressed_control": qc_pc_correlation(Z0, A.obs), "regressed": qc_pc_correlation(Z, A.obs)}
    A.obsm["X_pca_unregressed"], A.obsm["X_pca"] = Z0, Z
    import harmonypy
    ho = harmonypy.run_harmony(Z, A.obs[["sample_id"]], ["sample_id"], random_state=0)
    Zh = np.asarray(ho.Z_corr)
    Zh = Zh.T if Zh.shape[0] != A.n_obs else Zh
    A.obsm["X_pca_harmony"] = Zh.astype(np.float32)
    # per-cell AUCell over the full gene space (computed once; used by the gates and the evidence stage)
    ident, states, req, missing = marker_sets(m, A.var_names)
    if missing:
        say(f"marker genes not in the mouse matrix (check orthologs): {missing}")
    col = {g: i for i, g in enumerate(A.var_names)}
    S = aucell(A.layers["lognorm"], {u: [col[g] for g in ident[u]] for u in UNITS})
    S.index = A.obs_names
    for u in UNITS:
        A.obs[f"auc_{u}"] = S[u].values
    zc = (S - S.mean()) / S.std(ddof=0)
    srt = np.sort(zc.values, 1)
    A.obs["cell_unit"] = zc.idxmax(axis=1).values
    A.obs["cell_unit_margin"] = srt[:, -1] - srt[:, -2]
    conf = A.obs["cell_unit_margin"].values >= 1.0
    # gates: G1 iLISI(sample) ; G4 known-pair kNN cross-fraction among confident cells
    pairs = [("bplasma", "tnk"), ("myeloid", "tnk"), ("endothelial", "stromal"), ("epithelial", "stromal")]
    gates = {}
    for nm, Zx in [("pca", Z), ("harmony", Zh)]:
        idx = knn_indices(Zx[:, :N_PCS], 30)
        lab = A.obs["cell_unit"].values
        g4 = {}
        for a_, b_ in pairs:
            va = np.where(conf & (lab == a_))[0]
            vb = np.where(conf & (lab == b_))[0]
            if len(va) < 20 or len(vb) < 20:
                g4[f"{a_}-{b_}"] = np.nan
                continue
            fa = (np.isin(idx[va], vb)).mean()
            fb = (np.isin(idx[vb], va)).mean()
            g4[f"{a_}-{b_}"] = float(max(fa, fb))
        gates[nm] = {"iLISI_sample_median": ilisi(idx, A.obs["sample_id"].values), "G4_pair_crossfrac": g4,
                     "G4_worst": float(np.nanmax(list(g4.values()))) if any(np.isfinite(list(g4.values()))) else np.nan}
    g1 = gates["harmony"]["iLISI_sample_median"] > gates["pca"]["iLISI_sample_median"]
    g4ok = gates["harmony"]["G4_worst"] <= max(0.05, 1.5 * gates["pca"]["G4_worst"])
    choice = LATENT if LATENT != "auto" else ("harmony" if (g1 and g4ok) else "pca")
    why = f"G1 {'pass' if g1 else 'fail'}, G4 {'pass' if g4ok else 'fail'}" if LATENT == "auto" else "set in config"
    rep = {"harmony": "X_pca_harmony", "pca": "X_pca"}[choice]
    say(f"gates: {json.dumps(gates, default=float)} -> latent {choice} ({why})")
    sc.pp.neighbors(A, n_neighbors=N_NEIGHBORS, use_rep=rep, n_pcs=N_PCS, random_state=0)
    sc.tl.umap(A, random_state=0)
    for r in LEIDEN_LADDER:
        sc.tl.leiden(A, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, directed=False,
                     random_state=0)
        say(f"leiden {r}: {A.obs[f'leiden_{r}'].nunique()} clusters")
    A.uns["L1_latent"] = {"choice": choice, "why": why}
    A.uns["L1_round"] = rnd
    A.X = A.layers["counts"]
    A.write_h5ad(emb_file(rnd), compression="gzip")
    log_event("embed", round=rnd, **(extra_log or {}), n_cells=A.n_obs, n_genes=A.n_vars, n_hvg=len(hvg),
              hvg_method=hvg_method, pc_qc_corr=corr, gates=gates,
              latent=choice, latent_why=why, missing_marker_genes=missing,
              leiden={str(r): int(A.obs[f"leiden_{r}"].nunique()) for r in LEIDEN_LADDER},
              note="sample_id == condition here (one sample per condition): Harmony on sample_id can also "
                   "remove condition signal; compare L1 calls between latents before trusting either")


# ----------------------------------------------------------------------------- stage: recluster
def stage_recluster(args):
    """Drop L1 clusters of a round (Leiden L1_RES ids) and re-embed the remaining cells from raw counts."""
    global ROUND
    ids = [x.strip() for x in (args.remove or "").split(",") if x.strip()]
    if not ids:
        raise SystemExit("--stage recluster needs --remove <cluster ids>, e.g. --remove 3,7")
    prev = args.from_round or latest_l1_round()
    new = prev + 1
    key = f"leiden_{args.from_res or L1_RES}"
    P = sc.read_h5ad(emb_file(prev))
    if key not in P.obs:
        raise SystemExit(f"{key} not in L1 round {prev}; available {[c for c in P.obs if c.startswith('leiden_')]}")
    pc = P.obs[key].astype(str)
    bad = [c for c in ids if c not in set(pc)]
    if bad:
        raise SystemExit(f"clusters {bad} not in L1 round {prev} {key} (has {sorted(set(pc), key=int)})")
    stamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    for old in [emb_file(new), review_dir(new)]:  # never delete an earlier attempt
        if old.exists():
            old.rename(old.with_name(f"{old.stem}_prev_{stamp}{old.suffix}"))
    rm = pc.isin(ids).values
    removed = pd.DataFrame({"cell_id": P.obs_names[rm], "sample_id": P.obs["sample_id"].astype(str).values[rm],
                            "from_round": prev, "cluster_key": key, "cluster": pc.values[rm],
                            "reason": [f"L1 round{prev} {key} cluster {c}: {args.reason}" for c in pc.values[rm]]})
    removed.to_csv(OUT_DIR / "work" / f"L1_round{new}_removed_cells.csv", index=False)
    say(f"L1 round {new}: removing clusters {ids} of round {prev} {key} ({rm.sum():,} cells "
        f"{removed['cluster'].value_counts().to_dict()}); {(~rm).sum():,} cells re-embedded")
    drop = [c for c in P.obs.columns if c.startswith(("leiden_", "auc_")) or c in ("cell_unit", "cell_unit_margin")]
    obs = P.obs.loc[~rm].drop(columns=drop)
    prev_col = f"round{prev}_{key}"
    obs[prev_col] = pc.values[~rm]
    A = ad.AnnData(X=P.layers["counts"][~rm].copy(), obs=obs, var=P.var[[c for c in P.var.columns
                                                                         if c != "highly_variable"]].copy())
    A.layers["counts"] = A.X.copy()
    del P
    ROUND = new
    embed_cells(A, markers_table(), new, extra_log={"from_round": prev, "removed_clusters": ids,
                                                     "n_removed": int(rm.sum()), "reason": args.reason})
    A = sc.read_h5ad(emb_file(new))
    ct = pd.crosstab(A.obs[prev_col], A.obs[f"leiden_{L1_RES}"])
    review_dir(new).mkdir(parents=True, exist_ok=True)
    ct.to_csv(review_dir(new) / "L1_prev_vs_new_clusters.csv")
    say(f"previous clusters ({prev_col}, rows) x new clusters (leiden_{L1_RES}, columns):\n{ct.to_string()}")
    say(f"next: --stage evidence, --stage figures (round {new} is now the latest)")


# ----------------------------------------------------------------------------- stage: evidence
def stage_evidence(args):
    p = paths()
    m = markers_table()
    A = sc.read_h5ad(emb_file(ROUND))
    key = f"leiden_{L1_RES}"
    A.obs[key] = A.obs[key].astype(str)
    cl = A.obs[key]
    ident, states, req, missing = marker_sets(m, A.var_names)
    X = A.layers["lognorm"].tocsr()
    C = A.layers["counts"].tocsr()
    col = {g: i for i, g in enumerate(A.var_names)}
    gene2unit = {g: u for u in UNITS for g in ident[u]}

    # --- DE: full-data Wilcoxon cluster-vs-rest, ranked by scores
    A.X = X
    say("Wilcoxon DE")
    sc.tl.rank_genes_groups(A, key, method="wilcoxon", use_raw=False, pts=True, key_added="de")
    de = sc.get.rank_genes_groups_df(A, group=None, key="de")
    de = de.rename(columns={"names": "gene", "logfoldchanges": "logFC", "pvals_adj": "padj",
                            "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
    de["group"] = de["group"].astype(str)
    up = de[(de["padj"] < 0.05) & (de["logFC"] > 0)].sort_values(["group", "scores"], ascending=[True, False])
    rp = up["gene"].str.match(RP_DE_PAT)
    raw_top10 = up.groupby("group").head(10)
    rp_frac = raw_top10.assign(rp=rp.loc[raw_top10.index]).groupby("group")["rp"].mean()
    up = up[~rp]
    top30 = up.groupby("group").head(30).copy()
    top30["rank"] = top30.groupby("group").cumcount() + 1
    top30["marker_unit"] = top30["gene"].map(gene2unit).fillna("")
    top30[["group", "rank", "gene", "scores", "logFC", "padj", "pct_in", "pct_out", "marker_unit"]].to_csv(
        p["review"] / "L1_DE_top30.csv", index=False)

    # --- cluster-level marker scores / calls
    auc = A.obs[[f"auc_{u}" for u in UNITS]].copy()
    auc.columns = UNITS
    cm = auc.groupby(cl.values).mean()
    z = (cm - cm.mean()) / cm.std(ddof=0).replace(0, 1)
    det = {u: cluster_detection(C, cl.values, [col[g] for g in ident[u]]) for u in UNITS}
    for u in UNITS:
        det[u].columns = ident[u]
    clusters = sorted(cl.unique(), key=lambda x: (-int((cl == x).sum()), x))
    rows = []
    for c in clusters:
        order = z.loc[c].sort_values(ascending=False)
        call, path, veto = "unresolved", "", []
        for u in order.index:
            best = det[u].loc[c].max()
            if u in req and det[u].loc[c, req[u]].max() < SUPPORT_FLOOR:
                veto.append(f"{u} (required {'/'.join(req[u])} < {SUPPORT_FLOOR:.0%})")
                continue
            if best < SUPPORT_FLOOR:
                continue
            nxt = [v for v in order.index if v != u]
            margin = float(order[u] - order[nxt].max())
            if best >= DETECT_STRONG:
                call, path = u, "detection"
            elif margin >= MARGIN_MIN:
                call, path = u, "margin"
            break
        u1, u2 = order.index[0], order.index[1]
        margin12 = float(order.iloc[0] - order.iloc[1])
        rows.append({"cluster": c, "call": call, "call_path": path, "top_set": u1, "runner_up": u2,
                     "margin": round(margin12, 3), "veto": "; ".join(veto),
                     "best_gene": det[u1].loc[c].idxmax(), "best_gene_det": round(float(det[u1].loc[c].max()), 3),
                     "set_detection": ", ".join(f"{g} {v:.0%}" for g, v in det[u1].loc[c].sort_values(
                         ascending=False).head(5).items())})
    T = pd.DataFrame(rows).set_index("cluster")

    # --- size, composition, QC
    n = cl.value_counts()
    T["n_cells"] = n.reindex(T.index).values
    T["size_rank"] = np.arange(1, len(T) + 1)
    T["ratio_to_prev"] = np.r_[np.nan, T["n_cells"].values[1:] / T["n_cells"].values[:-1]]
    drop_pos = np.where(T["ratio_to_prev"].values < SIZE_DROP_RATIO)[0]
    drop_at = int(drop_pos[0]) if len(drop_pos) else len(T)
    T["below_size_drop"] = np.arange(len(T)) >= drop_at
    comp = pd.crosstab(cl, A.obs["sample_id"], normalize="index")
    for s_ in comp.columns:
        T[f"frac_{s_}"] = comp.reindex(T.index)[s_].round(3).values
    T["top_sample"] = comp.idxmax(axis=1).reindex(T.index).values
    T["top_sample_share"] = comp.max(axis=1).reindex(T.index).round(3).values
    g = A.obs.groupby(cl.values)
    T["median_genes"] = g["n_genes_by_counts"].median().reindex(T.index).values
    T["median_counts"] = g["total_counts"].median().reindex(T.index).values
    T["median_pct_mt"] = g["pct_counts_mt"].median().reindex(T.index).round(2).values
    T["scrublet_median"] = g["scrublet_score"].median().reindex(T.index).round(3).values
    T["frac_doublet_candidate"] = g["doublet_candidate"].mean().reindex(T.index).round(3).values
    T["frac_low_info_top10"] = rp_frac.reindex(T.index).fillna(0).round(2).values

    # --- per-cell dual-lineage test (>=2 detected markers of each of two units)
    cnt = pd.DataFrame({u: np.asarray((C[:, [col[x] for x in ident[u]]] > 0).sum(1)).ravel() for u in UNITS},
                       index=A.obs_names)
    has2 = cnt >= 2
    dual_ratio, dual_frac, second = [], [], []
    for c in T.index:
        rr = has2[cl.values == c]
        prev = rr.mean()
        own = T.loc[c, "call"] if T.loc[c, "call"] in UNITS else T.loc[c, "top_set"]
        oth = prev.drop(own).idxmax()
        second.append(oth)
        dual_ratio.append(float(prev[oth] / max(prev[own], 1e-9)))
        dual_frac.append(float((rr[own] & rr[oth]).mean()))
    T["second_lineage"], T["dual_ratio"], T["dual_frac"] = second, np.round(dual_ratio, 3), np.round(dual_frac, 3)
    dual_med = float(np.median(dual_frac))
    T["dual_flag"] = (T["dual_ratio"] >= DUAL_RATIO) & (T["dual_frac"] >= max(DUAL_FRAC, 2 * dual_med))

    # --- DE-based contamination / own-identity genes
    def de_counts(c, k):
        t = top30[top30["group"] == c].head(k)
        own = T.loc[c, "call"] if T.loc[c, "call"] in UNITS else T.loc[c, "top_set"]
        mu = t["marker_unit"]
        return int((mu == own).sum()), int(((mu != "") & (mu != own)).sum()), \
            ", ".join(f"{gg}({uu})" for gg, uu in zip(t["gene"], mu) if uu and uu != own)
    for k in (10, 30):
        res = [de_counts(c, k) for c in T.index]
        T[f"own_markers_top{k}"] = [r[0] for r in res]
        T[f"other_markers_top{k}"] = [r[1] for r in res]
        if k == 30:
            T["other_marker_genes_top30"] = [r[2] for r in res]
    T["top10_DE"] = [", ".join(f"{r.gene} {r.pct_in:.0%}/{r.pct_out:.0%}" for r in
                               top30[top30["group"] == c].head(10).itertuples()) for c in T.index]

    # --- complexity vs same-call clusters (ambient vs doublet diagnosis)
    clean = ~T["dual_flag"]
    rel = []
    for c in T.index:
        same = T[(T["top_set"] == T.loc[c, "top_set"]) & clean & (T.index != c)]
        ref = same["median_genes"].median() if len(same) else T["median_genes"].median()
        rel.append(float(T.loc[c, "median_genes"] / ref))
    T["rel_complexity"] = np.round(rel, 2)
    dbl_hi = T["scrublet_median"] > T["scrublet_median"].quantile(0.9)
    T["diagnosis"] = ""
    T.loc[T["dual_flag"] & (T["rel_complexity"] < 0.7) & ~dbl_hi, "diagnosis"] = "ambient (low complexity)"
    T.loc[T["dual_flag"] & (T["rel_complexity"] >= 0.9) & dbl_hi, "diagnosis"] = "doublet (normal complexity, high score)"
    T.loc[T["dual_flag"] & (T["diagnosis"] == ""), "diagnosis"] = "dual-lineage, mechanism unclear"

    # --- states (second axis; relative to the lineage's own median)
    T["proposed_state"], T["state_evidence"] = "", ""
    if states:
        sdet = {s_: cluster_detection(C, cl.values, [col[x] for x in gs]).median(axis=1) for s_, gs in states.items()}
        sdet = pd.DataFrame(sdet).reindex(T.index)
        for c in T.index:
            peers = T.index[T["top_set"] == T.loc[c, "top_set"]]
            base = sdet.loc[peers].median()
            fold = sdet.loc[c] / base.replace(0, np.nan)
            ok = [s_ for s_ in sdet.columns if sdet.loc[c, s_] >= STATE_FLOOR and fold[s_] >= STATE_FOLD.get(s_, 2.5)]
            if ok:
                best = max(ok, key=lambda s_: fold[s_])
                T.loc[c, "proposed_state"] = best
                T.loc[c, "state_evidence"] = f"{best}: median det {sdet.loc[c, best]:.0%}, {fold[best]:.1f}x lineage"
    mt_base = T.groupby("top_set")["median_pct_mt"].transform("median")
    mt_hi = (T["median_pct_mt"] >= 2 * mt_base) & (T["median_pct_mt"] >= 5) & (T["proposed_state"] == "")
    T.loc[mt_hi, "proposed_state"] = "Mitochondrial-high"
    T.loc[mt_hi, "state_evidence"] = "median %mt >= 2x lineage median"

    # --- proposals
    recs, confs, reasons = [], [], []
    for c, r in T.iterrows():
        why, rec = [], "keep"
        ident_txt = f"{r['top_set']} ({r['set_detection']}); margin {r['margin']}"
        if r["call"] == "unresolved":
            rec = "review"
            why.append(f"unresolved: best set {r['top_set']} best gene {r['best_gene']} {r['best_gene_det']:.0%}"
                       + (f"; veto {r['veto']}" if r["veto"] else ""))
        if r["dual_flag"]:
            rec = "remove"
            why.append(f"{r['diagnosis']}: dual {r['dual_frac']:.0%}, ratio {r['dual_ratio']:.2f} vs "
                       f"{r['second_lineage']} (dataset median dual {dual_med:.0%})")
        elif r["other_markers_top10"] >= 3 and r["other_markers_top10"] > r["own_markers_top10"]:
            rec = "remove"
            why.append(f"contamination: top-10 DE carries other-lineage markers {r['other_marker_genes_top30']}")
        elif r["frac_low_info_top10"] >= 0.5:
            rec = "remove"
            why.append(f"low information: {r['frac_low_info_top10']:.0%} of top-10 DE are RP/MT genes")
        elif r["below_size_drop"] and r["own_markers_top30"] < 2 and r["proposed_state"] == "":
            rec = "remove"
            why.append(f"below size drop (rank {r['size_rank']}, ratio {r['ratio_to_prev']:.2f}) and no own markers "
                       "in top-30 DE")
        elif r["below_size_drop"]:
            why.append(f"small (below size drop, rank {r['size_rank']}) but top-30 DE carry "
                       f"{r['own_markers_top30']} {r['top_set']} markers - check for a rare type")
        if r["top_sample_share"] >= 0.9:
            why.append(f"{r['top_sample_share']:.0%} {r['top_sample']} (sample == condition: may be condition biology)")
        if r["proposed_state"]:
            why.append(f"state {r['state_evidence']}")
        if rec == "keep" and r["call_path"] == "detection" and len(why) == 0:
            cf = "high"
        elif rec in ("keep", "remove") and r["call"] != "unresolved":
            cf = "medium" if len(why) <= 1 else "low"
        else:
            cf = "low"
        recs.append(rec)
        confs.append(cf)
        reasons.append(f"{ident_txt}. " + "; ".join(why))
    T["proposed_L1"] = np.where(T["call"] == "unresolved", "unresolved", T["call"])
    T["recommendation"], T["confidence"], T["reason"] = recs, confs, reasons
    for k in ["user_label", "user_state", "user_action", "user_notes"]:
        T[k] = ""
    lead = ["n_cells", "size_rank", "ratio_to_prev", "below_size_drop", "proposed_L1", "proposed_state",
            "recommendation", "confidence", "reason", "user_label", "user_state", "user_action", "user_notes"]
    T = T[lead + [c for c in T.columns if c not in lead]]
    T.index.name = "cluster"
    T.to_csv(p["review"] / "L1_review_table.csv")
    try:
        T.to_excel(p["review"] / "L1_review_table.xlsx")
    except Exception as e:  # openpyxl missing
        say(f"xlsx not written ({e})")
    # README
    lat = A.uns.get("L1_latent", {})
    log = json.loads((p["logs"] / "pipeline_log.json").read_text())
    emb = [x for x in log if x["stage"] == "embed" and x.get("round", 1) == ROUND][-1]
    dbl = [x for x in log if x["stage"] == "doublets"][-1]
    rc = T["recommendation"].value_counts().to_dict()
    lines = [
        f"# Mouse 2405072 - L1 lineage review (round {ROUND})", "",
        *([f"Round {ROUND}: re-embedded after removing clusters {emb.get('removed_clusters')} of round "
           f"{emb.get('from_round')} ({emb.get('n_removed')} cells: {emb.get('reason')}).", ""] if ROUND > 1 else []),
        f"Cells after QC + consensus-doublet removal: {A.n_obs:,}; samples: "
        + ", ".join(f"{k} {v:,}" for k, v in A.obs['sample_id'].value_counts().items()), "",
        f"- Ambient correction: skipped (filtered matrices only).",
        f"- Doublets: methods {dbl['methods']}, rule {dbl['consensus_rule']}.",
        f"- Latent: **{lat.get('choice')}** ({lat.get('why')}). Gates: `{json.dumps(emb['gates'], default=float)}`",
        f"- One sample per condition, so sample_id (Harmony key) == condition.",
        f"- L1 resolution: Leiden {L1_RES} ({len(T)} clusters); ladder {LEIDEN_LADDER} kept as obs columns.",
        f"- Size rule: clusters sorted by size, first n/n_prev < {SIZE_DROP_RATIO} at rank {drop_at + 1}; "
        "below-drop clusters proposed for removal unless their top-30 DE carry their lineage's markers.",
        "- DE: full-data Wilcoxon cluster-vs-rest on lognorm, ranked by scores, padj<0.05, logFC>0; "
        "RPS/RPL and mt- genes removed from DE tables and plots.",
        f"- Marker calls: AUCell (top {AUCELL_TOP_FRAC:.0%} genes) per cell, cluster mean z across clusters; "
        f"support floor {SUPPORT_FLOOR:.0%}, resolved by detection >= {DETECT_STRONG:.0%} or margin >= {MARGIN_MIN}.",
        f"- Recommendations: {rc}.", "",
        "Clusters needing judgement:", "",
        "```", T[(T["confidence"] != "high")][["n_cells", "proposed_L1", "proposed_state", "recommendation",
                                                "confidence", "reason"]].to_string(), "```",
    ]
    (p["review"] / "README_L1.md").write_text("\n".join(lines))
    A.obs[[key]].to_parquet(p["work"] / ("04_clusters.parquet" if ROUND == 1 else f"04_clusters_round{ROUND}.parquet"))
    say(f"review table: {len(T)} clusters, recommendations {rc}")
    print(T[["n_cells", "proposed_L1", "proposed_state", "recommendation", "confidence"]].to_string())
    log_event("evidence", round=ROUND, resolution=L1_RES, n_clusters=len(T), recommendations=rc, size_drop_rank=drop_at + 1,
              dual_median=dual_med)


# ----------------------------------------------------------------------------- stage: figures
UNIT_COLORS = {"epithelial": "#2a9d8f", "endothelial": "#e76f51", "stromal": "#8d6e63", "myeloid": "#457b9d",
               "tnk": "#e9c46a", "bplasma": "#9c6ade", "unresolved": "#d62828", "removed": "#cfcfcf"}


def dotplot(X, groups, rows, row_labels, row_colors, blocks, fname, title, outline=None, share=None):
    genes = [g for _, gs in blocks for g in gs]
    gidx = [X["col"][g] for g in genes]
    M = X["X"][:, gidx]
    mean, frac = np.zeros((len(rows), len(genes))), np.zeros((len(rows), len(genes)))
    gv = np.asarray(groups)
    for i, r in enumerate(rows):
        sub = M[gv == r]
        mean[i] = np.asarray(sub.mean(0)).ravel()
        frac[i] = np.asarray((sub > 0).mean(0)).ravel()
    sd = mean.std(0)
    zz = (mean - mean.mean(0)) / np.where(sd == 0, 1, sd)
    w = max(6, 0.22 * len(genes) + 3.5)
    h = max(3, 0.28 * len(rows) + 2.2)
    fig, ax = plt.subplots(figsize=(w, h))
    yy, xx = np.meshgrid(np.arange(len(rows)), np.arange(len(genes)), indexing="ij")
    ec = np.full(mean.shape, "none", dtype=object)
    if outline:
        for i, r in enumerate(rows):
            for j, g in enumerate(genes):
                if g in outline.get(r, set()):
                    ec[i, j] = "black"
    sca = ax.scatter(xx.ravel(), yy.ravel(), s=(frac.ravel() * 110) + 0.5, c=zz.ravel(), cmap="RdBu_r",
                     vmin=-2.5, vmax=2.5, edgecolors=ec.ravel(), linewidths=0.8)
    xl = [f"{g} ({share[g]})" if share else g for g in genes]
    ax.set_xticks(range(len(genes)), xl, rotation=90, fontsize=7)
    ax.set_yticks(range(len(rows)), row_labels, fontsize=7)
    for t, c in zip(ax.get_yticklabels(), row_colors):
        t.set_color(c)
    x0 = 0
    for nm, gs in blocks:
        ax.text(x0 + (len(gs) - 1) / 2, -1.0, nm, ha="center", va="bottom", fontsize=7, rotation=0 if len(gs) > 2 else 45)
        x0 += len(gs)
        ax.axvline(x0 - 0.5, color="grey", lw=0.5)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -1.6)
    cb = fig.colorbar(sca, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("z of mean (rows shown)", fontsize=7)
    for f_ in (0.2, 0.5, 1.0):
        ax.scatter([], [], s=f_ * 110, c="grey", label=f"{f_:.0%}")
    ax.legend(title="% detected", fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.06, 1))
    ax.set_title(title, fontsize=9, pad=24)
    fig.tight_layout()
    fig.savefig(fname, dpi=160)
    plt.close(fig)


def umap_plot(U, labels, fname, title, colors=None, text=True, order=None, legend_counts=True):
    lab = pd.Series(np.asarray(labels).astype(str))
    cats = order or sorted(lab.unique(), key=lambda x: (len(x), x))
    pal = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(8.5, 7))
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(lab))
    cmap = {c: (colors.get(c) if colors and c in colors else pal(i % 20)) for i, c in enumerate(cats)}
    ax.scatter(U[perm, 0], U[perm, 1], s=max(0.3, 30000 / len(lab) * 0.1), c=[cmap[x] for x in lab.values[perm]],
               linewidths=0, rasterized=True)
    for c in cats:
        sel = lab.values == c
        if text:
            ax.text(*np.median(U[sel], 0), c, fontsize=8, ha="center", va="center", weight="bold",
                    bbox=dict(facecolor="white", alpha=0.6, lw=0, pad=0.5))
        ax.scatter([], [], c=[cmap[c]], s=30, label=f"{c} ({sel.sum():,})" if legend_counts else c)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, markerscale=1)
    ax.set_xticks([]), ax.set_yticks([])
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(fname, dpi=160)
    plt.close(fig)


def stage_figures(args):
    p = paths()
    m = markers_table()
    A = sc.read_h5ad(emb_file(ROUND))
    key = f"leiden_{L1_RES}"
    cl = A.obs[key].astype(str).values
    T = pd.read_csv(p["review"] / "L1_review_table.csv", dtype={"cluster": str}).set_index("cluster")
    de = pd.read_csv(p["review"] / "L1_DE_top30.csv", dtype={"group": str})
    ident, states, _, _ = marker_sets(m, A.var_names)
    U = A.obsm["X_umap"]
    umap_plot(U, cl, p["figs"] / "umap_L1_clusters.png", f"L1 Leiden {L1_RES} ({len(T)} clusters)")
    umap_plot(U, A.obs["sample_id"].values, p["figs"] / "umap_L1_samples.png", "samples", text=False)
    lab = pd.Series(cl).map(lambda c: "removed" if T.loc[c, "recommendation"] == "remove" else T.loc[c, "proposed_L1"])
    umap_plot(U, lab.values, p["figs"] / "umap_L1_proposed_labels.png",
              "proposed L1 lineage (removal recommendations grey)", colors=UNIT_COLORS, text=False,
              order=[u for u in UNITS + ["unresolved", "removed"] if u in set(lab)])
    # rows ordered by proposed label then size
    unit_rank = {u: i for i, u in enumerate(UNITS + ["unresolved"])}
    rows = sorted(T.index, key=lambda c: (unit_rank.get(T.loc[c, "proposed_L1"], 9), -T.loc[c, "n_cells"]))
    rl = [f"{c} | {T.loc[c, 'n_cells']:,} | {T.loc[c, 'proposed_L1']}{' [' + T.loc[c, 'proposed_state'] + ']' if isinstance(T.loc[c, 'proposed_state'], str) and T.loc[c, 'proposed_state'] else ''} | m={T.loc[c, 'margin']:.2f} | {T.loc[c, 'recommendation']}" for c in rows]
    rcol = ["#d62828" if T.loc[c, "proposed_L1"] == "unresolved" else
            ("#888888" if T.loc[c, "recommendation"] == "remove" else "black") for c in rows]
    Xd = {"X": A.layers["lognorm"].tocsc(), "col": {g: i for i, g in enumerate(A.var_names)}}
    dotplot(Xd, cl, rows, rl, rcol, [(u, ident[u]) for u in UNITS], p["figs"] / "L1_marker_dotplot.png",
            "L1 panel markers by unit (dot = % detected; colour = z of mean across rows shown)")
    top = de[de["rank"] <= 10]
    share = top.groupby("gene")["group"].nunique().to_dict()
    blocks, seen = [], set()
    for c in rows:
        gs = [g for g in top.loc[top["group"] == c, "gene"] if g not in seen]
        seen.update(gs)
        if gs:
            blocks.append((c, gs))
    outline = {c: set(top.loc[top["group"] == c, "gene"]) for c in rows}
    dotplot(Xd, cl, rows, rl, rcol, blocks, p["figs"] / "L1_top_DEG_dotplot.png",
            "top-10 Wilcoxon DEGs per cluster (score-ranked, padj<0.05, RP/MT removed; own-row outlined; "
            "(k) = clusters sharing the gene in their top 10)", outline=outline, share=share)
    if states:
        dotplot(Xd, cl, rows, rl, rcol, list(states.items()), p["figs"] / "L1_state_genes_dotplot.png",
                "state genes (Proliferating / IFN / stress / metallothionein)")
    log_event("figures", files=sorted(x.name for x in p["figs"].glob("*.png")))
    say(f"figures in {p['figs']}")


# ----------------------------------------------------------------------------- stage: apply
UNIT_DISPLAY = {"epithelial": "Epithelial", "endothelial": "Endothelial", "stromal": "Stromal",
                "myeloid": "Myeloid", "tnk": "T/NK", "bplasma": "B/plasma"}
ROUTE_MARGIN = 1.0  # per-cell routing of state clusters: min z-margin of the best unit over the runner-up


def load_decisions(args, p):
    """Agent review table + user decisions -> one row per cluster with final_action/label/call/state/route."""
    rv = p["review"] / "L1_review_table.csv"
    T = pd.read_csv(rv, dtype=str, keep_default_na=False).set_index("cluster")
    f = Path(args.decisions) if args.decisions else \
        SCRIPT_DIR / ("user_decisions_L1.csv" if ROUND == 1 else f"user_decisions_L1_round{ROUND}.csv")
    if f.exists():
        U = pd.read_csv(f, dtype=str, keep_default_na=False).set_index("cluster")
        miss, extra = sorted(set(T.index) - set(U.index), key=int), sorted(set(U.index) - set(T.index), key=int)
        if miss or extra:
            raise SystemExit(f"{f.name} does not match the review table: missing clusters {miss}, unknown {extra}")
        for c in ["user_label", "user_call", "user_state", "user_action", "route", "user_notes"]:
            T[c] = U[c].reindex(T.index).fillna("") if c in U else ""
        src = str(f)
    else:  # fall back to the user_* columns filled in the review table itself
        f = Path(args.review) if args.review else p["review"] / "L1_review_table.xlsx"
        rd = pd.read_excel if f.suffix == ".xlsx" else pd.read_csv
        R = rd(f, dtype=str, keep_default_na=False).set_index("cluster")
        for c in ["user_label", "user_state", "user_action", "user_notes"]:
            T[c] = R[c].reindex(T.index).fillna("")
        T["user_call"], T["route"] = "", ""
        src = str(f)
    T["final_action"] = [a.strip().lower() or r for a, r in zip(T["user_action"], T["recommendation"])]
    T["route"] = [r.strip() or "cluster" for r in T["route"]]
    T["final_label"] = [ul.strip() or (pl if rt == "cluster" else "")
                        for ul, pl, rt in zip(T["user_label"], T["proposed_L1"], T["route"])]
    T["final_state"] = [us.strip() or ps for us, ps in zip(T["user_state"], T["proposed_state"])]
    T["final_call"] = [uc.strip() or UNIT_DISPLAY.get(fl, "") for uc, fl in zip(T["user_call"], T["final_label"])]
    T["user_notes"] = [n or ("" if a else "no user entry: agent proposal taken")
                       for n, a in zip(T["user_notes"], T["user_action"])]
    bad = T[~T["final_action"].isin(["keep", "remove"])
            | ((T["final_action"] == "keep") & (T["route"] == "cluster") & ~T["final_label"].isin(UNITS))
            | ~T["route"].isin(["cluster", "per_cell"])]
    if len(bad):
        raise SystemExit("clusters without a valid decision (action keep/remove; keep needs a unit label or "
                         f"route=per_cell): {list(bad.index)}")
    return T, src


def stage_apply(args):
    p = paths()
    T, src = load_decisions(args, p)
    today = str(date.today())
    T.assign(decided_on=today).to_csv(p["review"] / f"user_decisions_L1_{L1_RES}.csv")
    say(f"applying L1 round {ROUND} decisions from {src}")
    A = sc.read_h5ad(emb_file(ROUND))
    key = f"leiden_{L1_RES}"
    cl = A.obs[key].astype(str)
    qc = sc.read_h5ad(p["work"] / "01_qc.h5ad", backed="r").obs
    D = pd.read_parquet(p["work"] / "02_doublets.parquet").set_index("cell_id")
    L = pd.DataFrame({"sample_id": qc["sample_id"].values, "condition": qc["condition"].values}, index=qc.index)
    L.index.name = "cell_id"
    L["L1_cluster"] = cl.reindex(L.index).fillna("").astype(str)
    for c in ["L1_lineage", "L1_call", "L1_state", "L1_route", "L1_removal_reason", "L1_label_before_removal"]:
        L[c] = ""
    L["L1_status"] = "kept"
    L["L1_round"] = ""
    L.loc[cl.index, "L1_round"] = str(ROUND)
    L["L1_route_margin"] = np.nan
    noqc = ~qc["qc_pass"].values
    L.loc[noqc, "L1_status"] = "removed_qc"
    L.loc[noqc, "L1_removal_reason"] = "QC: " + qc.loc[noqc, ["qc_fail_counts", "qc_fail_genes", "qc_fail_mt",
                                                             "qc_fail_mingenes"]].apply(
        lambda r: "/".join(k.replace("qc_fail_", "") for k, v in r.items() if v), axis=1)
    dbl = D.index[D["doublet_consensus"].astype(bool)]
    L.loc[dbl, "L1_status"] = "removed_doublet"
    L.loc[dbl, "L1_removal_reason"] = "doublet consensus " + D.loc[dbl, "n_doublet_calls"].astype(str) + "/" + \
        D.loc[dbl, "n_doublet_methods"].astype(str)
    for j in range(2, ROUND + 1):  # cells dropped by earlier L1 recluster rounds
        rmc = pd.read_csv(p["work"] / f"L1_round{j}_removed_cells.csv", dtype=str)
        L.loc[rmc["cell_id"], "L1_status"] = "removed"
        L.loc[rmc["cell_id"], "L1_round"] = rmc["from_round"].values
        L.loc[rmc["cell_id"], "L1_cluster"] = rmc["cluster"].values
        L.loc[rmc["cell_id"], "L1_removal_reason"] = rmc["reason"].values
    unacc = L.index[(L["L1_status"] == "kept") & ~L.index.isin(cl.index)]
    assert len(unacc) == 0, f"{len(unacc)} QC/doublet-passing cells are in no L1 round and no removed-cells file"
    routing = []
    for c, r in T.iterrows():
        cells = cl.index[cl.values == c]
        tag = f"L1 round{ROUND} leiden_{L1_RES}_cluster_{c}"
        L.loc[cells, "L1_route"] = r["route"]
        L.loc[cells, "L1_state"] = r["final_state"]
        if r["final_action"] == "remove":
            L.loc[cells, "L1_status"] = "removed"
            L.loc[cells, "L1_removal_reason"] = f"{tag}: {r['user_notes'] or r['reason']}"
            L.loc[cells, "L1_label_before_removal"] = r["final_call"] or r["proposed_L1"]
            L.loc[cells, "L1_state"] = ""
        elif r["route"] == "cluster":
            L.loc[cells, "L1_lineage"] = r["final_label"]
            L.loc[cells, "L1_call"] = r["final_call"]
        else:  # state cluster: lineage per cell from the L1 panel AUCell z (embed stage)
            u = A.obs.loc[cells, "cell_unit"].astype(str)
            mg = A.obs.loc[cells, "cell_unit_margin"].astype(float)
            ok = mg >= ROUTE_MARGIN
            L.loc[cells, "L1_route_margin"] = mg.values
            L.loc[cells[ok.values], "L1_lineage"] = u[ok].values
            L.loc[cells[ok.values], "L1_call"] = u[ok].map(UNIT_DISPLAY).values
            low = cells[~ok.values]
            L.loc[low, "L1_status"] = "removed"
            L.loc[low, "L1_removal_reason"] = [f"{tag} ({r['final_state']}): per-cell lineage unresolved, "
                                               f"margin {m:.2f} < {ROUTE_MARGIN}" for m in mg[~ok]]
            L.loc[low, "L1_label_before_removal"] = u[~ok].map(UNIT_DISPLAY).values
            L.loc[low, "L1_state"] = ""
            vc = u[ok].value_counts()
            routing.append({"cluster": c, "state": r["final_state"], "n_cells": len(cells),
                            "n_unresolved": int((~ok).sum()), **{k: int(vc.get(k, 0)) for k in UNITS}})
    L["decided_on"] = today
    assert L.index.is_unique
    kept = L["L1_status"] == "kept"
    assert L.loc[kept, "L1_lineage"].isin(UNITS).all(), "kept cell without a unit"
    assert (L.loc[kept, "L1_call"] != "").all(), "kept cell without an L1 call"
    assert (L.loc[~kept, "L1_removal_reason"] != "").all(), "removed cell without a reason"
    assert int(kept.sum()) + int((~kept).sum()) == qc.shape[0]
    L.reset_index().to_parquet(p["annot"] / "L1_labels.parquet", index=False)
    if routing:
        R = pd.DataFrame(routing)
        R.to_csv(p["annot"] / "L1_state_cluster_routing.csv", index=False)
        print("state clusters routed per cell:\n" + R.to_string(index=False))
    for u in UNITS:
        L[kept & (L["L1_lineage"] == u)].reset_index()[["cell_id", "sample_id", "condition", "L1_cluster",
                                                        "L1_call", "L1_state"]].to_parquet(
            p["units"] / f"{u}_cells.parquet", index=False)
    emb_kept = A.obs_names[L.loc[A.obs_names, "L1_status"].values == "kept"]
    U_ = A[emb_kept].obsm["X_umap"]
    umap_plot(U_, L.loc[emb_kept, "L1_lineage"].values, p["figs"] / "umap_L1_final.png", "L1 lineages (agreed)",
              colors=UNIT_COLORS, text=False, order=[u for u in UNITS if u in set(L.loc[emb_kept, "L1_lineage"])])
    umap_plot(U_, L.loc[emb_kept, "L1_call"].values, p["figs"] / "umap_L1_final_calls.png", "L1 calls (agreed)")
    st = L.loc[emb_kept, "L1_state"].replace("", "none").values
    umap_plot(U_, st, p["figs"] / "umap_L1_final_states.png", "L1 states", text=False,
              colors={"none": "#e0e0e0", "Proliferating": "#d62828", "IFN response": "#1d3557"})
    ct = pd.crosstab([L["L1_lineage"].replace("", "(removed)"), L["L1_call"].replace("", "-")], L["sample_id"])
    ct.to_csv(p["annot"] / "L1_counts_by_sample.csv")
    print(L["L1_status"].value_counts().to_string())
    print(ct.to_string())
    log_event("apply", round=ROUND, decisions=src, route_margin=ROUTE_MARGIN, status=L["L1_status"].value_counts().to_dict(),
              lineage_counts=L.loc[kept, "L1_lineage"].value_counts().to_dict(),
              call_counts=L.loc[kept, "L1_call"].value_counts().to_dict(), routing=routing)


# ----------------------------------------------------------------------------- stage: split
def stage_split(args):
    """One h5ad per L1 lineage for L2: kept cells, raw counts (X and layers['counts']), all measured genes,
    L1 metadata + QC + doublet columns, the L1 UMAP / PCA / Harmony latents (L2 re-derives HVGs and latents)."""
    p = paths()
    L = pd.read_parquet(p["annot"] / "L1_labels.parquet").set_index("cell_id")
    kept = L.index[L["L1_status"] == "kept"]
    rnd = sorted(set(L.loc[kept, "L1_round"].astype(str))) if "L1_round" in L else ["1"]
    assert len(rnd) == 1, f"kept cells come from several L1 rounds {rnd}"
    rnd = int(rnd[0])
    A = sc.read_h5ad(emb_file(rnd))
    say(f"splitting L1 round {rnd} ({len(kept):,} kept cells)")
    A = A[A.obs_names.isin(kept)].copy()
    counts = A.layers["counts"].astype(np.float32)
    A = ad.AnnData(X=counts, obs=A.obs, var=A.var, obsm=dict(A.obsm), layers={"counts": counts.copy()})
    keep_obs = ["barcode", "sample_id", "condition", "batch", "n_genes_by_counts", "total_counts", "pct_counts_mt",
                "pct_counts_ribo", "pct_counts_hb", "scrublet_score", "scdbl_score", "df_pann", "n_doublet_calls",
                "n_doublet_methods", "doublet_candidate"] + [f"leiden_{r}" for r in LEIDEN_LADDER] + \
               [f"auc_{u}" for u in UNITS] + ["cell_unit", "cell_unit_margin"]
    obs = A.obs[[c for c in keep_obs if c in A.obs]].copy()
    for c in ["L1_round", "L1_cluster", "L1_lineage", "L1_call", "L1_state", "L1_route", "L1_route_margin",
              "decided_on"]:
        if c in L:
            obs[c] = L.loc[A.obs_names, c].values
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].astype(str).astype("category")
    A.obs = obs
    for k in [k for k in A.obsm if k not in ("X_umap", "X_pca", "X_pca_harmony")]:
        del A.obsm[k]
    A.obsm["X_umap_L1"] = A.obsm.pop("X_umap")
    A.uns = {"L1_info": {"source": "l1_mouse.py", "L1_round": rnd, "L1_resolution": L1_RES, "latent": "see logs/pipeline_log.json",
                         "X": "raw counts (float32); layers['counts'] identical",
                         "decided_on": str(L["decided_on"].iloc[0])}}
    A.obsp = {}
    out = OUT_DIR / "L1" / "lineages"
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    for u in UNITS:
        a = A[A.obs["L1_lineage"] == u].copy()
        if a.n_obs == 0:
            say(f"{u}: no cells, skipped")
            continue
        for c in a.obs.columns:
            if hasattr(a.obs[c], "cat"):
                a.obs[c] = a.obs[c].cat.remove_unused_categories()
        assert np.allclose(a.X.data, np.round(a.X.data)), "X is not raw counts"
        f = out / f"mouse_2405072_{u}_L1.h5ad"
        a.write_h5ad(f, compression="gzip")
        written[u] = {"file": str(f), "n_cells": int(a.n_obs), "n_genes": int(a.n_vars),
                      "calls": a.obs["L1_call"].value_counts().to_dict(),
                      "states": a.obs["L1_state"].astype(str).replace("", "none").value_counts().to_dict(),
                      "by_sample": a.obs["sample_id"].value_counts().to_dict()}
        say(f"{u}: {a.n_obs:,} cells -> {f.name}  calls {written[u]['calls']}")
    total = sum(v["n_cells"] for v in written.values())
    assert total == len(kept), f"split cells {total} != kept {len(kept)}"
    log_event("split", round=rnd, files=written)


# ----------------------------------------------------------------------------- main
def main():
    global DATA_DIR, OUT_DIR, SAMPLE_SHEET, ROUND
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True, choices=["markers", "qc", "doublets", "embed", "evidence", "figures",
                                                        "apply", "split", "recluster"])
    ap.add_argument("--data-dir", default=str(DATA_DIR))
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--sample-sheet", default=str(SAMPLE_SHEET))
    ap.add_argument("--marker-md", help="mouse panel in the lab markdown format (default: marker_panel_mouse_v2.md "
                    "next to this script, else <repo>/references/; used when no other panel is given)")
    ap.add_argument("--human-panel", help="lab human L1 marker panel CSV (gene + compartment columns)")
    ap.add_argument("--mouse-panel", help="mouse marker CSV with columns gene,set,source[,required]")
    ap.add_argument("--orthologs", help="optional human->mouse ortholog CSV (columns containing 'human' and 'mouse')")
    ap.add_argument("--review", help="filled review table for --stage apply (used when no decisions CSV exists)")
    ap.add_argument("--decisions", help="user decisions CSV for --stage apply (default: user_decisions_L1.csv next "
                    "to this script); columns cluster,user_label,user_call,user_state,user_action,route,user_notes")
    ap.add_argument("--round", type=int, help="L1 round for evidence / figures / apply (default: latest)")
    ap.add_argument("--remove", help="recluster: comma-separated L1 cluster ids to drop, e.g. 3,7")
    ap.add_argument("--from-round", type=int, help="recluster: round the ids refer to (default: latest)")
    ap.add_argument("--from-res", type=float, help=f"recluster: Leiden resolution of the ids (default {L1_RES})")
    ap.add_argument("--reason", default="removed by user", help="recluster: reason recorded for removed cells")
    ap.add_argument("--force", action="store_true", help="continue past a failed gate (logged)")
    a = ap.parse_args()
    DATA_DIR, OUT_DIR, SAMPLE_SHEET = Path(a.data_dir), Path(a.out_dir), Path(a.sample_sheet)
    sc.settings.n_jobs = 8
    if a.stage in ("evidence", "figures", "apply"):
        ROUND = a.round or latest_l1_round()
    say(f"stage {a.stage}: data {DATA_DIR} -> out {OUT_DIR}" + (f" (L1 round {ROUND})" if ROUND > 1 else ""))
    {"markers": stage_markers, "qc": stage_qc, "doublets": stage_doublets, "embed": stage_embed,
     "evidence": stage_evidence, "figures": stage_figures, "apply": stage_apply,
     "split": stage_split, "recluster": stage_recluster}[a.stage](a)
    say("done")


if __name__ == "__main__":
    main()
