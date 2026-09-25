#!/usr/bin/env python3
"""
substitution_ranker.py
======================

Copyright (c) 2026 Christien Dykstra
SPDX-License-Identifier: MIT

Activation-category chart for amino-acid substitutions in GPCRs, as used in a delta-opioid
receptor (OPRD1) directed-evolution screen (fig5f).

Each substitution is followed across ALL genotypes it appears in (substitutions
sit inside larger, complex genotypes -- that is the point). Every carrier
genotype is placed in exactly ONE activation category, by the LOWEST agonist
concentration at which it clears the fold-activation gate:

    4 uM   -> activates under LIMITING agonist        (strong signaller)
    100 uM -> activates only under SATURATING agonist (needs a high dose)
    0 uM   -> never activates at any tested dose      (dead in that genotype)

Categories are mutually exclusive, so the three columns sum to the pool (total
carriers). Colour and marker size encode how carriers distribute across them.

Reading a row:
    weight at 4 uM          -> carriers work even at limiting agonist  -> driver
    weight at 0 uM          -> carriers never fire                     -> passenger
    spread across all three -> presence does not decide the outcome    -> backbone

Only substitutions carried by >= --min-occurrence genotypes are shown.


INPUT FORMAT
------------
A single-sheet .xlsx with one row per well/genotype. Required columns:

    AA_substitutions_nonsynonymous   comma- or semicolon-separated substitutions
                                     in single-letter form (e.g. "A51T, N131S").
                                     Entries not matching <AA><pos><AA> are skipped.
    Colony_barcode                   well identifier; rows whose barcode contains
                                     --wt-pattern are treated as wild-type controls

Plus, for each concentration in --concs, the columns required by the chosen gate:

    gate=threshold (default)   {conc}_meets_fold_threshold    boolean
    gate=mu_max                {conc}_mu_max_per_hr           float
    gate=cutoff                {conc}_{--metric}              float

`{conc}_mu_max_per_hr` is also used for dot colour in the default --style dot,
and is optional otherwise (missing values render grey).

Concentration labels come from --concs and must match the column prefixes
exactly, e.g. "0uM 4uM 100uM" for columns like `4uM_meets_fold_threshold`.


OUTPUT
------
A timestamped folder created next to the input file (or at --outdir) holding:
the figure in each requested format, a CSV of the count matrix with per-category
mean mu_max, and a text dump of every parameter used for that run.


USAGE
-----
    python substitution_ranker.py data.xlsx
    python substitution_ranker.py data.xlsx --value fraction
    python substitution_ranker.py data.xlsx --min-occurrence 5 --sort enrichment
    python substitution_ranker.py data.xlsx --gate cutoff --metric AUC_fc --cutoff 35
    python substitution_ranker.py data.xlsx --style gradient

Requires: pandas numpy matplotlib openpyxl
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, PowerNorm

# ---------------------------------------------------------------- CONFIG -----
# Defaults for every option. Anything marked * can be overridden on the command
# line; run with --help for the full list.
CONFIG = dict(
    input=None,                # * path to the .xlsx (required)
    sheet=0,                   # * sheet name or index
    concs=("0uM", "4uM", "100uM"),   # * low -> high; first is the "never" bin

    min_occurrence=6,          # * show substitutions carried by >= N genotypes

    # column ordering -----------------------------------------------------
    order="sensitivity",       # * sensitivity (inactive -> 100uM -> 4uM) | concentration
    include_wt=True,           # * add a wild-type reference row
    wt_pattern="WT",           # * Colony_barcode substring identifying WT controls

    # how "activates"/"fires" is defined for the 4 / 100 uM gates ----------
    gate="threshold",          # * threshold (FC_max>=2 flag) | mu_max | cutoff
    metric="AUC_fc",           # * used only when gate == "cutoff"
    cutoff=35.0,               # * used only when gate == "cutoff"
    mu_threshold=0.08,         # * mu_max gate (per hour); used when gate == "mu_max"

    # cell value & colour -------------------------------------------------
    value="fraction",          # * count | fraction  (fraction = carriers / pool)
    cmap="magma",              # *
    norm="power",              # * linear | power
    gamma=0.5,                 # * PowerNorm lifts small counts into visible range
    vmax_pct=90.0,             # * clip colour scale so the biggest pool doesn't wash out

    style="dot",               # * dot (size=fraction, colour=mu_max) | blocks | gradient
    dot_size_max=1100.0,       # scatter area (pts^2) at fraction = 1.0
    dot_size_min=30.0,         # floor so small fractions stay visible

    # ordering (top -> bottom) --------------------------------------------
    sort="enrichment",         # * enrichment(4uM/pool) | active(4uM) | pool | position | label
    descending=True,           # *

    # annotations ---------------------------------------------------------
    annotate=True,             # * print the count (or fraction) inside each cell
    show_pool=True,            # * print total carrier pool on the right
    lane_frac=0.90,

    out_stem="substitution_ranker",   # *
    outdir=None,               # * default: alongside the input file
    dpi=300,                   # *
    formats=("png", "svg", "pdf"),  # *
)

RUN_STAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

SUB_RE = re.compile(r"^[A-Z]\d+[A-Z]$")
POS_RE = re.compile(r"(\d+)")


# ---------------------------------------------------------------- data -------
def split_subs(cell):
    if pd.isna(cell):
        return []
    return [p.strip() for p in str(cell).replace(";", ",").split(",") if p.strip()]


def activates(row, conc, cfg):
    """True if this well clears the activation gate at the given concentration."""
    if cfg["gate"] == "mu_max":
        v = row[f"{conc}_mu_max_per_hr"]
        return bool(pd.notna(v) and v >= cfg["mu_threshold"])
    if cfg["gate"] == "cutoff":
        v = row[f"{conc}_{cfg['metric']}"]
        return bool(pd.notna(v) and v >= cfg["cutoff"])
    v = row[f"{conc}_meets_fold_threshold"]
    return bool(v) if pd.notna(v) else False


def fire_rule(cfg):
    """Short human description of the current activation gate."""
    if cfg["gate"] == "mu_max":
        return f"fire = \u00b5\u2098\u2090\u2093 \u2265 {cfg['mu_threshold']:g} h\u207b\u00b9"
    if cfg["gate"] == "cutoff":
        return f"fire = {cfg['metric']} \u2265 {cfg['cutoff']:g}"
    return "fire = FC_max \u2265 2 (fold gate)"


def make_mu_norm(mumat, cfg):
    finite = mumat[np.isfinite(mumat)]
    vmax = np.percentile(finite, cfg["vmax_pct"]) if finite.size else 1.0
    return Normalize(vmin=0.0, vmax=max(vmax, 1e-9))   # vmin=0 -> mu_max~0 is black


def build(df, cfg):
    concs = list(cfg["concs"])
    lo, hi = concs[1], concs[2]

    if cfg["gate"] == "cutoff":
        need = [f"{c}_{cfg['metric']}" for c in (lo, hi)]
    elif cfg["gate"] == "mu_max":
        need = [f"{c}_mu_max_per_hr" for c in (lo, hi)]
    else:
        need = [f"{c}_meets_fold_threshold" for c in (lo, hi)]
    missing = [c for c in need if c not in df.columns]
    if missing:
        sys.exit(f"[error] missing columns for gate='{cfg['gate']}': {missing}")

    for required in ("AA_substitutions_nonsynonymous",):
        if required not in df.columns:
            sys.exit(f"[error] missing required column: {required}")

    df = df.copy()
    df["_subs"] = df["AA_substitutions_nonsynonymous"].apply(split_subs)

    def bin_of(r):
        if activates(r, lo, cfg):
            return lo
        if activates(r, hi, cfg):
            return hi
        return concs[0]

    def mu_of(r, cat):
        if cat == concs[0]:                 # never activates -> no signal -> black
            return 0.0
        col = f"{cat}_mu_max_per_hr"
        return r[col] if col in df.columns else np.nan

    recs = []
    for _, r in df.iterrows():
        cat = bin_of(r)
        mu = mu_of(r, cat)
        for s in r["_subs"]:
            if SUB_RE.match(s):
                recs.append((s, cat, mu))
    if not recs:
        sys.exit("[error] no substitutions parsed.")

    wt = None
    if cfg["include_wt"]:
        if "Colony_barcode" not in df.columns:
            print("[warn] no Colony_barcode column; skipping WT reference row")
        else:
            mask = df["Colony_barcode"].astype(str).str.contains(
                cfg["wt_pattern"], case=False, na=False)
            wt_wells = df[mask]
            if len(wt_wells):
                wt_cats = wt_wells.apply(bin_of, axis=1)
                wt_c = np.array([float((wt_cats == c).sum()) for c in concs])
                wt_mu = np.array([
                    float(wt_wells.loc[wt_cats == c]
                          .apply(lambda rr: mu_of(rr, c), axis=1).mean())
                    if (wt_cats == c).any() else np.nan
                    for c in concs])
                wt = ("WT", wt_c, int(wt_c.sum()), wt_mu)
            else:
                print(f"[warn] no rows matched wt_pattern={cfg['wt_pattern']!r}")

    e = pd.DataFrame(recs, columns=["sub", "cat", "mu"])
    ct = (e.pivot_table(index="sub", columns="cat", aggfunc="size", fill_value=0)
            .reindex(columns=concs, fill_value=0))
    mm = (e.pivot_table(index="sub", columns="cat", values="mu", aggfunc="mean")
            .reindex(index=ct.index, columns=concs))
    pool = ct.sum(axis=1)
    keep = pool >= cfg["min_occurrence"]
    ct, mm, pool = ct[keep], mm[keep], pool[keep]
    if ct.empty:
        sys.exit(f"[error] nothing reaches min_occurrence={cfg['min_occurrence']}.")

    cmat = ct.to_numpy(float)
    mumat = mm.to_numpy(float)
    pool_v = pool.to_numpy(float)
    vmat = cmat / pool_v[:, None] if cfg["value"] == "fraction" else cmat
    labels = list(ct.index)

    def position(l):
        m = POS_RE.search(l)
        return int(m.group(1)) if m else 10 ** 9

    key = {
        "enrichment": lambda i: cmat[i, 1] / pool_v[i],
        "active":     lambda i: cmat[i, 1],
        "pool":       lambda i: pool_v[i],
        "position":   lambda i: position(labels[i]),
        "label":      lambda i: labels[i],
    }[cfg["sort"]]
    idx = sorted(range(len(labels)), key=key, reverse=cfg["descending"])
    order = np.array(idx)
    return ([labels[i] for i in order], vmat[order], cmat[order],
            mumat[order], pool_v[order].astype(int), wt)


# ---------------------------------------------------------------- plot -------
def make_norm(vmat, cfg):
    finite = vmat[np.isfinite(vmat)]
    vmax = np.percentile(finite, cfg["vmax_pct"]) if finite.size else 1.0
    vmax = max(vmax, 1e-9)
    if cfg["norm"] == "power":
        return PowerNorm(gamma=cfg["gamma"], vmin=0.0, vmax=vmax)
    return Normalize(vmin=0.0, vmax=vmax)


def fmt(v, cfg):
    return f"{v:.2f}" if cfg["value"] == "fraction" else f"{int(round(v))}"


def plot(labels, vmat, cmat, mumat, pool, wt, cfg):
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "DejaVu Sans"

    has_wt = wt is not None
    if has_wt:
        wt_lab, wt_c, wt_pool, wt_mu = wt
        wt_v = wt_c / wt_pool if cfg["value"] == "fraction" else wt_c
        labels = [wt_lab] + list(labels)
        vmat = np.vstack([wt_v, vmat])
        cmat = np.vstack([wt_c, cmat])
        mumat = np.vstack([wt_mu, mumat])
        pool = np.concatenate([[wt_pool], pool])

    if cfg["order"] == "sensitivity":
        perm = [0, 2, 1]
        col_labels = ["inactive\n(no activation)", "low\n(needs 100 \u00b5M)",
                      "high\n(fires at 4 \u00b5M)"]
        xtitle = "Agonist sensitivity  \u2192"
    else:
        perm = [0, 1, 2]
        col_labels = ["0 \u00b5M\n(no activation)", "4 \u00b5M\n(limiting)",
                      "100 \u00b5M\n(saturating)"]
        xtitle = "Lowest activating agonist concentration"
    vmat, cmat, mumat = vmat[:, perm], cmat[:, perm], mumat[:, perm]

    n, ncol = vmat.shape
    norm = make_norm(vmat, cfg)
    munorm = None
    cmap = plt.get_cmap(cfg["cmap"]).copy()
    cmap.set_bad("0.85")

    fig_h = max(2.5, 0.40 * n + 1.7)
    fig, ax = plt.subplots(figsize=(6.6, fig_h))
    lane = cfg["lane_frac"]
    node_x = np.arange(ncol)
    xmin, xmax = -0.5, ncol - 0.5

    if cfg["style"] == "gradient":
        K = 400
        xpix = np.linspace(xmin, xmax, K)
        for i, row in enumerate(vmat):
            strip = np.interp(xpix, node_x, row)[None, :]
            ax.imshow(strip, aspect="auto", cmap=cmap, norm=norm,
                      extent=[xmin, xmax, i + lane / 2, i - lane / 2],
                      interpolation="bilinear", zorder=1)
    elif cfg["style"] == "dot":
        munorm = make_mu_norm(mumat, cfg)
        s_max, s_min = cfg["dot_size_max"], cfg["dot_size_min"]
        ok = dict(x=[], y=[], s=[], c=[])
        na = dict(x=[], y=[], s=[])
        for i in range(n):
            for j in range(ncol):
                cnt = cmat[i, j]
                if cnt <= 0:
                    continue
                s = s_min + (s_max - s_min) * (cnt / pool[i])
                if np.isfinite(mumat[i, j]):
                    ok["x"].append(node_x[j]); ok["y"].append(i)
                    ok["s"].append(s); ok["c"].append(mumat[i, j])
                else:
                    na["x"].append(node_x[j]); na["y"].append(i); na["s"].append(s)
        if na["x"]:
            ax.scatter(na["x"], na["y"], s=na["s"], facecolor="0.85",
                       edgecolor="0.5", linewidth=0.6, zorder=3)
        ax.scatter(ok["x"], ok["y"], s=ok["s"], c=ok["c"], cmap=cmap,
                   norm=munorm, edgecolor="0.2", linewidth=0.6, zorder=3)
    else:  # blocks
        for i, row in enumerate(vmat):
            for j, v in enumerate(row):
                col = cmap(norm(v))
                ax.add_patch(plt.Rectangle((node_x[j] - lane / 2, i - lane / 2),
                                           lane, lane, facecolor=col,
                                           edgecolor="white", lw=1.3, zorder=1))
                if cfg["annotate"]:
                    lum = 0.299 * col[0] + 0.587 * col[1] + 0.114 * col[2]
                    ax.text(node_x[j], i, fmt(v, cfg), ha="center", va="center",
                            fontsize=7.5, color="white" if lum < 0.55 else "black",
                            zorder=3)

    if cfg["show_pool"]:
        for i, p in enumerate(pool):
            ax.text(xmax + 0.14, i, f"n={int(p)}", ha="left", va="center",
                    fontsize=8, color="0.35")

    if has_wt:
        ax.axhline(0.5, color="0.6", lw=0.8, ls=(0, (4, 3)), zorder=2)

    ax.set_xlim(xmin - 0.10, xmax + (0.95 if cfg["show_pool"] else 0.1))
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(node_x)
    ax.set_xticklabels(col_labels)
    ax.set_yticks(range(n))
    yt = ax.set_yticklabels(labels)
    if has_wt:
        yt[0].set_color("0.40")
        yt[0].set_fontstyle("italic")
    ax.tick_params(length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    ax.set_xlabel(xtitle, fontsize=10, labelpad=8)
    ax.set_ylabel("Amino-acid substitution", fontsize=10)
    what = ("mean \u00b5\u2098\u2090\u2093 / carrier fraction" if cfg["style"] == "dot"
            else "carrier fraction" if cfg["value"] == "fraction" else "carrier genotypes")
    ax.set_title(f"OPRD1 substitution activation category  \u00b7  {what}",
                 fontsize=11, pad=18)
    ax.text(0.5, 1.008, fire_rule(cfg), transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8.5, color="0.45")

    if cfg["style"] == "dot":
        cbar_norm, cbar_label = munorm, "mean \u00b5\u2098\u2090\u2093 (h\u207b\u00b9)"
    else:
        cbar_norm = norm
        cbar_label = ("carrier fraction" if cfg["value"] == "fraction"
                      else "carrier genotypes (count)")
    sm = ScalarMappable(norm=cbar_norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label(cbar_label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_visible(False)

    if cfg["style"] == "dot":
        from matplotlib.lines import Line2D
        s_max, s_min = cfg["dot_size_max"], cfg["dot_size_min"]
        handles = []
        for f in (0.25, 0.5, 0.75, 1.0):
            s = s_min + (s_max - s_min) * f
            handles.append(Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor="0.55", markeredgecolor="0.2",
                                  markersize=2 * np.sqrt(s / np.pi),
                                  label=f"{int(f * 100)}%"))
        ax.legend(handles=handles, title="carriers\n(fraction of pool)",
                  loc="upper left", bbox_to_anchor=(1.16, 1.0), frameon=False,
                  labelspacing=1.6, borderpad=1.0, fontsize=8, title_fontsize=8)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- cli --------
def parse_args():
    p = argparse.ArgumentParser(
        description="Activation-category chart for amino-acid substitutions in a "
                    "directed-evolution screen.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Input: single-sheet .xlsx with AA_substitutions_nonsynonymous, "
               "Colony_barcode, and per-concentration gate columns. "
               "See the module docstring for the full format.")
    p.add_argument("input", nargs="?", default=CONFIG["input"],
                   help="Path to the input .xlsx file")
    p.add_argument("--sheet", default=CONFIG["sheet"],
                   help="Sheet name or index (default: first sheet)")
    p.add_argument("--concs", nargs=3, default=list(CONFIG["concs"]),
                   metavar=("NEVER", "LOW", "HIGH"),
                   help="Column prefixes for the three categories, low to high "
                        f"(default: {' '.join(CONFIG['concs'])})")
    p.add_argument("--min-occurrence", type=int, default=CONFIG["min_occurrence"])
    p.add_argument("--order", default=CONFIG["order"],
                   choices=["sensitivity", "concentration"])
    p.add_argument("--no-wt", dest="include_wt", action="store_false",
                   default=CONFIG["include_wt"])
    p.add_argument("--wt-pattern", default=CONFIG["wt_pattern"])
    p.add_argument("--gate", default=CONFIG["gate"],
                   choices=["threshold", "mu_max", "cutoff"])
    p.add_argument("--metric", default=CONFIG["metric"])
    p.add_argument("--cutoff", type=float, default=CONFIG["cutoff"])
    p.add_argument("--mu-threshold", type=float, default=CONFIG["mu_threshold"])
    p.add_argument("--value", default=CONFIG["value"], choices=["count", "fraction"])
    p.add_argument("--cmap", default=CONFIG["cmap"])
    p.add_argument("--norm", default=CONFIG["norm"], choices=["linear", "power"])
    p.add_argument("--gamma", type=float, default=CONFIG["gamma"])
    p.add_argument("--vmax-pct", type=float, default=CONFIG["vmax_pct"])
    p.add_argument("--style", default=CONFIG["style"],
                   choices=["dot", "blocks", "gradient"])
    p.add_argument("--sort", default=CONFIG["sort"],
                   choices=["enrichment", "active", "pool", "position", "label"])
    p.add_argument("--ascending", dest="descending", action="store_false",
                   default=CONFIG["descending"])
    p.add_argument("--no-annotate", dest="annotate", action="store_false",
                   default=CONFIG["annotate"])
    p.add_argument("--no-pool", dest="show_pool", action="store_false",
                   default=CONFIG["show_pool"])
    p.add_argument("--out-stem", default=CONFIG["out_stem"])
    p.add_argument("--outdir", default=CONFIG["outdir"])
    p.add_argument("--dpi", type=int, default=CONFIG["dpi"])
    p.add_argument("--formats", nargs="+", default=list(CONFIG["formats"]))
    return p.parse_args()


def main():
    cfg = dict(CONFIG)
    cfg.update(vars(parse_args()))
    cfg["concs"] = tuple(cfg["concs"])

    if not cfg["input"]:
        sys.exit("[error] no input file given. Usage: "
                 "python substitution_ranker.py <data.xlsx>  (see --help)")

    path = Path(cfg["input"]).expanduser()
    if not path.exists():
        sys.exit(f"[error] input not found: {path.resolve()}")

    try:
        sheet = int(cfg["sheet"])
    except (TypeError, ValueError):
        sheet = cfg["sheet"]
    df = pd.read_excel(path, sheet_name=sheet)

    labels, vmat, cmat, mumat, pool, wt = build(df, cfg)

    # run folder: timestamp + source filename, created alongside the input file
    # (or at --outdir)
    if cfg["outdir"]:
        run_dir = Path(cfg["outdir"]).expanduser().resolve()
    else:
        run_dir = path.resolve().parent / f"substitution_ranker_{RUN_STAMP}_{path.stem}"
    run_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{cfg['out_stem']}_{cfg['value']}_{cfg['style']}"

    print(f"[info] {len(labels)} substitutions | {fire_rule(cfg)} | value={cfg['value']}"
          f" | order={cfg['order']} | categories mutually exclusive (sum to pool)")
    print(f"[info] output dir: {run_dir}")
    print(f"       {'sub':>6} {'none':>5} {'4uM':>4} {'100uM':>5} {'pool':>5}  "
          f"sens%(4uM/pool)")
    rows = ([wt] if wt else []) + list(zip(labels, cmat, pool))
    for item in rows:
        l, c, p = item[0], item[1], item[2]
        print(f"       {l:>6} {int(c[0]):>5} {int(c[1]):>4} {int(c[2]):>5} {int(p):>5}"
              f"   {c[1] / p * 100:5.0f}%")

    # --- CSV of the count matrix (WT row first if present) ------------------
    csv_rows = ([wt[0]] if wt else []) + list(labels)
    csv_mat = np.vstack([wt[1], cmat]) if wt else cmat
    csv_pool = np.concatenate([[wt[2]], pool]) if wt else pool
    csv_df = pd.DataFrame(csv_mat, index=csv_rows, columns=list(cfg["concs"]))
    csv_df["pool"] = csv_pool
    csv_mu = np.vstack([wt[3], mumat]) if wt else mumat
    for k, c in enumerate(cfg["concs"]):
        csv_df[f"mumax_{c}"] = csv_mu[:, k]
    csv_path = run_dir / f"{stem}.csv"
    csv_df.to_csv(csv_path, index_label="substitution")
    print(f"[saved] {csv_path}")

    # --- config params dump -------------------------------------------------
    cfg_path = run_dir / f"{stem}_config.txt"
    with open(cfg_path, "w", encoding="utf-8") as fh:
        fh.write("substitution_ranker.py run parameters\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"timestamp : {RUN_STAMP}\n")
        fh.write(f"source    : {path.resolve()}\n")
        fh.write(f"gate rule : {fire_rule(cfg)}\n")
        fh.write(f"n_subs    : {len(labels)}\n\n")
        for k, v in cfg.items():
            fh.write(f"{k:>16} : {v}\n")
    print(f"[saved] {cfg_path}")

    # --- figures ------------------------------------------------------------
    fig = plot(labels, vmat, cmat, mumat, pool, wt, cfg)
    for fmt_ in cfg["formats"]:
        out = run_dir / f"{stem}.{fmt_}"
        fig.savefig(out, dpi=cfg["dpi"], bbox_inches="tight")
        print(f"[saved] {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
