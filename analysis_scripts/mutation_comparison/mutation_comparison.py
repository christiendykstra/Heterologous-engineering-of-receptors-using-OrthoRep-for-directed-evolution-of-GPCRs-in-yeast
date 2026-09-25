#!/usr/bin/env python3
"""
mutation_comparison.py
======================

Copyright (c) 2026 Christien Dykstra
SPDX-License-Identifier: MIT

Compare per-clone mutation counts between two conditions.

Splits clones into two groups by regex on a label column, then tests whether
mutation counts differ between them. Synonymous substitutions are the primary
readout,
since it is expected to be neutral and therefore reports on the underlying
mutation rate rather than on selection.

Tests performed:

    Mann-Whitney U          total / synonymous / nonsynonymous counts per clone
    rank-biserial r         effect size for each Mann-Whitney comparison
    Kolmogorov-Smirnov      synonymous count distributions
    Poisson rate ratio      exact conditional test (statsmodels, method
                            'exact-cond') with score-method 95% CI, on total
                            and synonymous event counts
    Fisher exact            dN/dS contingency between conditions
    minimum detectable RR   synonymous rate ratio detectable at 80% power,
                            alpha 0.05, given the observed event counts

The Poisson rate ratio and its CI are the basis for any claim that two
conditions have equivalent mutation rates; the minimum detectable rate ratio
states how large a difference the data could have resolved.


INPUT FORMAT
------------
One .xlsx or .csv, one row per clone. Required columns:

    barcode                                 label used for condition assignment
    NT_substitutions_count                  integer
    AA_substitutions_nonsynonymous_count    integer

Optional:

    AA_substitutions_synonymous   comma-separated; synonymous counts are taken
                                  from here when present. If absent, synonymous
                                  count falls back to
                                  (NT_count - nonsynonymous_count), clipped at 0.
    condition                     if present, used directly and the regex
                                  patterns are ignored
    Count                         read count, used only when --collapse is set

Conditions are assigned by matching --pattern-a and --pattern-b against the
`barcode` column. Each row must match exactly one pattern; rows matching both
or neither are reported and dropped.


OUTPUT
------
A timestamped folder next to the input (or --outdir) containing:

    mutation_comparison.svg / .pdf / .png
                                4-panel figure: three violin panels (total,
                                synonymous, nonsynonymous) and an ECDF of
                                synonymous substitutions per clone
    mutation_comparison_per_clone.csv   per-clone counts and assigned condition
    mutation_comparison_stats.csv       every statistic above, as a flat table


USAGE
-----
    python mutation_comparison.py data.xlsx
    python mutation_comparison.py data.xlsx --pattern-a 3AT --pattern-b "^SNC-80(?!.*3AT)" \
        --label-a "+3AT" --label-b "-3AT" --title "OPRD1 6 passages"
    python mutation_comparison.py data.xlsx --collapse --outdir results/

Requires: pandas numpy scipy matplotlib openpyxl statsmodels
(statsmodels is optional; without it the Poisson tests are skipped.)
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.stats import fisher_exact, ks_2samp, mannwhitneyu

try:
    from statsmodels.stats.rates import confint_poisson_2indep, test_poisson_2indep
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


# ============================================================================
# CONFIG — defaults; all of these are also command-line options
# ============================================================================

# Regex patterns splitting the label column into two groups.
# Each label must match exactly one pattern, or that row is dropped.
PATTERN_A = "3AT"                   # condition A: any row containing 3AT
PATTERN_B = r"^SNC-80(?!.*3AT)"     # condition B: SNC-80 rows WITHOUT 3AT

LABEL_A = "+3AT"                    # display label for condition A
LABEL_B = "-3AT"                    # display label for condition B

TITLE = "Mutation count comparison"

# Whether to collapse several haplotype rows sharing one barcode into a single
# read-weighted clone. Leave off when 'barcode' is a GROUP label and each row
# is already one clone; turn on only when 'barcode' is a real per-clone
# barcode that can carry several haplotype rows.
COLLAPSE_HAPLOTYPES = False

# Figure colours, condition A then B
COLOR_A = "#d95f02"
COLOR_B = "#1b9e77"

JITTER_SEED = 7                     # reproducible point jitter in the violins


# -------------------------- data loading --------------------------

def count_subs_from_string(s):
    """Count comma-separated substitutions in a string field. WT/NaN -> 0."""
    if pd.isna(s) or str(s).strip().upper() == "WT" or str(s).strip() == "":
        return 0
    return len([x for x in str(s).split(",") if x.strip() and x.strip().upper() != "WT"])


def load_and_label(path, pattern_a, pattern_b, label_a, label_b, collapse):
    path = str(path)
    if path.lower().endswith((".xlsx", ".xls", ".xlsm")):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    required = {"barcode", "NT_substitutions_count",
                "AA_substitutions_nonsynonymous_count"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: missing required columns: {sorted(missing)}\n"
                 f"       found: {list(df.columns)}")

    if "AA_substitutions_synonymous" in df.columns:
        df["syn_count"] = df["AA_substitutions_synonymous"].apply(count_subs_from_string)
    else:
        df["syn_count"] = (df["NT_substitutions_count"]
                           - df["AA_substitutions_nonsynonymous_count"]).clip(lower=0)
        print("NOTE: no AA_substitutions_synonymous column; using total-nonsyn subtraction.")

    df["nonsyn_count"] = df["AA_substitutions_nonsynonymous_count"]

    if "condition" in df.columns:
        print("NOTE: 'condition' column present; regex patterns ignored.")
    else:
        try:
            re_a, re_b = re.compile(pattern_a), re.compile(pattern_b)
        except re.error as exc:
            sys.exit(f"ERROR: invalid regex pattern: {exc}")

        def assign(bc):
            a_hit = bool(re_a.search(str(bc)))
            b_hit = bool(re_b.search(str(bc)))
            if a_hit and not b_hit:
                return label_a
            if b_hit and not a_hit:
                return label_b
            return None

        df["condition"] = df["barcode"].apply(assign)
        ambiguous = int(df["condition"].isna().sum())
        if ambiguous:
            unmatched = sorted(set(df.loc[df["condition"].isna(), "barcode"].astype(str)))[:8]
            print(f"WARNING: {ambiguous} rows matched neither / both patterns; dropping.")
            print(f"         e.g. {unmatched}")
            df = df.dropna(subset=["condition"])

    if df.empty:
        sys.exit("ERROR: no rows left after condition assignment. "
                 "Check --pattern-a / --pattern-b against the barcode column.")

    has_repeats = df.groupby("barcode").size().max() > 1
    if collapse and "Count" in df.columns and has_repeats:
        print("NOTE: collapsing multiple haplotypes per barcode by read-weighted mean.")

        def wmean(g):
            return pd.Series({
                "NT_substitutions_count": np.average(g["NT_substitutions_count"], weights=g["Count"]),
                "syn_count":               np.average(g["syn_count"],               weights=g["Count"]),
                "nonsyn_count":            np.average(g["nonsyn_count"],            weights=g["Count"]),
                "condition": g["condition"].iloc[0],
            })

        df = df.groupby("barcode").apply(wmean, include_groups=False).reset_index()
    elif collapse and has_repeats:
        print("NOTE: --collapse requested but no 'Count' column; "
              "treating each row as one clone.")
    elif has_repeats:
        print("NOTE: 'barcode' repeats but --collapse is off; "
              "treating each row as one clone (no read-weighting).")

    # A condition with no clones makes every downstream test meaningless.
    counts = df["condition"].value_counts()
    for lbl in (label_a, label_b):
        if counts.get(lbl, 0) == 0:
            sys.exit(f"ERROR: no clones assigned to condition '{lbl}'. "
                     f"Assigned: {counts.to_dict()}")

    return df


# -------------------------- statistics --------------------------

def run_stats(df, label_a, label_b):
    """Run all comparisons. Returns dict of results."""
    a = df[df.condition == label_a]
    b = df[df.condition == label_b]
    results = {"n_a": len(a), "n_b": len(b), "label_a": label_a, "label_b": label_b}

    for col, key in [("NT_substitutions_count", "total"),
                     ("syn_count", "syn"),
                     ("nonsyn_count", "nonsyn")]:
        results[f"mean_a_{key}"] = float(a[col].mean())
        results[f"mean_b_{key}"] = float(b[col].mean())
        results[f"median_a_{key}"] = float(a[col].median())
        results[f"median_b_{key}"] = float(b[col].median())
        results[f"sum_a_{key}"] = float(a[col].sum())
        results[f"sum_b_{key}"] = float(b[col].sum())
        U, p = mannwhitneyu(a[col], b[col], alternative="two-sided")
        results[f"mwu_U_{key}"] = float(U)
        results[f"mwu_p_{key}"] = float(p)
        if len(a) and len(b):
            results[f"rank_biserial_{key}"] = float(1 - (2 * U) / (len(a) * len(b)))

    D, p_ks = ks_2samp(a["syn_count"], b["syn_count"])
    results["ks_D_syn"] = float(D)
    results["ks_p_syn"] = float(p_ks)

    if HAS_STATSMODELS:
        for col, key in [("NT_substitutions_count", "total"), ("syn_count", "syn")]:
            ev_a, ev_b = int(round(a[col].sum())), int(round(b[col].sum()))
            if ev_a == 0 and ev_b == 0:
                continue
            res = test_poisson_2indep(ev_a, len(a), ev_b, len(b),
                                      method="exact-cond", value=1)
            ci = confint_poisson_2indep(ev_a, len(a), ev_b, len(b),
                                        method="score", compare="ratio")
            rr = (ev_a / len(a)) / (ev_b / len(b)) if ev_b else float("nan")
            results[f"poisson_RR_{key}"] = float(rr)
            results[f"poisson_RR_ci_lo_{key}"] = float(ci[0])
            results[f"poisson_RR_ci_hi_{key}"] = float(ci[1])
            results[f"poisson_p_{key}"] = float(res.pvalue)

    syn_a, ns_a = a["syn_count"].sum(), a["nonsyn_count"].sum()
    syn_b, ns_b = b["syn_count"].sum(), b["nonsyn_count"].sum()
    results["dnds_a"] = float(ns_a / syn_a) if syn_a else float("nan")
    results["dnds_b"] = float(ns_b / syn_b) if syn_b else float("nan")
    table = [[int(round(ns_a)), int(round(syn_a))],
             [int(round(ns_b)), int(round(syn_b))]]
    if all(v >= 0 for rowv in table for v in rowv) and sum(sum(r) for r in table) > 0:
        OR, p_f = fisher_exact(table)
        results["fisher_OR_dnds"] = float(OR)
        results["fisher_p_dnds"] = float(p_f)

    ev_a_syn = a["syn_count"].sum()
    ev_b_syn = b["syn_count"].sum()
    if ev_a_syn > 0 and ev_b_syn > 0:
        z = st.norm.ppf(0.975) + st.norm.ppf(0.80)
        mde_log = z * np.sqrt(1 / ev_a_syn + 1 / ev_b_syn)
        results["mde_rate_ratio_syn"] = float(np.exp(mde_log))

    return results


# -------------------------- plotting --------------------------

def make_figure(df, label_a, label_b, results, title, out_prefix):
    """Generate 4-panel figure: three violins + ECDF of syn counts."""
    mpl.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(1, 4, figsize=(12, 3.3), constrained_layout=True)
    colors = {label_a: COLOR_A, label_b: COLOR_B}
    order = [label_b, label_a]

    def violin_panel(ax, col, title_text, ylabel):
        data = [df[df.condition == c][col].values for c in order]
        parts = ax.violinplot(data, positions=[0, 1], widths=0.7,
                              showmeans=False, showmedians=False, showextrema=False)
        for pc, c in zip(parts["bodies"], order):
            pc.set_facecolor(colors[c])
            pc.set_alpha(0.35)
            pc.set_edgecolor("none")
        rng = np.random.default_rng(JITTER_SEED)
        for i, (vals, c) in enumerate(zip(data, order)):
            jitter = rng.uniform(-0.13, 0.13, size=len(vals))
            ax.scatter(i + jitter, vals, s=14, color=colors[c],
                       edgecolor="white", linewidth=0.4, alpha=0.85, zorder=3)
        for i, vals in enumerate(data):
            if len(vals) == 0:
                continue
            m = np.mean(vals)
            sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
            ax.errorbar(i, m, yerr=sem, fmt="_", color="black",
                        markersize=20, markeredgewidth=2, capsize=4, zorder=5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(order)
        ax.set_ylabel(ylabel)
        ax.set_title(title_text, fontsize=10)
        U, p = mannwhitneyu(data[0], data[1], alternative="two-sided")
        ymax = max(np.max(data[0]), np.max(data[1])) if (len(data[0]) and len(data[1])) else 1
        ax.set_ylim(-0.3, ymax * 1.30 + 0.3)
        ax.plot([0, 1], [ymax * 1.10] * 2, color="black", lw=0.8)
        ax.text(0.5, ymax * 1.14, f"p = {p:.3f}", ha="center", va="bottom", fontsize=8.5)

    violin_panel(axes[0], "NT_substitutions_count",
                 "Total NT substitutions", "NT substitutions / clone")
    violin_panel(axes[1], "syn_count",
                 "Synonymous", "Syn NT substitutions / clone")
    violin_panel(axes[2], "nonsyn_count",
                 "Non-synonymous", "Non-syn NT substitutions / clone")

    ax = axes[3]
    for c in order:
        vals = np.sort(df[df.condition == c]["syn_count"].values)
        if len(vals) == 0:
            continue
        y = np.arange(1, len(vals) + 1) / len(vals)
        x_plot = np.concatenate([[vals[0]], vals])
        y_plot = np.concatenate([[0], y])
        ax.step(x_plot, y_plot, where="post", color=colors[c], lw=2, label=c)
    ax.set_xlabel("Syn NT substitutions / clone")
    ax.set_ylabel("Cumulative fraction of clones")
    ax.set_title("Synonymous substitutions ECDF", fontsize=10)
    ax.legend(frameon=False, loc="lower right", fontsize=8.5)
    if "ks_D_syn" in results:
        ax.text(0.02, 0.97,
                f"K-S D = {results['ks_D_syn']:.3f}\np = {results['ks_p_syn']:.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=8.5)
    ax.set_ylim(-0.02, 1.05)

    n_a, n_b = results["n_a"], results["n_b"]
    fig.suptitle(f"{title}  (n = {n_b} {label_b}, n = {n_a} {label_a})",
                 fontsize=10.5, y=1.03)

    for ext in ("svg", "pdf", "png"):
        path = f"{out_prefix}.{ext}"
        fig.savefig(path, bbox_inches="tight",
                    dpi=200 if ext == "png" else None)
        print(f"  wrote: {path}")
    plt.close(fig)


# -------------------------- output tables --------------------------

def save_outputs(df, results, out_prefix):
    """Write per-clone CSV and stats summary CSV."""
    per_clone_path = f"{out_prefix}_per_clone.csv"
    cols = [c for c in ["barcode", "condition", "Count",
                        "NT_substitutions", "NT_substitutions_count",
                        "syn_count", "nonsyn_count",
                        "AA_substitutions_nonsynonymous",
                        "AA_substitutions_synonymous"]
            if c in df.columns]
    df[cols].to_csv(per_clone_path, index=False)
    print(f"  wrote: {per_clone_path}")

    la, lb = results["label_a"], results["label_b"]
    rows = [
        ("n clones", la, results["n_a"]),
        ("n clones", lb, results["n_b"]),
        ("Mean total NT/clone", la, f"{results['mean_a_total']:.3f}"),
        ("Mean total NT/clone", lb, f"{results['mean_b_total']:.3f}"),
        ("Mean syn/clone", la, f"{results['mean_a_syn']:.3f}"),
        ("Mean syn/clone", lb, f"{results['mean_b_syn']:.3f}"),
        ("Mean non-syn/clone", la, f"{results['mean_a_nonsyn']:.3f}"),
        ("Mean non-syn/clone", lb, f"{results['mean_b_nonsyn']:.3f}"),
        ("Total syn events", la, f"{results['sum_a_syn']:.0f}"),
        ("Total syn events", lb, f"{results['sum_b_syn']:.0f}"),
        ("Mann-Whitney U (syn)", "\u2014", f"{results['mwu_U_syn']:.0f}"),
        ("Mann-Whitney p (syn)", "\u2014", f"{results['mwu_p_syn']:.4f}"),
        ("Rank-biserial r (syn)", "\u2014", f"{results['rank_biserial_syn']:+.3f}"),
        ("K-S D (syn)", "\u2014", f"{results['ks_D_syn']:.3f}"),
        ("K-S p (syn)", "\u2014", f"{results['ks_p_syn']:.4f}"),
        ("Mann-Whitney p (total)", "\u2014", f"{results['mwu_p_total']:.4f}"),
        ("Mann-Whitney p (non-syn)", "\u2014", f"{results['mwu_p_nonsyn']:.4f}"),
        ("dN/dS", la, f"{results['dnds_a']:.3f}"),
        ("dN/dS", lb, f"{results['dnds_b']:.3f}"),
    ]
    if "poisson_RR_syn" in results:
        rows += [
            ("Poisson RR (syn, A/B)", "\u2014", f"{results['poisson_RR_syn']:.3f}"),
            ("Poisson RR 95% CI (syn)", "\u2014",
             f"[{results['poisson_RR_ci_lo_syn']:.3f}, {results['poisson_RR_ci_hi_syn']:.3f}]"),
            ("Poisson exact p (syn)", "\u2014", f"{results['poisson_p_syn']:.4f}"),
        ]
    if "poisson_RR_total" in results:
        rows += [
            ("Poisson RR (total, A/B)", "\u2014", f"{results['poisson_RR_total']:.3f}"),
            ("Poisson RR 95% CI (total)", "\u2014",
             f"[{results['poisson_RR_ci_lo_total']:.3f}, {results['poisson_RR_ci_hi_total']:.3f}]"),
            ("Poisson exact p (total)", "\u2014", f"{results['poisson_p_total']:.4f}"),
        ]
    if "fisher_OR_dnds" in results:
        rows += [
            ("Fisher OR (dN/dS A vs B)", "\u2014", f"{results['fisher_OR_dnds']:.3f}"),
            ("Fisher p (dN/dS)", "\u2014", f"{results['fisher_p_dnds']:.4f}"),
        ]
    if "mde_rate_ratio_syn" in results:
        rows.append(("Min detectable RR (syn, 80% power)", "\u2014",
                     f"{results['mde_rate_ratio_syn']:.2f}x"))
    if not HAS_STATSMODELS:
        rows.append(("Poisson tests", "\u2014", "skipped (statsmodels not installed)"))

    stats_path = f"{out_prefix}_stats.csv"
    pd.DataFrame(rows, columns=["Statistic", "Group", "Value"]).to_csv(stats_path, index=False)
    print(f"  wrote: {stats_path}")


# -------------------------- cli --------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compare per-clone mutation counts between two conditions "
                    "(Mann-Whitney, K-S, Poisson rate ratio, dN/dS Fisher).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Input: one .xlsx/.csv, one row per clone, with barcode, "
               "NT_substitutions_count and AA_substitutions_nonsynonymous_count. "
               "See the module docstring for the full format.")
    p.add_argument("input", help="Input .xlsx or .csv file")
    p.add_argument("-o", "--outdir", default=None,
                   help="Output directory (default: timestamped folder next to the input)")
    p.add_argument("--pattern-a", default=PATTERN_A,
                   help=f"Regex selecting condition A (default: {PATTERN_A!r})")
    p.add_argument("--pattern-b", default=PATTERN_B,
                   help=f"Regex selecting condition B (default: {PATTERN_B!r})")
    p.add_argument("--label-a", default=LABEL_A,
                   help=f"Display label for condition A (default: {LABEL_A!r})")
    p.add_argument("--label-b", default=LABEL_B,
                   help=f"Display label for condition B (default: {LABEL_B!r})")
    p.add_argument("--title", default=TITLE, help="Figure title")
    p.add_argument("--collapse", dest="collapse", action="store_true",
                   default=COLLAPSE_HAPLOTYPES,
                   help="Collapse haplotype rows sharing a barcode by read-weighted mean")
    return p.parse_args()


# -------------------------- main --------------------------

def main():
    args = parse_args()

    csv_path = Path(args.input).expanduser().resolve()
    if not csv_path.exists():
        sys.exit(f"ERROR: input file not found: {csv_path}")

    if not HAS_STATSMODELS:
        print("WARNING: statsmodels not installed. Poisson rate-ratio tests will "
              "be skipped.\n         Install with:  pip install statsmodels")

    if args.outdir:
        out_dir = Path(args.outdir).expanduser().resolve()
    else:
        ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out_dir = csv_path.parent / f"mutation_comparison_{csv_path.stem}_{ts}"

    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = str(out_dir / "mutation_comparison")

    print(f"Loading {csv_path}...")
    df = load_and_label(csv_path, args.pattern_a, args.pattern_b,
                        args.label_a, args.label_b, args.collapse)
    print(f"  {len(df)} clones loaded.")
    print(f"  conditions: {df['condition'].value_counts().to_dict()}")

    print("\nRunning statistics...")
    results = run_stats(df, args.label_a, args.label_b)

    print(f"\n  A ({args.label_a}, n={results['n_a']}):  "
          f"total={results['mean_a_total']:.2f}  syn={results['mean_a_syn']:.2f}  "
          f"nonsyn={results['mean_a_nonsyn']:.2f}")
    print(f"  B ({args.label_b}, n={results['n_b']}):  "
          f"total={results['mean_b_total']:.2f}  syn={results['mean_b_syn']:.2f}  "
          f"nonsyn={results['mean_b_nonsyn']:.2f}")
    print(f"\n  Mann-Whitney p (syn):    {results['mwu_p_syn']:.4f}")
    print(f"  K-S p (syn):             {results['ks_p_syn']:.4f}")
    if "poisson_RR_syn" in results:
        print(f"  Poisson RR (syn):        {results['poisson_RR_syn']:.3f}  "
              f"[{results['poisson_RR_ci_lo_syn']:.3f}, {results['poisson_RR_ci_hi_syn']:.3f}]  "
              f"p={results['poisson_p_syn']:.4f}")
    print(f"  dN/dS:                   {args.label_a}={results['dnds_a']:.3f}  "
          f"{args.label_b}={results['dnds_b']:.3f}")
    if "mde_rate_ratio_syn" in results:
        print(f"  Min detectable RR (syn): {results['mde_rate_ratio_syn']:.2f}x")

    print(f"\nWriting outputs to:\n  {out_dir}")
    make_figure(df, args.label_a, args.label_b, results, args.title, out_prefix)
    save_outputs(df, results, out_prefix)
    print("\nDone.")


if __name__ == "__main__":
    main()
