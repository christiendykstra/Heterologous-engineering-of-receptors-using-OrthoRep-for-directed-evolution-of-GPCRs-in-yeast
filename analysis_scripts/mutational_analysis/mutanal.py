#!/usr/bin/env python3
"""
mutanal.py
==========

Copyright (c) 2026 Christien Dykstra
SPDX-License-Identifier: MIT

Mutation analysis for nanopore-sequenced directed-evolution populations demultiplexed by pipeline.sh.

Parses per-clone nucleotide and amino-acid substitution calls and produces:

  * mutation frequency bar chart (with optional y-axis break and WT reference)
  * mutation burden histograms (AA and NT, mutations per clone)
  * position distribution histogram
  * mutation complexity distribution
  * substitution matrix heatmap (AA 20x20 or NT)
  * dedicated 4x4 NT substitution matrix annotated with the Ts/Tv ratio
  * gene map of nonsynonymous and synonymous mutation density along the ORF
  * a text report and a dump of every parameter used for the run


INPUT FORMAT
------------
One or more .xlsx or .csv files, each with one row per sequenced clone.
Columns are read by name; missing optional columns degrade gracefully with a
warning rather than failing.

    AA_substitutions_nonsynonymous    comma-separated, e.g. "F116S, T189A"
    AA_substitutions_synonymous       comma-separated (optional; enables the
                                      WT-synonymous bar and the gene map's
                                      synonymous track)
    NT_substitutions                  comma-separated, e.g. "A347G, C912T"

    AA_substitutions_nonsynonymous_count   integer (optional)
    NT_substitutions_count                 integer (optional)

The *_count columns are used for the burden histograms when present; otherwise
the substitution strings are parsed and counted. A clone with empty/NA
substitution fields counts as wild-type.

Substitutions must be <original><position><new> with a single-character
original and new (e.g. A51T, or A347G for NT). Entries that don't parse are
skipped silently.

Multiple input files are concatenated, so a run can pool several plates or
timepoints. A consolidated copy of the merged table is written to the run
folder.


OUTPUT
------
A timestamped subfolder (default `mutstats_<date>_<time>/`) created under
--outdir, containing every figure in each requested format, the consolidated
CSV, per-type text reports, and run_parameters.txt recording the full config.


USAGE
-----
    python mutanal.py P10.xlsx
    python mutanal.py P5.xlsx P10.xlsx --outdir results/
    python mutanal.py data.csv --mutation-type BOTH
    python mutanal.py P10.xlsx --formats png --dpi 150

Every constant in the CONFIG sections below can be edited directly; the most
commonly changed ones are also exposed as command-line flags (see --help).

Requires: pandas numpy matplotlib seaborn openpyxl
"""

import argparse
import os
import random
import sys
from collections import Counter
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, PowerNorm
from matplotlib.patches import Rectangle

# ============================================================================
# 1. OUTPUT
# ============================================================================
# Input files and output directory are given on the command line; see --help.

CSV_FILENAME = 'consolidated_input.csv'   # merged input written into the run folder

# Each run gets its own timestamped subfolder (set to None to write directly
# into the output directory). A run_parameters.txt is written into this folder.
OUTPUT_SUBFOLDER = 'mutstats'          # -> "mutstats_2026-04-06_14-30-00/"

# What to analyze and how to export.
MUTATION_TYPE = 'AA'                   # 'AA', 'NT', or 'BOTH'
OUTPUT_FORMATS = ['png', 'pdf', 'svg']
DPI = 300
SVG_FONTTYPE = 'none'                  # 'none' = editable text, 'path' = outlines
PDF_FONTTYPE = 42                      # 42 = editable TrueType, 3 = Type3 outlines

# ============================================================================
# 2. GLOBAL COLORS & FONTS
# ============================================================================

COLOR_SINGLE = '#2E86AB'               # default single-color bars
COLOR_WT = '#808080'                   # WT bars
COLOR_WT_SYN = '#A9A9A9'               # WT-synonymous bars
COLOR_BURDEN_AA = '#2E86AB'            # AA burden histogram
COLOR_BURDEN_NT = '#E85D04'            # NT burden histogram

FONTSIZE_TITLE = 6
FONTSIZE_LABELS = 4
FONTSIZE_TICKS = 4
FONTSIZE_LEGEND = 4

# Resulting-amino-acid palette (used when COLOR_BY_AMINO_ACID = True)
AA_COLOR_PALETTE = {
    'A': '#F84514', 'C': '#79B1DC', 'D': '#7143CC', 'E': '#8F8BD9',
    'F': '#F93D4A', 'G': '#D43206', 'H': '#FFBE92', 'I': '#E89F88',
    'K': '#FA9A5A', 'L': '#FB9275', 'M': '#D89272', 'N': '#9DD3E4',
    'P': '#B46546', 'Q': '#B7EBEA', 'R': '#F97B27', 'S': '#054CB7',
    'T': '#4281CB', 'V': '#F96239', 'W': '#BD575E', 'Y': '#6E2987',
}

# ============================================================================
# 3. FIGURE SIZES
# ============================================================================

FIGSIZE_MUTATION_FREQUENCY = (3, 3)
FIGSIZE_POSITION_DISTRIBUTION = (10, 10)
FIGSIZE_COMPLEXITY = (3, 3)
FIGSIZE_HEATMAP = (7, 5)               # AA / general substitution matrix
FIGSIZE_NT_HEATMAP = (5, 5)            # dedicated 4x4 NT matrix
FIGSIZE_BURDEN = (3, 3)
GENE_MAP_FIGSIZE = (14, 4)

# ============================================================================
# 4. BAR CHART APPEARANCE
# ============================================================================

ALPHA_BARS = 1
EDGECOLOR_BARS = 'black'
LINEWIDTH_BARS = 0
BAR_WIDTH = 1.0                        # 1.0 = no gaps

YLIM_MUTATION_FREQUENCY = (0, 0.3)     # None = auto
YLIM_POSITION_DISTRIBUTION = None
YLIM_COMPLEXITY = None
YLIM_BURDEN = None

# Optional y-axis break on the frequency plot
USE_AXIS_BREAK = True
AXIS_BREAK_LOWER_RANGE = (0, 0.05)
AXIS_BREAK_UPPER_RANGE = (0.70, 1)
AXIS_BREAK_LOWER_HEIGHT = 0.60
AXIS_BREAK_UPPER_HEIGHT = 0.40

# ============================================================================
# 5. MUTATION FREQUENCY BAR CHART
# ============================================================================

# Limit how many mutation bars are shown. Of DISPLAY_N_MUTATIONS, the first
# DISPLAY_TOP_ABUNDANT are the most frequent; the rest are randomly sampled
# from the long tail to show diversity. Set DISPLAY_N_MUTATIONS = None for all.
DISPLAY_N_MUTATIONS = 20
DISPLAY_TOP_ABUNDANT = 10
DISPLAY_RANDOM_SEED = 42               # None = non-reproducible random fill

WT_GAP_WIDTH = 1.5                     # visual gap between WT group and mutations
COLOR_BY_AMINO_ACID = False            # color AA bars by resulting amino acid
SHOW_NUMBERS_ON_BARS = False

# ============================================================================
# 6. MUTATION BURDEN HISTOGRAMS
# ============================================================================

GENERATE_BURDEN_HISTOGRAMS = True
BURDEN_MAX_MUTATIONS_AA = 8            # last bin is "N+"
BURDEN_MAX_MUTATIONS_NT = 11
BURDEN_SHOW_PROPORTION = True          # y as proportion of sequences vs raw count
SHOW_NUMBERS_ON_BURDEN = True
BURDEN_BAR_WIDTH = 0.8

# ============================================================================
# 7. COMPLEXITY
# ============================================================================

SHOW_NUMBERS_ON_COMPLEXITY = True

# ============================================================================
# 8. SUBSTITUTION HEATMAP (AA 20x20 / general)  ---  gradient output
# ============================================================================
# Applies to plot_substitution_heatmap. Color precedence is:
#     GAMMA (PowerNorm)  >  LOG (LogNorm)  >  linear (VMIN/VMAX)
# i.e. if HEATMAP_GAMMA is set it wins; log is ignored. Log scale blanks
# zero cells (they render white with no annotation).

COLORMAP_HEATMAP = 'magma_r'

# Cell value scaling:
#   'count'        -> raw substitution counts
#   'proportion'   -> fraction of all substitutions (cells sum to 1)
#   'per_sequence' -> count / n_sequences
#   'per_residue'  -> count / (n_sequences * length)   [uses HEATMAP_LENGTH_*]
HEATMAP_NORMALIZE = 'per_sequence'

HEATMAP_GAMMA = 0.45                   # PowerNorm gamma (<1 brightens low cells); None to disable
HEATMAP_LOG_SCALE = False              # only used when HEATMAP_GAMMA is None
HEATMAP_VMIN = None                    # None = auto from data
HEATMAP_VMAX = None

HEATMAP_LENGTH_AA = 400                # residues, used only by 'per_residue'
HEATMAP_LENGTH_NT = 1200               # bp, used only by 'per_residue' (NT)

SHOW_NUMBERS_ON_HEATMAP = False
REMOVE_EMPTY_HEATMAP_ROWS = True       # drop originals that never mutate
REMOVE_EMPTY_HEATMAP_COLS = True       # drop results that never appear

# ============================================================================
# 9. DEDICATED 4x4 NT SUBSTITUTION HEATMAP (Ts/Tv)  ---  gradient output
# ============================================================================
# Always shows the full 4x4 ACGT grid (empty cells kept), can outline
# transitions/transversions, and annotates the Ts/Tv ratio. Runs regardless
# of MUTATION_TYPE. Same color precedence as section 8.

GENERATE_NT_SUBSTITUTION_HEATMAP = True
COLORMAP_NT_HEATMAP = 'magma_r'
SHOW_NUMBERS_ON_NT_HEATMAP = True

# Cell value scaling: 'count' | 'proportion' | 'per_sequence' | 'per_base'
NT_HEATMAP_NORMALIZE = 'proportion'
NT_HEATMAP_GENE_LENGTH_BP = 1200       # used only by 'per_base'

NT_HEATMAP_GAMMA = None                # PowerNorm gamma; None to disable
NT_HEATMAP_LOG_SCALE = True            # only used when NT_HEATMAP_GAMMA is None
NT_HEATMAP_VMIN = 0.001                # None = auto
NT_HEATMAP_VMAX = 1

# Highlight transition cells (A<->G, C<->T)
NT_HEATMAP_HIGHLIGHT_TRANSITIONS = False
NT_HEATMAP_TRANSITION_EDGECOLOR = '#00A651'
NT_HEATMAP_TRANSITION_LINEWIDTH = 3
# Highlight transversion cells (all other off-diagonal)
NT_HEATMAP_HIGHLIGHT_TRANSVERSIONS = False
NT_HEATMAP_TRANSVERSION_EDGECOLOR = '#E63946'
NT_HEATMAP_TRANSVERSION_LINEWIDTH = 3

NT_HEATMAP_SHOW_TSTV_RATIO = True

# ============================================================================
# 10. GENE MAP
# ============================================================================
# NOTE: the gene map deliberately uses two *semantic* gradients (nonsynonymous
# yellow->red, synonymous gray->black, white = zero) so the two tracks stay
# distinguishable. It is intentionally NOT switched to magma_r; change the
# anchor colors below if you want a different scheme.

GENERATE_GENE_MAP = True
GENE_MAP_BIN_SIZE = 1                   # amino acids per bin (1 = per-position)
GENE_MAP_GENE_LENGTH_AA = 400
GENE_MAP_GENE_LENGTH_NT = 1200
GENE_MAP_SPACING = 1
GENE_MAP_MAX_FREQUENCY = 0.80           # colorbar max (None = auto)

GENE_MAP_NONSYN_ZERO_COLOR = '#FFFFFF'
GENE_MAP_NONSYN_LOW_COLOR = '#ffde17'
GENE_MAP_NONSYN_HIGH_COLOR = '#d23727'
GENE_MAP_SYN_ZERO_COLOR = '#FFFFFF'
GENE_MAP_SYN_LOW_COLOR = '#a7a9ac'
GENE_MAP_SYN_HIGH_COLOR = '#000000'


# ============================================================================
# SHARED HELPERS
# ============================================================================

def build_color_norm(display_matrix, *, gamma=None, log_scale=False,
                     vmin=None, vmax=None):
    """
    Return (data_for_plot, norm, label_prefix) for a heatmap.

    Color precedence: gamma (PowerNorm) > log (LogNorm) > linear (Normalize).
    Log scale blanks zero cells (NaN) so they render empty with no annotation.
    """
    data = display_matrix.astype(float).copy()
    prefix = ''

    if gamma is not None:
        norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
    elif log_scale:
        data = data.replace(0, np.nan)
        finite = data.values[np.isfinite(data.values)]
        lo = vmin if vmin is not None else (finite.min() if finite.size else None)
        hi = vmax if vmax is not None else (finite.max() if finite.size else None)
        norm = LogNorm(vmin=lo, vmax=hi)
        prefix = 'log '
    elif vmin is not None or vmax is not None:
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = None

    return data, norm, prefix


def write_run_parameters(output_dir, input_files=None):
    """
    Write every ALL-CAPS config constant to run_parameters.txt in the run
    folder, so each figure batch is reproducible from its own record.
    """
    path = os.path.join(output_dir, 'run_parameters.txt')
    g = globals()

    with open(path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("RUN PARAMETERS\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        if input_files:
            f.write(f"Input files: {', '.join(str(p) for p in input_files)}\n")
        f.write(f"Output folder: {output_dir}\n")
        f.write("=" * 70 + "\n\n")

        for name in sorted(g):
            if not name.isupper():
                continue
            val = g[name]
            if isinstance(val, (str, int, float, bool, tuple, list, dict)) or val is None:
                f.write(f"{name} = {val!r}\n")

    print(f"Run parameters saved: {path}")
    return path


# ============================================================================
# DATA LOADING & PARSING
# ============================================================================

def load_input_files(paths):
    """Load and concatenate the given .csv / .xlsx files."""
    dfs = []
    for filepath in paths:
        if not os.path.isfile(filepath):
            sys.exit(f"[error] input file not found: {filepath}")
        ext = os.path.splitext(filepath)[1].lower()
        if ext == '.csv':
            df = pd.read_csv(filepath)
        elif ext in ('.xlsx', '.xls', '.xlsm'):
            df = pd.read_excel(filepath)
        else:
            sys.exit(f"[error] unsupported file type '{ext}': {filepath}")
        dfs.append(df)
        print(f"  Loaded: {os.path.basename(filepath)} ({len(df)} rows)")

    combined = pd.concat(dfs, ignore_index=True)
    if len(paths) > 1:
        print(f"  Combined: {len(combined)} rows from {len(paths)} files")
    return combined


def parse_mutations(df, mutation_type='AA'):
    """
    Parse comma-separated mutation strings (e.g. 'F116S, T189A') into
    one row per individual mutation with original/position/new fields.
    """
    if mutation_type == 'AA':
        mutation_col = 'AA_substitutions_nonsynonymous'
    elif mutation_type == 'NT':
        mutation_col = 'NT_substitutions'
    else:
        raise ValueError("mutation_type must be 'AA' or 'NT'")

    if mutation_col not in df.columns:
        print(f"Warning: Column '{mutation_col}' not found. "
              f"Available columns: {list(df.columns)}")
        return pd.DataFrame()

    print(f"Parsing {mutation_type} mutations from column: '{mutation_col}'")

    mutations = []
    for idx, mut_string in enumerate(df[mutation_col]):
        if pd.isna(mut_string):
            continue
        for mut in [m.strip() for m in str(mut_string).split(',')]:
            if len(mut) < 3:
                continue
            try:
                position = int(mut[1:-1])
            except ValueError:
                continue
            mutations.append({
                'mutation': mut,
                'original': mut[0],
                'position': position,
                'new': mut[-1],
                'sequence_id': idx,
            })

    result_df = pd.DataFrame(mutations)
    print(f"  Parsed {len(result_df)} total {mutation_type} mutations "
          f"from {len(df)} sequences")
    return result_df


def calculate_statistics(df, total_sequences):
    """Summary statistics from parsed mutation data."""
    if df.empty:
        return {
            'total_mutations': 0,
            'unique_mutations': 0,
            'unique_positions': 0,
            'mean_mutations_per_sequence': 0,
            'total_sequences': total_sequences,
        }

    has_seq = 'sequence_id' in df.columns
    n_seq_with_muts = df['sequence_id'].nunique() if has_seq else 0
    return {
        'total_mutations': len(df),
        'unique_mutations': df['mutation'].nunique(),
        'unique_positions': df['position'].nunique(),
        'unique_sequences_with_mutations': n_seq_with_muts,
        'mean_mutations_per_sequence': (len(df) / n_seq_with_muts
                                        if n_seq_with_muts > 0 else 0),
        'total_sequences': total_sequences,
    }


def count_wt_sequences(df, mutation_type='AA'):
    """Count WT (no mutations) and WT-synonymous sequences."""
    if mutation_type == 'AA':
        nonsyn_col = 'AA_substitutions_nonsynonymous'
        syn_col = 'AA_substitutions_synonymous'
    else:
        nonsyn_col = 'NT_substitutions'
        syn_col = None

    if nonsyn_col not in df.columns:
        print(f"Warning: '{nonsyn_col}' not found, cannot count WT sequences")
        return 0, 0

    if syn_col and syn_col in df.columns:
        wt_count = len(df[(df[nonsyn_col].isna() | (df[nonsyn_col] == '')) &
                          (df[syn_col].isna() | (df[syn_col] == ''))])
        wt_syn_count = len(df[(df[nonsyn_col].isna() | (df[nonsyn_col] == '')) &
                              (df[syn_col].notna() & (df[syn_col] != ''))])
    else:
        wt_count = len(df[df[nonsyn_col].isna() | (df[nonsyn_col] == '')])
        wt_syn_count = 0

    return wt_count, wt_syn_count


def count_mutations_per_sequence(df, count_col=None, mutation_col=None):
    """Mutations per sequence, preferring a pre-computed count column."""
    if count_col and count_col in df.columns:
        print(f"    Using pre-computed count column: '{count_col}'")
        return pd.to_numeric(df[count_col], errors='coerce').fillna(0).astype(int).values

    if mutation_col and mutation_col in df.columns:
        print(f"    Warning: Count column '{count_col}' not found, "
              f"falling back to parsing '{mutation_col}'")
        parsed_counts = []
        for mut_string in df[mutation_col]:
            if pd.isna(mut_string) or str(mut_string).strip() == '':
                parsed_counts.append(0)
            else:
                muts = [m.strip() for m in str(mut_string).split(',') if m.strip()]
                parsed_counts.append(len([m for m in muts if len(m) >= 3]))
        return np.array(parsed_counts)

    print(f"    Warning: Neither '{count_col}' nor '{mutation_col}' found")
    return np.zeros(len(df), dtype=int)


# ============================================================================
# PLOTTING
# ============================================================================

def plot_mutation_burden(df, output_dir, mutation_type='AA'):
    """Histogram of mutation burden per sequence (mutations per clone)."""
    if not GENERATE_BURDEN_HISTOGRAMS:
        return

    if mutation_type == 'AA':
        count_col = 'AA_substitutions_nonsynonymous_count'
        mutation_col = 'AA_substitutions_nonsynonymous'
        max_muts = BURDEN_MAX_MUTATIONS_AA
        bar_color = COLOR_BURDEN_AA
        label = 'Nonsynonymous AA'
    else:
        count_col = 'NT_substitutions_count'
        mutation_col = 'NT_substitutions'
        max_muts = BURDEN_MAX_MUTATIONS_NT
        bar_color = COLOR_BURDEN_NT
        label = 'NT'

    if count_col not in df.columns and mutation_col not in df.columns:
        print(f"  Warning: Neither '{count_col}' nor '{mutation_col}' found, "
              f"skipping {mutation_type} burden histogram")
        return

    print(f"  - Plotting {mutation_type} mutation burden histogram...")

    mut_counts = count_mutations_per_sequence(df, count_col=count_col,
                                              mutation_col=mutation_col)
    total_sequences = len(mut_counts)

    bins = list(range(0, max_muts + 1))
    bin_counts = np.array([
        np.sum(mut_counts == b) if b < max_muts else np.sum(mut_counts >= b)
        for b in bins
    ])
    x_labels = [str(b) for b in bins[:-1]] + [f'{max_muts}+']

    if BURDEN_SHOW_PROPORTION:
        y_values = bin_counts / total_sequences
        y_label = 'Proportion of sequences'
    else:
        y_values = bin_counts.astype(float)
        y_label = 'Number of sequences'

    fig, ax = plt.subplots(figsize=FIGSIZE_BURDEN)
    bars = ax.bar(range(len(bins)), y_values, width=BURDEN_BAR_WIDTH,
                  color=bar_color, alpha=ALPHA_BARS,
                  edgecolor=EDGECOLOR_BARS, linewidth=0.5)

    if SHOW_NUMBERS_ON_BURDEN:
        for bar, yv, count in zip(bars, y_values, bin_counts):
            if count > 0:
                txt = f'{yv:.2f}\n({int(count)})' if BURDEN_SHOW_PROPORTION else f'{int(count)}'
                ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                        txt, ha='center', va='bottom', fontsize=FONTSIZE_TICKS - 2)

    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(x_labels, fontsize=FONTSIZE_TICKS)
    ax.set_xlabel(f'Number of {label} mutations per sequence', fontsize=FONTSIZE_LABELS)
    ax.set_ylabel(y_label, fontsize=FONTSIZE_LABELS)
    ax.set_title(f'{label} Mutation Burden Distribution',
                 fontsize=FONTSIZE_TITLE, fontweight='bold')
    ax.tick_params(axis='y', labelsize=FONTSIZE_TICKS)

    ax.text(0.97, 0.95,
            f'N = {total_sequences}\nMean = {np.mean(mut_counts):.1f}\n'
            f'Median = {np.median(mut_counts):.0f}',
            transform=ax.transAxes, fontsize=FONTSIZE_LABELS - 1,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.8))

    if YLIM_BURDEN is not None:
        ax.set_ylim(YLIM_BURDEN)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'{mutation_type.lower()}_mutation_burden.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_mutation_frequency(mutation_counts, output_dir, title, filename_prefix,
                            color_by_aa=False, display_n=None, top_abundant=None,
                            show_numbers=True, total_sequences=None,
                            wt_count=0, wt_syn_count=0):
    """
    Bar chart of mutation frequencies (proportion of total sequences), sorted
    by position. When display_n is set, shows the top_abundant most frequent
    plus a random long-tail sample, with a visual gap after the WT group.
    """
    if not mutation_counts and wt_count == 0 and wt_syn_count == 0:
        print(f"No data for {title}")
        return
    if total_sequences is None or total_sequences == 0:
        print("Warning: total_sequences not provided or zero, cannot calculate proportions")
        return

    def extract_position(mutation):
        try:
            return int(mutation[1:-1])
        except (ValueError, IndexError):
            return 0

    # --- selection: top abundant + random fill ---
    all_mutations = list(mutation_counts.items())
    if display_n is not None and len(all_mutations) > display_n:
        if top_abundant is None:
            top_abundant = display_n
        top_abundant = min(top_abundant, display_n, len(all_mutations))

        by_freq = sorted(all_mutations, key=lambda x: x[1], reverse=True)
        winners = by_freq[:top_abundant]
        winner_names = {m[0] for m in winners}
        remaining = [m for m in all_mutations if m[0] not in winner_names]
        n_random = min(display_n - top_abundant, len(remaining))

        rng = random.Random(DISPLAY_RANDOM_SEED)
        random_picks = rng.sample(remaining, n_random) if n_random > 0 else []

        selected = sorted(winners + random_picks, key=lambda x: extract_position(x[0]))
        subtitle_note = (f" (top {top_abundant} + {n_random} random, "
                         f"{len(all_mutations)} total unique)")
    else:
        selected = sorted(all_mutations, key=lambda x: extract_position(x[0]))
        subtitle_note = ""

    mutations = [x[0] for x in selected]
    counts = [x[1] for x in selected]

    # --- build WT group + gap + mutation x-positions ---
    wt_labels, wt_counts_list = [], []
    if wt_count > 0:
        wt_labels.append('WT'); wt_counts_list.append(wt_count)
    if wt_syn_count > 0:
        wt_labels.append('WT-syn'); wt_counts_list.append(wt_syn_count)
    n_wt = len(wt_labels)

    x_positions, all_labels, all_counts, all_proportions, colors_list = [], [], [], [], []

    for i, (lbl, cnt) in enumerate(zip(wt_labels, wt_counts_list)):
        x_positions.append(i)
        all_labels.append(lbl)
        all_counts.append(cnt)
        all_proportions.append(cnt / total_sequences)
        colors_list.append(COLOR_WT if lbl == 'WT' else COLOR_WT_SYN)

    gap = WT_GAP_WIDTH if n_wt > 0 else 0
    mut_start = n_wt + gap
    for i, (mut, cnt) in enumerate(zip(mutations, counts)):
        x_positions.append(mut_start + i)
        all_counts.append(cnt)
        all_proportions.append(cnt / total_sequences)
        if color_by_aa and len(mut) >= 3 and mut[-1].isalpha():
            colors_list.append(AA_COLOR_PALETTE.get(mut[-1], COLOR_SINGLE))
            all_labels.append(mut[:-1])
        else:
            colors_list.append(COLOR_SINGLE)
            all_labels.append(mut)

    full_title = title + subtitle_note

    def add_aa_legend(target_ax):
        if not color_by_aa:
            return
        aa_mutations = [m for m in mutations if len(m) >= 3 and m[-1].isalpha()]
        if not aa_mutations:
            return
        legend_elements = []
        if wt_count > 0:
            legend_elements.append(plt.Rectangle((0, 0), 1, 1, fc=COLOR_WT,
                                                 edgecolor='black', label='WT'))
        if wt_syn_count > 0:
            legend_elements.append(plt.Rectangle((0, 0), 1, 1, fc=COLOR_WT_SYN,
                                                 edgecolor='black', label='WT-syn'))
        legend_elements.extend(
            plt.Rectangle((0, 0), 1, 1, fc=AA_COLOR_PALETTE.get(aa, COLOR_SINGLE),
                          edgecolor='black', label=aa)
            for aa in sorted(set(m[-1] for m in aa_mutations)))
        target_ax.legend(handles=legend_elements, title='Resulting AA',
                         loc='upper right', fontsize=FONTSIZE_LEGEND)

    if USE_AXIS_BREAK:
        fig = plt.figure(figsize=FIGSIZE_MUTATION_FREQUENCY)
        gs = fig.add_gridspec(2, 1,
                              height_ratios=[AXIS_BREAK_UPPER_HEIGHT, AXIS_BREAK_LOWER_HEIGHT],
                              hspace=0.05)
        ax_upper = fig.add_subplot(gs[0])
        ax_lower = fig.add_subplot(gs[1], sharex=ax_upper)

        for ax in (ax_upper, ax_lower):
            ax.bar(x_positions, all_proportions, width=BAR_WIDTH, color=colors_list,
                   alpha=ALPHA_BARS, edgecolor=EDGECOLOR_BARS, linewidth=LINEWIDTH_BARS)
        ax_upper.set_ylim(AXIS_BREAK_UPPER_RANGE)
        ax_lower.set_ylim(AXIS_BREAK_LOWER_RANGE)
        ax_upper.spines['bottom'].set_visible(False)
        ax_lower.spines['top'].set_visible(False)
        ax_upper.xaxis.tick_top()
        ax_upper.tick_params(labeltop=False, top=False, bottom=False)
        ax_lower.xaxis.tick_bottom()

        d = 0.01
        kwargs = dict(color='k', clip_on=False, linewidth=1.5)
        ax_upper.plot([0, d], [0, -d], transform=ax_upper.transAxes, **kwargs)
        ax_upper.plot([1 - d, 1], [-d, 0], transform=ax_upper.transAxes, **kwargs)
        ax_lower.plot([0, d], [1 + d, 1], transform=ax_lower.transAxes, **kwargs)
        ax_lower.plot([1 - d, 1], [1, 1 + d], transform=ax_lower.transAxes, **kwargs)

        if show_numbers:
            for xp, proportion, count in zip(x_positions, all_proportions, all_counts):
                if AXIS_BREAK_LOWER_RANGE[0] <= proportion <= AXIS_BREAK_LOWER_RANGE[1]:
                    ax = ax_lower
                elif AXIS_BREAK_UPPER_RANGE[0] <= proportion <= AXIS_BREAK_UPPER_RANGE[1]:
                    ax = ax_upper
                else:
                    continue
                ax.text(xp, proportion, f'{proportion:.3f}\n({int(count)})',
                        ha='center', va='bottom', fontsize=FONTSIZE_TICKS - 2)

        add_aa_legend(ax_upper)
        ax_lower.set_xlabel("Mutation (ordered N\u2192C or 5'\u21923')", fontsize=FONTSIZE_LABELS)
        fig.text(0.04, 0.5, 'Frequency (proportion of total sequences)',
                 va='center', rotation='vertical', fontsize=FONTSIZE_LABELS)
        ax_upper.set_title(full_title, fontsize=FONTSIZE_TITLE, fontweight='bold', pad=20)
        ax_lower.set_xticks(x_positions)
        ax_lower.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=FONTSIZE_TICKS)
        ax_lower.tick_params(axis='y', labelsize=FONTSIZE_TICKS)
        ax_upper.tick_params(axis='y', labelsize=FONTSIZE_TICKS)
        ax_upper.text(0.02, 0.95, f'N = {total_sequences}', transform=ax_upper.transAxes,
                      fontsize=FONTSIZE_LABELS, verticalalignment='top',
                      bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.8))
    else:
        fig, ax = plt.subplots(figsize=FIGSIZE_MUTATION_FREQUENCY)
        bars = ax.bar(x_positions, all_proportions, width=BAR_WIDTH, color=colors_list,
                      alpha=ALPHA_BARS, edgecolor=EDGECOLOR_BARS, linewidth=LINEWIDTH_BARS)
        if show_numbers:
            for bar, proportion, count in zip(bars, all_proportions, all_counts):
                ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                        f'{proportion:.3f}\n({int(count)})',
                        ha='center', va='bottom', fontsize=FONTSIZE_TICKS - 2)
        add_aa_legend(ax)
        ax.set_xlabel("Mutation (ordered N\u2192C or 5'\u21923')", fontsize=FONTSIZE_LABELS)
        ax.set_ylabel('Frequency (proportion of total sequences)', fontsize=FONTSIZE_LABELS)
        ax.set_title(full_title, fontsize=FONTSIZE_TITLE, fontweight='bold')
        ax.set_xticks(x_positions)
        ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=FONTSIZE_TICKS)
        ax.tick_params(axis='y', labelsize=FONTSIZE_TICKS)
        if YLIM_MUTATION_FREQUENCY is not None:
            ax.set_ylim(YLIM_MUTATION_FREQUENCY)
        ax.text(0.02, 0.98, f'N = {total_sequences}', transform=ax.transAxes,
                fontsize=FONTSIZE_LABELS, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.8))
        plt.tight_layout()

    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'{filename_prefix}.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_position_distribution(df, output_dir, title, filename_prefix):
    """Histogram of mutation positions."""
    if df.empty:
        print(f"No data for {title}")
        return

    fig, ax = plt.subplots(figsize=FIGSIZE_POSITION_DISTRIBUTION)
    ax.hist(df['position'].values, bins=50, color=COLOR_SINGLE, alpha=ALPHA_BARS,
            edgecolor=EDGECOLOR_BARS, linewidth=LINEWIDTH_BARS)
    ax.set_xlabel('Position', fontsize=FONTSIZE_LABELS)
    ax.set_ylabel('Frequency', fontsize=FONTSIZE_LABELS)
    ax.set_title(title, fontsize=FONTSIZE_TITLE, fontweight='bold')
    ax.tick_params(labelsize=FONTSIZE_TICKS)
    if YLIM_POSITION_DISTRIBUTION is not None:
        ax.set_ylim(YLIM_POSITION_DISTRIBUTION)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'{filename_prefix}.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_mutation_complexity(df, output_dir, show_numbers=True):
    """Distribution of number of mutations per sequence."""
    if df.empty or 'sequence_id' not in df.columns:
        print("No data for mutation complexity")
        return

    complexity_counts = Counter(df.groupby('sequence_id').size())
    complexities = sorted(complexity_counts.keys())
    counts = [complexity_counts[c] for c in complexities]

    fig, ax = plt.subplots(figsize=FIGSIZE_COMPLEXITY)
    bars = ax.bar(complexities, counts, color=COLOR_SINGLE, alpha=ALPHA_BARS,
                  edgecolor=EDGECOLOR_BARS, linewidth=LINEWIDTH_BARS)
    if show_numbers:
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                    f'{int(count)}', ha='center', va='bottom', fontsize=FONTSIZE_TICKS)
    ax.set_xlabel('Number of Mutations per Sequence', fontsize=FONTSIZE_LABELS)
    ax.set_ylabel('Frequency', fontsize=FONTSIZE_LABELS)
    ax.set_title('Mutation Complexity Distribution', fontsize=FONTSIZE_TITLE, fontweight='bold')
    ax.tick_params(labelsize=FONTSIZE_TICKS)
    if YLIM_COMPLEXITY is not None:
        ax.set_ylim(YLIM_COMPLEXITY)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'mutation_complexity.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_substitution_heatmap(df, output_dir, mutation_type='AA',
                              show_numbers=True, remove_empty_rows=True,
                              remove_empty_cols=True, total_sequences=None):
    """
    Heatmap of original vs new residue/base. Supports count / proportion /
    per_sequence / per_residue normalization with gamma > log > linear color
    scaling (see section 8 config).
    """
    if df.empty:
        print(f"No data for {mutation_type} substitution heatmap")
        return

    if mutation_type == 'AA':
        bases = list('ACDEFGHIKLMNPQRSTVWY')
        seq_length = HEATMAP_LENGTH_AA
    else:
        bases = list('ACGT')
        seq_length = HEATMAP_LENGTH_NT

    matrix = pd.DataFrame(0, index=bases, columns=bases)
    for _, row in df.iterrows():
        if row['original'] in bases and row['new'] in bases:
            matrix.loc[row['original'], row['new']] += 1

    if remove_empty_rows:
        matrix = matrix.loc[matrix.sum(axis=1) > 0]
    if remove_empty_cols:
        matrix = matrix.loc[:, matrix.sum(axis=0) > 0]

    total_subs = int(matrix.values.sum())
    if matrix.empty or total_subs == 0:
        print(f"No substitutions found for {mutation_type} heatmap")
        return

    # --- normalize cell values ---
    if HEATMAP_NORMALIZE == 'proportion':
        display_matrix = matrix / total_subs
        cbar_label = f'Proportion of {mutation_type} substitutions'
        annot_fmt = '.3f'
    elif HEATMAP_NORMALIZE == 'per_sequence':
        display_matrix = matrix / total_sequences if total_sequences else matrix.astype(float)
        cbar_label = 'Proportion of total sequences'
        annot_fmt = '.3f'
    elif HEATMAP_NORMALIZE == 'per_residue':
        total_sites = (total_sequences or 0) * seq_length
        display_matrix = matrix / total_sites if total_sites > 0 else matrix.astype(float)
        cbar_label = ('Substitutions per residue' if mutation_type == 'AA'
                      else 'Substitutions per base')
        annot_fmt = '.2e'
    else:  # 'count'
        display_matrix = matrix
        cbar_label = 'Count'
        annot_fmt = 'g'

    # --- single, consolidated color-norm (gamma > log > linear) ---
    heatmap_data, norm, prefix = build_color_norm(
        display_matrix, gamma=HEATMAP_GAMMA, log_scale=HEATMAP_LOG_SCALE,
        vmin=HEATMAP_VMIN, vmax=HEATMAP_VMAX)
    cbar_label = prefix + cbar_label

    fig, ax = plt.subplots(figsize=FIGSIZE_HEATMAP)
    heatmap_kwargs = dict(
        annot=show_numbers, fmt=annot_fmt, cmap=COLORMAP_HEATMAP,
        cbar_kws={'label': cbar_label}, linewidths=0.5, linecolor='gray',
        ax=ax, annot_kws={'fontsize': FONTSIZE_TICKS - 2},
    )
    if norm is not None:
        heatmap_kwargs['norm'] = norm
    sns.heatmap(heatmap_data, **heatmap_kwargs)

    ax.set_xlabel(f'New {mutation_type}', fontsize=FONTSIZE_LABELS)
    ax.set_ylabel(f'Original {mutation_type}', fontsize=FONTSIZE_LABELS)
    filter_note = " (Empty rows/columns removed)" if (remove_empty_rows or remove_empty_cols) else ""
    ax.set_title(f'{mutation_type} Substitution Matrix{filter_note}',
                 fontsize=FONTSIZE_TITLE, fontweight='bold')
    ax.tick_params(labelsize=FONTSIZE_TICKS)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'{mutation_type.lower()}_substitution_heatmap.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_nt_substitution_heatmap_4x4(df, output_dir):
    """
    Dedicated 4x4 NT substitution heatmap (ACGT x ACGT). Always shows the full
    grid, optionally outlines transitions/transversions, and annotates Ts/Tv.
    Uses the same gamma > log > linear color scaling as section 9 config.
    """
    if not GENERATE_NT_SUBSTITUTION_HEATMAP:
        return

    print("  - Plotting 4x4 NT substitution heatmap (Ts/Tv)...")

    parsed_df = parse_mutations(df, mutation_type='NT')
    bases = list('ACGT')
    matrix = pd.DataFrame(0, index=bases, columns=bases)
    if not parsed_df.empty:
        for _, row in parsed_df.iterrows():
            if row['original'] in bases and row['new'] in bases:
                matrix.loc[row['original'], row['new']] += 1

    total_subs = int(matrix.values.sum())
    if total_subs == 0:
        print("    No NT substitutions found, skipping 4x4 heatmap")
        return

    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    ts_count = tv_count = 0
    for orig in bases:
        for new in bases:
            if orig == new:
                continue
            val = int(matrix.loc[orig, new])
            if (orig, new) in transitions:
                ts_count += val
            else:
                tv_count += val
    tstv_ratio = (ts_count / tv_count) if tv_count > 0 else float('inf')

    total_sequences = len(df)
    if NT_HEATMAP_NORMALIZE == 'proportion':
        display_matrix = matrix / total_subs
        cbar_label = 'Proportion of NT substitutions'
        annot_fmt = '.3f'
    elif NT_HEATMAP_NORMALIZE == 'per_sequence':
        display_matrix = matrix / total_sequences if total_sequences > 0 else matrix.astype(float)
        cbar_label = 'Proportion of total sequences'
        annot_fmt = '.3f'
    elif NT_HEATMAP_NORMALIZE == 'per_base':
        total_bases = total_sequences * NT_HEATMAP_GENE_LENGTH_BP
        display_matrix = matrix / total_bases if total_bases > 0 else matrix.astype(float)
        cbar_label = 'Substitutions per base'
        annot_fmt = '.2e'
    else:  # 'count'
        display_matrix = matrix
        cbar_label = 'Count'
        annot_fmt = 'g'

    heatmap_data, norm, prefix = build_color_norm(
        display_matrix, gamma=NT_HEATMAP_GAMMA, log_scale=NT_HEATMAP_LOG_SCALE,
        vmin=NT_HEATMAP_VMIN, vmax=NT_HEATMAP_VMAX)
    cbar_label = prefix + cbar_label

    fig, ax = plt.subplots(figsize=FIGSIZE_NT_HEATMAP)
    heatmap_kwargs = dict(
        annot=SHOW_NUMBERS_ON_NT_HEATMAP, fmt=annot_fmt,
        cmap=COLORMAP_NT_HEATMAP, cbar_kws={'label': cbar_label},
        linewidths=0.5, linecolor='gray', ax=ax, square=True,
        annot_kws={'fontsize': FONTSIZE_TICKS},
    )
    if norm is not None:
        heatmap_kwargs['norm'] = norm
    sns.heatmap(heatmap_data, **heatmap_kwargs)

    if NT_HEATMAP_HIGHLIGHT_TRANSITIONS:
        for r, orig in enumerate(bases):
            for c, new in enumerate(bases):
                if (orig, new) in transitions:
                    ax.add_patch(Rectangle((c, r), 1, 1, fill=False,
                                           edgecolor=NT_HEATMAP_TRANSITION_EDGECOLOR,
                                           linewidth=NT_HEATMAP_TRANSITION_LINEWIDTH, zorder=5))
    if NT_HEATMAP_HIGHLIGHT_TRANSVERSIONS:
        for r, orig in enumerate(bases):
            for c, new in enumerate(bases):
                if orig != new and (orig, new) not in transitions:
                    ax.add_patch(Rectangle((c, r), 1, 1, fill=False,
                                           edgecolor=NT_HEATMAP_TRANSVERSION_EDGECOLOR,
                                           linewidth=NT_HEATMAP_TRANSVERSION_LINEWIDTH, zorder=5))

    ax.set_xlabel('New NT', fontsize=FONTSIZE_LABELS)
    ax.set_ylabel('Original NT', fontsize=FONTSIZE_LABELS)
    title = 'NT Substitution Matrix (4x4)'
    if NT_HEATMAP_SHOW_TSTV_RATIO:
        if tv_count > 0:
            title += f'\nTs/Tv = {tstv_ratio:.2f}  (Ts={ts_count}, Tv={tv_count})'
        else:
            title += f'\nTs={ts_count}, Tv={tv_count}'
    ax.set_title(title, fontsize=FONTSIZE_TITLE, fontweight='bold')
    ax.tick_params(labelsize=FONTSIZE_TICKS)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'nt_substitution_heatmap_4x4.{fmt}'),
                    dpi=DPI, bbox_inches='tight')
    plt.close()


def plot_gene_map(parsed_df, df_orig, output_dir, mutation_type='AA'):
    """
    Gene map of mutation density along the sequence: nonsynonymous (yellow->red)
    and synonymous (gray->black), white = zero. Frequency = proportion of
    total sequences. (Intentionally not magma_r; see section 10 note.)
    """
    if not GENERATE_GENE_MAP or mutation_type != 'AA':
        return

    print("  - Plotting gene map...")
    total_sequences = len(df_orig)
    protein_length = GENE_MAP_GENE_LENGTH_AA

    bins_aa = range(0, protein_length + GENE_MAP_BIN_SIZE, GENE_MAP_BIN_SIZE)
    nonsyn_counts, _ = np.histogram(parsed_df['position'], bins=bins_aa)
    nonsyn_freq = nonsyn_counts / total_sequences

    # Synonymous NT positions = NT subs whose codon is not hit by a nonsyn AA sub
    syn_nt_positions = []
    if 'NT_substitutions' in df_orig.columns:
        for idx, nt_mut_string in enumerate(df_orig['NT_substitutions']):
            if pd.isna(nt_mut_string) or nt_mut_string == '':
                continue
            nt_muts = [m.strip() for m in str(nt_mut_string).split(',')]
            nonsyn_nt_positions = set()
            if 'AA_substitutions_nonsynonymous' in df_orig.columns:
                aa_nonsyn_string = df_orig.iloc[idx]['AA_substitutions_nonsynonymous']
                if pd.notna(aa_nonsyn_string) and aa_nonsyn_string != '':
                    for aa_mut in [m.strip() for m in str(aa_nonsyn_string).split(',')]:
                        if len(aa_mut) >= 3:
                            try:
                                aa_pos = int(''.join(c for c in aa_mut if c.isdigit()))
                                nt_start = (aa_pos - 1) * 3 + 1
                                nonsyn_nt_positions.update([nt_start, nt_start + 1, nt_start + 2])
                            except (ValueError, IndexError):
                                continue
            for nt_mut in nt_muts:
                try:
                    position = int(''.join(c for c in nt_mut if c.isdigit()))
                    if position not in nonsyn_nt_positions:
                        syn_nt_positions.append(position)
                except (ValueError, IndexError):
                    continue

    bins_nt = range(0, GENE_MAP_GENE_LENGTH_NT + GENE_MAP_BIN_SIZE * 3, GENE_MAP_BIN_SIZE * 3)
    bin_centers_nt = [(bins_nt[i] + bins_nt[i + 1]) / 2 for i in range(len(bins_nt) - 1)]
    if syn_nt_positions:
        syn_counts, _ = np.histogram(syn_nt_positions, bins=bins_nt)
        syn_freq = syn_counts / total_sequences
    else:
        syn_freq = np.zeros(len(bin_centers_nt))

    nonsyn_cmap = LinearSegmentedColormap.from_list(
        'nonsyn_custom', [GENE_MAP_NONSYN_LOW_COLOR, GENE_MAP_NONSYN_HIGH_COLOR], N=256)
    syn_cmap = LinearSegmentedColormap.from_list(
        'syn_custom', [GENE_MAP_SYN_LOW_COLOR, GENE_MAP_SYN_HIGH_COLOR], N=256)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=GENE_MAP_FIGSIZE,
        gridspec_kw={'height_ratios': [1, 1], 'hspace': GENE_MAP_SPACING})

    ax1.add_patch(Rectangle((0, -0.15), protein_length, 0.3,
                            facecolor=GENE_MAP_NONSYN_ZERO_COLOR, edgecolor='black',
                            linewidth=0.5, zorder=0))
    colorbar_max_nonsyn = (GENE_MAP_MAX_FREQUENCY if GENE_MAP_MAX_FREQUENCY is not None
                           else (nonsyn_freq.max() if nonsyn_freq.max() > 0 else 1))
    for i, freq in enumerate(nonsyn_freq):
        if freq > 0:
            ax1.barh(0, GENE_MAP_BIN_SIZE, left=bins_aa[i], height=0.3,
                     color=nonsyn_cmap(min(freq / colorbar_max_nonsyn, 1.0)),
                     edgecolor='none', align='center', zorder=1)
    ax1.set_xlim(0, protein_length); ax1.set_ylim(-0.5, 0.5); ax1.set_yticks([])
    ax1.set_xlabel('Amino Acid Position', fontsize=FONTSIZE_LABELS)
    ax1.set_title('Nonsynonymous Mutation Density', fontsize=FONTSIZE_LABELS, fontweight='bold')
    for s in ('top', 'right', 'left'):
        ax1.spines[s].set_visible(False)
    sm1 = plt.cm.ScalarMappable(cmap=nonsyn_cmap,
                                norm=plt.Normalize(vmin=0, vmax=colorbar_max_nonsyn))
    sm1.set_array([])
    plt.colorbar(sm1, ax=ax1, orientation='vertical', pad=0.01, aspect=10
                 ).set_label('Frequency', fontsize=FONTSIZE_TICKS)

    ax2.add_patch(Rectangle((0, -0.15), GENE_MAP_GENE_LENGTH_NT, 0.3,
                            facecolor=GENE_MAP_SYN_ZERO_COLOR, edgecolor='black',
                            linewidth=0.5, zorder=0))
    colorbar_max_syn = (GENE_MAP_MAX_FREQUENCY if GENE_MAP_MAX_FREQUENCY is not None
                        else (syn_freq.max() if syn_freq.max() > 0 else 1))
    for i, freq in enumerate(syn_freq):
        if freq > 0:
            ax2.barh(0, GENE_MAP_BIN_SIZE * 3, left=bins_nt[i], height=0.3,
                     color=syn_cmap(min(freq / colorbar_max_syn, 1.0)),
                     edgecolor='none', align='center', zorder=1)
    ax2.set_xlim(0, GENE_MAP_GENE_LENGTH_NT); ax2.set_ylim(-0.5, 0.5); ax2.set_yticks([])
    ax2.set_xlabel('Nucleotide Position', fontsize=FONTSIZE_LABELS)
    ax2.set_title('Synonymous Mutation Density', fontsize=FONTSIZE_LABELS, fontweight='bold')
    for s in ('top', 'right', 'left'):
        ax2.spines[s].set_visible(False)
    sm2 = plt.cm.ScalarMappable(cmap=syn_cmap,
                                norm=plt.Normalize(vmin=0, vmax=colorbar_max_syn))
    sm2.set_array([])
    plt.colorbar(sm2, ax=ax2, orientation='vertical', pad=0.01, aspect=10
                 ).set_label('Frequency', fontsize=FONTSIZE_TICKS)

    plt.tight_layout()
    for fmt in OUTPUT_FORMATS:
        plt.savefig(os.path.join(output_dir, f'gene_map.{fmt}'), dpi=DPI, bbox_inches='tight')
    plt.close()


# ============================================================================
# REPORT
# ============================================================================

def generate_report(stats, mutation_counts, output_dir, wt_count=0,
                    wt_syn_count=0, filename='analysis_report.txt'):
    """Text report with summary statistics and the top 20 mutations."""
    report_path = os.path.join(output_dir, filename)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("MUTATION ANALYSIS REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")
        f.write("SUMMARY STATISTICS\n" + "-" * 70 + "\n")
        f.write(f"WT sequences (no mutations): {wt_count}\n")
        f.write(f"WT-synonymous sequences: {wt_syn_count}\n")
        for key, value in stats.items():
            f.write(f"{key.replace('_', ' ').title()}: {value}\n")
        f.write("\n" + "=" * 70 + "\n")
        f.write("TOP 20 MOST FREQUENT MUTATIONS\n" + "-" * 70 + "\n")
        for i, (mutation, count) in enumerate(mutation_counts.most_common(20), 1):
            proportion = count / stats['total_sequences'] if stats['total_sequences'] > 0 else 0
            f.write(f"{i:2d}. {mutation:10s} : {count:5d} ({proportion:.4f})\n")
        f.write("\n" + "=" * 70 + "\n")
    print(f"Report saved: {report_path}")
    return report_path


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Mutation analysis for directed-evolution sequencing data: "
                    "frequency, burden, substitution matrices, Ts/Tv, and gene maps.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Input: one or more .xlsx/.csv files with one row per clone and "
               "AA_substitutions_nonsynonymous / NT_substitutions columns. "
               "See the module docstring for the full format.")
    p.add_argument('input', nargs='+',
                   help="Input .xlsx or .csv file(s); multiple files are concatenated")
    p.add_argument('-o', '--outdir', default=None,
                   help="Output directory (default: alongside the first input file)")
    p.add_argument('-m', '--mutation-type', default=MUTATION_TYPE,
                   choices=['AA', 'NT', 'BOTH'],
                   help=f"Which mutations to analyze (default: {MUTATION_TYPE})")
    p.add_argument('--formats', nargs='+', default=OUTPUT_FORMATS,
                   help=f"Figure formats to write (default: {' '.join(OUTPUT_FORMATS)})")
    p.add_argument('--dpi', type=int, default=DPI,
                   help=f"Raster DPI (default: {DPI})")
    p.add_argument('--no-subfolder', dest='use_subfolder', action='store_false',
                   help="Write directly into --outdir instead of a timestamped subfolder")
    return p.parse_args()


# ============================================================================
# MAIN
# ============================================================================

def main(input_files, output_dir, csv_filename=CSV_FILENAME, use_subfolder=True):
    """Main analysis pipeline."""
    plt.rcParams['svg.fonttype'] = SVG_FONTTYPE
    plt.rcParams['pdf.fonttype'] = PDF_FONTTYPE

    if use_subfolder and OUTPUT_SUBFOLDER:
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        output_dir = os.path.join(output_dir, f'{OUTPUT_SUBFOLDER}_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "=" * 70)
    print("MUTATION ANALYSIS PIPELINE")
    print("=" * 70 + "\n")

    # Record the exact parameters used for this run.
    write_run_parameters(output_dir, input_files=input_files)

    print("Step 1: Loading data files...")
    df = load_input_files(input_files)
    total_sequences = len(df)
    print(f"\nTotal sequences loaded: {total_sequences}")

    df.to_csv(os.path.join(output_dir, csv_filename), index=False)
    print(f"Consolidated CSV saved: {os.path.join(output_dir, csv_filename)}\n")

    if MUTATION_TYPE == 'BOTH':
        mutation_types = ['AA', 'NT']
    elif MUTATION_TYPE in ('AA', 'NT'):
        mutation_types = [MUTATION_TYPE]
    else:
        raise ValueError("MUTATION_TYPE must be 'AA', 'NT', or 'BOTH'")

    all_results = {}
    for mut_type in mutation_types:
        print(f"\n{'=' * 70}\nANALYZING {mut_type} MUTATIONS\n{'=' * 70}\n")

        print("Step 2: Counting WT sequences...")
        wt_count, wt_syn_count = count_wt_sequences(df, mutation_type=mut_type)
        print(f"  WT sequences: {wt_count}")
        print(f"  WT-synonymous sequences: {wt_syn_count}")

        print(f"Step 3: Parsing {mut_type} mutations...")
        parsed_df = parse_mutations(df, mutation_type=mut_type)
        if parsed_df.empty and wt_count == 0 and wt_syn_count == 0:
            print(f"No {mut_type} mutations or WT sequences found, skipping...\n")
            continue

        print("Step 4: Calculating statistics...")
        stats = calculate_statistics(parsed_df, total_sequences)
        mutation_counts = Counter(parsed_df['mutation']) if not parsed_df.empty else Counter()

        print("Step 5: Generating visualizations...\n")
        print(f"  - Plotting {mut_type} mutation frequencies...")
        plot_mutation_frequency(
            mutation_counts, output_dir,
            f'{mut_type} Mutation Frequency Distribution',
            f'{mut_type.lower()}_mutation_frequency',
            color_by_aa=(COLOR_BY_AMINO_ACID and mut_type == 'AA'),
            display_n=DISPLAY_N_MUTATIONS, top_abundant=DISPLAY_TOP_ABUNDANT,
            show_numbers=SHOW_NUMBERS_ON_BARS, total_sequences=total_sequences,
            wt_count=wt_count, wt_syn_count=wt_syn_count)

        if not parsed_df.empty:
            print(f"  - Plotting {mut_type} position distribution...")
            plot_position_distribution(
                parsed_df, output_dir,
                f'{mut_type} Mutation Position Distribution',
                f'{mut_type.lower()}_position_distribution')

            print(f"  - Plotting {mut_type} mutation complexity...")
            plot_mutation_complexity(parsed_df, output_dir,
                                     show_numbers=SHOW_NUMBERS_ON_COMPLEXITY)

            print(f"  - Plotting {mut_type} substitution heatmap...")
            plot_substitution_heatmap(
                parsed_df, output_dir, mutation_type=mut_type,
                show_numbers=SHOW_NUMBERS_ON_HEATMAP,
                remove_empty_rows=REMOVE_EMPTY_HEATMAP_ROWS,
                remove_empty_cols=REMOVE_EMPTY_HEATMAP_COLS,
                total_sequences=total_sequences)

            plot_gene_map(parsed_df, df, output_dir, mutation_type=mut_type)

        print(f"\nStep 6: Generating {mut_type} report...")
        generate_report(stats, mutation_counts, output_dir,
                        wt_count=wt_count, wt_syn_count=wt_syn_count,
                        filename=f'{mut_type.lower()}_analysis_report.txt')

        all_results[mut_type] = {
            'parsed_df': parsed_df, 'mutation_counts': mutation_counts,
            'stats': stats, 'wt_count': wt_count, 'wt_syn_count': wt_syn_count,
        }

    if GENERATE_BURDEN_HISTOGRAMS:
        print(f"\n{'=' * 70}\nGENERATING MUTATION BURDEN HISTOGRAMS (AA + NT)\n{'=' * 70}\n")
        plot_mutation_burden(df, output_dir, mutation_type='AA')
        plot_mutation_burden(df, output_dir, mutation_type='NT')

    if GENERATE_NT_SUBSTITUTION_HEATMAP:
        print(f"\n{'=' * 70}\nGENERATING 4x4 NT SUBSTITUTION HEATMAP (Ts/Tv)\n{'=' * 70}\n")
        plot_nt_substitution_heatmap_4x4(df, output_dir)

    print(f"\n{'=' * 70}\nANALYSIS COMPLETE!\n{'=' * 70}")
    for mut_type in mutation_types:
        if mut_type in all_results:
            r = all_results[mut_type]
            print(f"\n{mut_type} Mutations:")
            print(f"  WT sequences: {r['wt_count']}")
            print(f"  WT-synonymous sequences: {r['wt_syn_count']}")
            print(f"  Total mutations: {r['stats']['total_mutations']}")
            print(f"  Unique mutations: {r['stats']['unique_mutations']}")
            print(f"  Unique positions: {r['stats']['unique_positions']}")
    print(f"\nAll files saved to: {output_dir}\n{'=' * 70}\n")
    return all_results


if __name__ == '__main__':
    args = parse_args()

    # CLI overrides for the module-level config, so write_run_parameters()
    # records the values actually used for this run.
    MUTATION_TYPE = args.mutation_type
    OUTPUT_FORMATS = list(args.formats)
    DPI = args.dpi

    input_paths = [os.path.abspath(os.path.expanduser(p)) for p in args.input]
    out_dir = (os.path.abspath(os.path.expanduser(args.outdir)) if args.outdir
               else os.path.dirname(input_paths[0]) or os.getcwd())

    results = main(input_paths, out_dir, use_subfolder=args.use_subfolder)
