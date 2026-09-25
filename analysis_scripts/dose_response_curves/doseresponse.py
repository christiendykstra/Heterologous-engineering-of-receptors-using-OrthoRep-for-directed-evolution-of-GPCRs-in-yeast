#!/usr/bin/env python3

"""
dose_response_fit.py — fit agonist dose-response curves and export figures + EC50 table.

Copyright (c) 2026 Christien Bret Dykstra
SPDX-License-Identifier: MIT

Part of: https://github.com/christiendykstra/Heterologous-engineering-of-receptors-using-OrthoRep-for-directed-evolution-of-GPCRs-in-yeast

Fits a four-parameter logistic (4PL/Hill) model to concentration-response data and
writes SVG/PDF/PNG figures plus a CSV of fitted parameters (EC50, Hill slope,
asymptotes, R^2, confidence intervals).

INPUT FORMAT
------------
A single-sheet .xlsx file laid out as:

    column 0        : log10[ligand] in molar (e.g. -9, -8.5, -8, ... -4)
    columns 1..N    : response values, grouped by sample in blocks of
                      REPLICATES_PER_SAMPLE consecutive columns

Sample names are taken from the first replicate column of each block, with
common replicate suffixes stripped (_Rep1, _R1, _1, .1, etc.). Blank cells are
ignored. A sample needs at least 4 usable concentration points to be fitted.

Example header for 2 samples x 4 replicates:

    logM | WT_Rep1 | WT_Rep2 | WT_Rep3 | WT_Rep4 | A51T_Rep1 | A51T_Rep2 | ...

USAGE
-----
    python dose_response_fit.py data.xlsx
    python dose_response_fit.py data.xlsx --outdir results/ --replicates 3
    python dose_response_fit.py data.xlsx --model 4PL 5PL --xlabel "Log10[DAMGO (M)]"

Run with --help for the full option list. Defaults for every plotting and model
parameter live in the CONFIG block below and can be edited directly.
"""

import argparse
import os
import re
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

# =============================================================================
# CONFIG — edit defaults here, or override the starred ones from the command line
# =============================================================================

# ----- Output -----
OUTPUT_DIR = None            # * None = write next to the input file
OUTPUT_FORMATS = ['svg', 'pdf', 'png']
DPI = 300

# ----- Data layout -----
REPLICATES_PER_SAMPLE = 4    # * consecutive columns per sample
MIN_POINTS_TO_FIT = 4        # samples with fewer usable points are skipped

# ----- Model -----
MODELS_TO_FIT = ['4PL']      # * '4PL', '3PL_bottom', '3PL_top', '5PL', 'biphasic'
BOTTOM_CONSTRAINT = None     # e.g. 0.0 to fix the baseline
TOP_CONSTRAINT = None        # e.g. 100.0 to fix the maximum
CI_LEVEL = 0.95

# ----- Axis labels -----
XLABEL = "Log10[ligand (M)]"   # Can change to different gradient label
YLABEL = "% Cell signalling"   # Can change to be FOLD-OD

# ----- Fonts -----
mpl.rcParams['svg.fonttype'] = 'none'   # keep SVG text editable in Illustrator
FONT_FAMILY = 'sans-serif'
FONT_NAME = 'Arial'
TITLE_FONTSIZE = 6
AXIS_LABEL_FONTSIZE = 6
LEGEND_FONTSIZE = 4
TICK_LABELSIZE = 4

# ----- Figure size (inches) -----
FIGURE_WIDTH = 3
FIGURE_HEIGHT = 3

# ----- Colour palette (cycled across samples in column order) -----
COLORS = [
    '#000000',  # black
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#17becf',  # cyan
    '#bcbd22',  # olive
    '#fbb040',  # tangerine
    '#9161a9',  # lavender
    '#2ca02c',  # green
    '#ff7f0e',  # orange
]

# ----- Points and error bars -----
SHOW_ERROR_BARS = True       # True = mean +/- error; False = individual replicates
ERROR_BAR_TYPE = 'SD'        # 'SD' or 'SEM'
ERROR_BAR_CAPSIZE = 1.5
ERROR_BAR_LINEWIDTH = 0.25
ERROR_BAR_COLOR = 'black'    # None = match sample colour
REPLICATE_POINT_SIZE = 2
REPLICATE_ALPHA = 0.4
AVERAGE_POINT_SIZE = 4
AVERAGE_MARKER = 'o'
MARKER_FILL = 'white'        # 'full' | 'white' | 'none'

# ----- Lines and axes -----
CURVE_LINE_WIDTH = 0.5
LINE_STYLE = '-'
AXIS_LINE_WIDTH = 0.5
TICK_WIDTH = 0.5
TICK_LENGTH = 2

# ----- Limits, ticks, grid -----
Y_MAX_LIMIT = 120            # * None = auto from data
X_MAJOR_TICKS = None         # None = auto, or e.g. [-9, -8, -7, -6, -5, -4]
Y_MAJOR_TICKS = [0, 20, 40, 60, 80, 100]
SHOW_GRID = False
GRID_ALPHA = 0.3
GRID_WHICH = 'major'

# ----- Zero-concentration break -----
SHOW_ZERO_BREAK = True       # relabel ZERO_TICK_VALUE as "0" and draw break marks
ZERO_TICK_VALUE = -9
BREAK_POSITION = -8.5

# ----- Legend -----
LEGEND_POSITION = (0, 1.14)
LEGEND_LOC = 'lower left'

# =============================================================================
# End of CONFIG
# =============================================================================


def apply_rcparams():
    mpl.rcParams['font.family'] = FONT_FAMILY
    mpl.rcParams['font.sans-serif'] = [FONT_NAME]
    mpl.rcParams['axes.labelsize'] = AXIS_LABEL_FONTSIZE
    mpl.rcParams['xtick.labelsize'] = TICK_LABELSIZE
    mpl.rcParams['ytick.labelsize'] = TICK_LABELSIZE
    mpl.rcParams['axes.linewidth'] = AXIS_LINE_WIDTH
    mpl.rcParams['xtick.major.width'] = TICK_WIDTH
    mpl.rcParams['ytick.major.width'] = TICK_WIDTH
    mpl.rcParams['xtick.major.size'] = TICK_LENGTH
    mpl.rcParams['ytick.major.size'] = TICK_LENGTH


# ===== MODEL DEFINITIONS =====

def sigmoid_4pl(x, bottom, top, ec50, hill):
    """4-Parameter Logistic (Hill equation)"""
    return bottom + (top - bottom) / (1 + 10 ** ((ec50 - x) * hill))


def sigmoid_5pl(x, bottom, top, ec50, hill, asymmetry):
    """5-Parameter Logistic (asymmetric sigmoid)"""
    return bottom + (top - bottom) / (1 + 10 ** ((ec50 - x) * hill)) ** asymmetry


def biphasic(x, bottom, top1, ec50_1, hill1, top2, ec50_2, hill2):
    """Biphasic dose-response (sum of two sigmoids)"""
    response1 = (top1 - bottom) / (1 + 10 ** ((ec50_1 - x) * hill1))
    response2 = (top2 - top1) / (1 + 10 ** ((ec50_2 - x) * hill2))
    return bottom + response1 + response2


def calculate_r_squared(y_actual, y_predicted):
    """Calculate R^2 (coefficient of determination)"""
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1 - (ss_res / ss_tot)


def calculate_confidence_intervals(popt, pcov, ci_level=0.95):
    """
    Calculate confidence intervals for fitted parameters.

    NOTE: degrees of freedom are fixed at 100 rather than (n - p). With a typical
    8-12 point curve this makes the reported intervals slightly narrower than a
    strict n-p treatment would give. Point estimates (EC50, Hill) are unaffected.
    Kept as-is so output matches the figures in the manuscript.
    """
    alpha = 1 - ci_level
    t_val = stats.t.ppf(1 - alpha / 2, df=100)

    ci_lower = []
    ci_upper = []

    for i in range(len(popt)):
        if pcov is not None and np.isfinite(pcov[i, i]) and pcov[i, i] >= 0:
            std_err = np.sqrt(pcov[i, i])
            ci_lower.append(popt[i] - t_val * std_err)
            ci_upper.append(popt[i] + t_val * std_err)
        else:
            ci_lower.append(np.nan)
            ci_upper.append(np.nan)

    return ci_lower, ci_upper


def extract_sample_name(columns):
    """Extract sample name from the first replicate column header"""
    first_col = str(columns[0])
    name = first_col
    for suffix in ['_Rep1', '_rep1', '_R1', '_r1', '-Rep1', '-rep1', '-R1', '-r1',
                   '_1', '-1', ' Rep1', ' rep1', ' R1', ' r1', ' 1']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    name = re.sub(r'\.\d+$', '', name)   # pandas duplicate-column suffixes
    name = name.rstrip('_-. ')
    return name if name else first_col


def fit_model(model_type, x_data, y_data, bottom_constraint=None, top_constraint=None):
    """Fit specified model to data and return results"""
    y_min, y_max = min(y_data), max(y_data)
    x_mid = np.median(x_data)
    increasing = y_data[-1] > y_data[0]

    results = {
        'model': model_type,
        'success': False,
        'popt': None,
        'pcov': None,
        'r_squared': np.nan,
        'ec50': np.nan,
        'ec50_ci_lower': np.nan,
        'ec50_ci_upper': np.nan,
        'hill': np.nan,
        'hill_ci_lower': np.nan,
        'hill_ci_upper': np.nan,
        'bottom': np.nan,
        'top': np.nan,
    }

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if model_type == '4PL':
                if increasing:
                    p0 = [y_min, y_max, x_mid, 1.0]
                else:
                    p0 = [y_max, y_min, x_mid, 1.0]

                bounds_lo = [-np.inf, -np.inf, -np.inf, 0.1]
                bounds_hi = [np.inf, np.inf, np.inf, 10]

                if bottom_constraint is not None:
                    p0[0] = bottom_constraint
                    bounds_lo[0] = bottom_constraint - 0.001
                    bounds_hi[0] = bottom_constraint + 0.001
                if top_constraint is not None:
                    p0[1] = top_constraint
                    bounds_lo[1] = top_constraint - 0.001
                    bounds_hi[1] = top_constraint + 0.001

                popt, pcov = curve_fit(sigmoid_4pl, x_data, y_data, p0=p0,
                                       bounds=(bounds_lo, bounds_hi),
                                       method='trf', maxfev=10000)

                results['popt'] = popt
                results['pcov'] = pcov
                results['bottom'] = popt[0]
                results['top'] = popt[1]
                results['ec50'] = popt[2]
                results['hill'] = popt[3]

                y_pred = sigmoid_4pl(x_data, *popt)

                ci_lower, ci_upper = calculate_confidence_intervals(popt, pcov, CI_LEVEL)
                results['ec50_ci_lower'] = ci_lower[2]
                results['ec50_ci_upper'] = ci_upper[2]
                results['hill_ci_lower'] = ci_lower[3]
                results['hill_ci_upper'] = ci_upper[3]

            elif model_type == '3PL_bottom':
                bottom_fixed = bottom_constraint if bottom_constraint is not None else y_min

                def model_func(x, top, ec50, hill):
                    return bottom_fixed + (top - bottom_fixed) / (1 + 10 ** ((ec50 - x) * hill))

                p0 = [y_max, x_mid, 1.0] if increasing else [y_min, x_mid, 1.0]

                popt, pcov = curve_fit(model_func, x_data, y_data, p0=p0,
                                       bounds=([-np.inf, -np.inf, 0.1], [np.inf, np.inf, 10]),
                                       method='trf', maxfev=10000)

                results['popt'] = np.array([bottom_fixed, popt[0], popt[1], popt[2]])
                results['pcov'] = pcov
                results['bottom'] = bottom_fixed
                results['top'] = popt[0]
                results['ec50'] = popt[1]
                results['hill'] = popt[2]

                y_pred = model_func(x_data, *popt)

                ci_lower, ci_upper = calculate_confidence_intervals(popt, pcov, CI_LEVEL)
                results['ec50_ci_lower'] = ci_lower[1]
                results['ec50_ci_upper'] = ci_upper[1]
                results['hill_ci_lower'] = ci_lower[2]
                results['hill_ci_upper'] = ci_upper[2]

            elif model_type == '3PL_top':
                top_fixed = top_constraint if top_constraint is not None else y_max

                def model_func(x, bottom, ec50, hill):
                    return bottom + (top_fixed - bottom) / (1 + 10 ** ((ec50 - x) * hill))

                p0 = [y_min, x_mid, 1.0] if increasing else [y_max, x_mid, 1.0]

                popt, pcov = curve_fit(model_func, x_data, y_data, p0=p0,
                                       bounds=([-np.inf, -np.inf, 0.1], [np.inf, np.inf, 10]),
                                       method='trf', maxfev=10000)

                results['popt'] = np.array([popt[0], top_fixed, popt[1], popt[2]])
                results['pcov'] = pcov
                results['bottom'] = popt[0]
                results['top'] = top_fixed
                results['ec50'] = popt[1]
                results['hill'] = popt[2]

                y_pred = model_func(x_data, *popt)

                ci_lower, ci_upper = calculate_confidence_intervals(popt, pcov, CI_LEVEL)
                results['ec50_ci_lower'] = ci_lower[1]
                results['ec50_ci_upper'] = ci_upper[1]
                results['hill_ci_lower'] = ci_lower[2]
                results['hill_ci_upper'] = ci_upper[2]

            elif model_type == '5PL':
                if increasing:
                    p0 = [y_min, y_max, x_mid, 1.0, 1.0]
                else:
                    p0 = [y_max, y_min, x_mid, 1.0, 1.0]

                popt, pcov = curve_fit(sigmoid_5pl, x_data, y_data, p0=p0,
                                       bounds=([-np.inf, -np.inf, -np.inf, 0.1, 0.1],
                                               [np.inf, np.inf, np.inf, 10, 10]),
                                       method='trf', maxfev=10000)

                results['popt'] = popt
                results['pcov'] = pcov
                results['bottom'] = popt[0]
                results['top'] = popt[1]
                results['ec50'] = popt[2]
                results['hill'] = popt[3]
                results['asymmetry'] = popt[4]

                y_pred = sigmoid_5pl(x_data, *popt)

                ci_lower, ci_upper = calculate_confidence_intervals(popt, pcov, CI_LEVEL)
                results['ec50_ci_lower'] = ci_lower[2]
                results['ec50_ci_upper'] = ci_upper[2]
                results['hill_ci_lower'] = ci_lower[3]
                results['hill_ci_upper'] = ci_upper[3]

            elif model_type == 'biphasic':
                x_range = max(x_data) - min(x_data)
                mid_y = (y_min + y_max) / 2

                if increasing:
                    p0 = [y_min, mid_y, x_mid - x_range / 4, 1.0, y_max, x_mid + x_range / 4, 1.0]
                else:
                    p0 = [y_max, mid_y, x_mid - x_range / 4, 1.0, y_min, x_mid + x_range / 4, 1.0]

                popt, pcov = curve_fit(biphasic, x_data, y_data, p0=p0,
                                       method='trf', maxfev=20000)

                results['popt'] = popt
                results['pcov'] = pcov
                results['bottom'] = popt[0]
                results['top'] = popt[4]
                results['ec50'] = popt[2]
                results['ec50_2'] = popt[5]
                results['hill'] = popt[3]
                results['hill_2'] = popt[6]

                y_pred = biphasic(x_data, *popt)

                ci_lower, ci_upper = calculate_confidence_intervals(popt, pcov, CI_LEVEL)
                results['ec50_ci_lower'] = ci_lower[2]
                results['ec50_ci_upper'] = ci_upper[2]
                results['hill_ci_lower'] = ci_lower[3]
                results['hill_ci_upper'] = ci_upper[3]

            results['r_squared'] = calculate_r_squared(y_data, y_pred)
            results['success'] = True

    except Exception as e:
        results['error'] = str(e)

    return results


def generate_fit_curve(model_type, popt, x_range, bottom_constraint=None):
    """Generate smooth curve for plotting"""
    x_fit = np.linspace(min(x_range), max(x_range), 200)

    if model_type == '4PL':
        y_fit = sigmoid_4pl(x_fit, *popt)
    elif model_type == '3PL_bottom':
        bottom_fixed = bottom_constraint if bottom_constraint is not None else popt[0]
        y_fit = bottom_fixed + (popt[1] - bottom_fixed) / (1 + 10 ** ((popt[2] - x_fit) * popt[3]))
    elif model_type == '3PL_top':
        y_fit = popt[0] + (popt[1] - popt[0]) / (1 + 10 ** ((popt[2] - x_fit) * popt[3]))
    elif model_type == '5PL':
        y_fit = sigmoid_5pl(x_fit, *popt)
    elif model_type == 'biphasic':
        y_fit = biphasic(x_fit, *popt)
    else:
        y_fit = np.zeros_like(x_fit)

    return x_fit, y_fit


def setup_axis_ticks(ax, x_data, y_data):
    """Configure axis ticks"""
    if X_MAJOR_TICKS is not None:
        ax.set_xticks(X_MAJOR_TICKS)
    else:
        x_min, x_max = min(x_data), max(x_data)
        x_ticks = np.arange(np.floor(x_min), np.ceil(x_max) + 1, 1)
        ax.set_xticks(x_ticks)

        if SHOW_ZERO_BREAK:
            labels = ['0' if t == ZERO_TICK_VALUE else str(int(t)) for t in x_ticks]
            ax.set_xticklabels(labels)

    if Y_MAJOR_TICKS is not None:
        ax.set_yticks(Y_MAJOR_TICKS)
    else:
        y_min, y_max = min(y_data), max(y_data)
        y_range = y_max - y_min

        if y_range <= 5:
            interval = 1
        elif y_range <= 10:
            interval = 2
        elif y_range <= 25:
            interval = 5
        elif y_range <= 50:
            interval = 10
        else:
            interval = 20

        y_start = np.floor(y_min / interval) * interval
        y_end = np.ceil(y_max / interval) * interval
        ax.set_yticks(np.arange(y_start, y_end + interval, interval))


def draw_axis_break(ax, x_pos):
    """Draw diagonal break lines directly on the axis spines (frame)"""
    xlim = ax.get_xlim()
    x_frac = (x_pos - xlim[0]) / (xlim[1] - xlim[0])
    d = 0.015

    ax.plot([x_frac - d, x_frac + d], [-d, d],
            transform=ax.transAxes, color='k', linewidth=0.8, clip_on=False)
    ax.plot([x_frac - d, x_frac + d], [1 - d, 1 + d],
            transform=ax.transAxes, color='k', linewidth=0.8, clip_on=False)


def calculate_error_bars(df, rep_columns):
    """Calculate SEM or SD for each concentration across replicates"""
    errors = []
    for _, row in df.iterrows():
        replicate_values = []
        for rep_col in rep_columns:
            y_val = row[rep_col]
            if pd.notna(y_val):
                try:
                    y_val = float(y_val)
                    if not np.isnan(y_val) and not np.isinf(y_val):
                        replicate_values.append(y_val)
                except (ValueError, TypeError):
                    pass

        if replicate_values:
            arr = np.array(replicate_values)
            if ERROR_BAR_TYPE.upper() == 'SEM':
                errors.append(stats.sem(arr) if len(arr) > 1 else 0.0)
            else:
                errors.append(np.std(arr, ddof=1) if len(arr) > 1 else 0.0)
        else:
            errors.append(0.0)

    return np.array(errors)


def create_plot_for_model(model_type, sample_data_list, input_filename, output_folder,
                          xlabel, ylabel, y_max_limit):
    """Create a complete plot for one model type"""
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

    all_x, all_y = [], []
    successful_fits = 0

    for sample_idx, sample_data in enumerate(sample_data_list):
        sample_name = sample_data['name']
        all_x_points = sample_data['all_x']
        all_y_points = sample_data['all_y']
        avg_x = sample_data['avg_x']
        avg_y = sample_data['avg_y']

        all_x.extend(all_x_points)
        all_y.extend(all_y_points)

        color = COLORS[sample_idx % len(COLORS)]

        if SHOW_ERROR_BARS:
            error_values = sample_data['errors']
            _mfc = color if MARKER_FILL == 'full' else MARKER_FILL
            _ecolor = ERROR_BAR_COLOR if ERROR_BAR_COLOR is not None else color
            ax.errorbar(avg_x, avg_y, yerr=error_values,
                        fmt=AVERAGE_MARKER, markersize=np.sqrt(AVERAGE_POINT_SIZE),
                        color=color, ecolor=_ecolor,
                        markerfacecolor=_mfc, markeredgecolor=color,
                        capsize=ERROR_BAR_CAPSIZE, elinewidth=ERROR_BAR_LINEWIDTH,
                        capthick=ERROR_BAR_LINEWIDTH, zorder=3)
        else:
            ax.scatter(all_x_points, all_y_points, s=REPLICATE_POINT_SIZE,
                       alpha=REPLICATE_ALPHA, color=color)
            _fc = color if MARKER_FILL == 'full' else MARKER_FILL
            ax.scatter(avg_x, avg_y, s=AVERAGE_POINT_SIZE, alpha=1.0,
                       facecolors=_fc, edgecolors=color, marker=AVERAGE_MARKER)

        # NOTE: the fit is performed on per-concentration means (avg_y), not on the
        # individual replicates. This is what produced the published figures.
        fit_result = fit_model(model_type, avg_x, avg_y, BOTTOM_CONSTRAINT, TOP_CONSTRAINT)

        if fit_result['success']:
            x_fit, y_fit = generate_fit_curve(model_type, fit_result['popt'], avg_x,
                                              BOTTOM_CONSTRAINT)
            ec50_um = 10 ** (fit_result['ec50']) * 1e6

            label = f"{sample_name} (EC50: {ec50_um:.2f} \u00b5M, nH={fit_result['hill']:.2f})"
            ax.plot(x_fit, y_fit, LINE_STYLE, linewidth=CURVE_LINE_WIDTH,
                    color=color, label=label)

            successful_fits += 1
            sample_data[f'{model_type}_fit'] = fit_result

    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(f"{input_filename} - {model_type} ({successful_fits}/{len(sample_data_list)} fits)",
                 fontsize=TITLE_FONTSIZE)

    if SHOW_GRID:
        ax.grid(True, which=GRID_WHICH, alpha=GRID_ALPHA)

    setup_axis_ticks(ax, all_x, all_y)

    ax.set_xlim(min(all_x) - 0.3, max(all_x) + 0.3)

    if y_max_limit is not None:
        ax.set_ylim(top=y_max_limit)

    if SHOW_ZERO_BREAK:
        draw_axis_break(ax, BREAK_POSITION)

    ax.legend(bbox_to_anchor=LEGEND_POSITION, loc=LEGEND_LOC, fontsize=LEGEND_FONTSIZE)
    ax.set_position([0.1, 0.1, 0.6, 0.6])   # square plot, legend outside

    for fmt in OUTPUT_FORMATS:
        output_path = os.path.join(output_folder, f'{input_filename}_{model_type}.{fmt}')
        fig.savefig(output_path, format=fmt, dpi=DPI if fmt == 'png' else None,
                    bbox_inches=None, pad_inches=0.1)
        print(f"Saved: {output_path}")

    plt.close(fig)
    return successful_fits


def load_samples(df, conc_column_name, replicates_per_sample):
    """Split the data columns into per-sample blocks and collect points"""
    all_columns = df.columns.tolist()[1:]
    num_samples = len(all_columns) // replicates_per_sample

    leftover = len(all_columns) % replicates_per_sample
    if leftover:
        print(f"Warning: {leftover} trailing column(s) ignored "
              f"({len(all_columns)} data columns is not a multiple of "
              f"{replicates_per_sample})")

    print(f"Detected {num_samples} samples with {replicates_per_sample} replicates each")

    sample_data_list = []

    for sample_idx in range(num_samples):
        start_col = sample_idx * replicates_per_sample
        rep_columns = all_columns[start_col:start_col + replicates_per_sample]
        sample_name = extract_sample_name(rep_columns)
        print(f"\nSample {sample_idx + 1}: '{sample_name}' from columns {rep_columns}")

        all_x_points, all_y_points = [], []
        avg_x_values, avg_y_values = [], []

        for _, row in df.iterrows():
            x_val = row[conc_column_name]
            replicate_values = []

            for rep_col in rep_columns:
                y_val = row[rep_col]
                if pd.notna(y_val):
                    try:
                        y_val = float(y_val)
                        if not np.isnan(y_val) and not np.isinf(y_val):
                            replicate_values.append(y_val)
                            all_x_points.append(x_val)
                            all_y_points.append(y_val)
                    except (ValueError, TypeError):
                        pass

            if replicate_values:
                avg_x_values.append(x_val)
                avg_y_values.append(np.mean(replicate_values))

        errors = calculate_error_bars(df, rep_columns) if SHOW_ERROR_BARS else None

        if len(avg_x_values) >= MIN_POINTS_TO_FIT:
            sample_data_list.append({
                'name': sample_name,
                'columns': rep_columns,
                'all_x': all_x_points,
                'all_y': all_y_points,
                'avg_x': np.array(avg_x_values),
                'avg_y': np.array(avg_y_values),
                'errors': errors,
            })
        else:
            print(f"  Skipping: insufficient data points ({len(avg_x_values)})")

    return sample_data_list


def build_results_table(sample_data_list, models):
    """Assemble the per-sample fitted-parameter rows"""
    all_results = []

    for sample_data in sample_data_list:
        result_row = {
            'Sample': sample_data['name'],
            'N_Points': len(sample_data['avg_x']),
        }

        for model_type in models:
            fit_key = f'{model_type}_fit'
            if fit_key in sample_data and sample_data[fit_key]['success']:
                fit = sample_data[fit_key]
                ec50_log = fit['ec50']

                result_row[f'{model_type}_EC50_log10M'] = ec50_log
                result_row[f'{model_type}_EC50_uM'] = 10 ** ec50_log * 1e6

                if not np.isnan(fit['ec50_ci_lower']):
                    result_row[f'{model_type}_EC50_CI_lower_uM'] = 10 ** fit['ec50_ci_lower'] * 1e6
                    result_row[f'{model_type}_EC50_CI_upper_uM'] = 10 ** fit['ec50_ci_upper'] * 1e6

                result_row[f'{model_type}_Hill'] = fit['hill']
                result_row[f'{model_type}_Hill_CI_lower'] = fit.get('hill_ci_lower', np.nan)
                result_row[f'{model_type}_Hill_CI_upper'] = fit.get('hill_ci_upper', np.nan)
                result_row[f'{model_type}_Bottom'] = fit['bottom']
                result_row[f'{model_type}_Top'] = fit['top']
                result_row[f'{model_type}_R2'] = fit['r_squared']

                if model_type == '5PL' and 'asymmetry' in fit:
                    result_row[f'{model_type}_Asymmetry'] = fit['asymmetry']
                if model_type == 'biphasic':
                    result_row[f'{model_type}_EC50_2_log10M'] = fit.get('ec50_2', np.nan)
                    result_row[f'{model_type}_Hill_2'] = fit.get('hill_2', np.nan)

        all_results.append(result_row)

    return all_results


def print_summary(all_results, models):
    print(f"\n{'=' * 80}")
    print("FIT RESULTS SUMMARY")
    print(f"{'=' * 80}")

    for result in all_results:
        print(f"\n{result['Sample']}:")
        for model in models:
            ec50_key = f'{model}_EC50_uM'
            if ec50_key in result and not np.isnan(result.get(ec50_key, np.nan)):
                ec50 = result[ec50_key]
                hill = result.get(f'{model}_Hill', np.nan)
                r2 = result.get(f'{model}_R2', np.nan)
                ci_lo = result.get(f'{model}_EC50_CI_lower_uM', np.nan)
                ci_hi = result.get(f'{model}_EC50_CI_upper_uM', np.nan)

                ci_str = f"[{ci_lo:.2f} - {ci_hi:.2f}]" if not np.isnan(ci_lo) else "[CI unavailable]"
                print(f"  {model}: EC50 = {ec50:.3f} \u00b5M {ci_str}, "
                      f"nH = {hill:.3f}, R\u00b2 = {r2:.4f}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit dose-response curves (4PL and related models) and export "
                    "figures plus an EC50 table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Input: single-sheet .xlsx, column 0 = log10[ligand] in M, then "
               "N samples x --replicates response columns.",
    )
    parser.add_argument('input', help="Path to the input .xlsx file")
    parser.add_argument('-o', '--outdir', default=OUTPUT_DIR,
                        help="Output directory (default: alongside the input file)")
    parser.add_argument('-r', '--replicates', type=int, default=REPLICATES_PER_SAMPLE,
                        help=f"Replicate columns per sample (default: {REPLICATES_PER_SAMPLE})")
    parser.add_argument('-m', '--model', nargs='+', default=MODELS_TO_FIT,
                        choices=['4PL', '3PL_bottom', '3PL_top', '5PL', 'biphasic'],
                        help=f"Model(s) to fit (default: {' '.join(MODELS_TO_FIT)})")
    parser.add_argument('--sheet', default=0,
                        help="Sheet name or index to read (default: first sheet)")
    parser.add_argument('--xlabel', default=XLABEL, help="X axis label")
    parser.add_argument('--ylabel', default=YLABEL, help="Y axis label")
    parser.add_argument('--ymax', type=float, default=Y_MAX_LIMIT,
                        help="Y axis upper limit (default from CONFIG; use -1 for auto)")
    return parser.parse_args()


def main():
    args = parse_args()
    apply_rcparams()

    input_path = os.path.abspath(os.path.expanduser(args.input))
    if not os.path.isfile(input_path):
        raise SystemExit(f"Input file not found: {input_path}")

    sheet = args.sheet
    if isinstance(sheet, str) and sheet.isdigit():
        sheet = int(sheet)

    df = pd.read_excel(input_path, sheet_name=sheet)
    input_filename = os.path.splitext(os.path.basename(input_path))[0]

    output_folder = args.outdir or os.path.dirname(input_path) or os.getcwd()
    os.makedirs(output_folder, exist_ok=True)

    y_max_limit = None if (args.ymax is not None and args.ymax < 0) else args.ymax

    print(f"Excel file loaded: {input_path}")
    print(f"Shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")

    conc_column_name = df.columns[0]
    print(f"\nConcentration column: '{conc_column_name}'")

    sample_data_list = load_samples(df, conc_column_name, args.replicates)

    if not sample_data_list:
        raise SystemExit("No samples had enough usable data points to fit.")

    print(f"\n{'=' * 60}")
    print(f"Processing {len(sample_data_list)} samples with {len(args.model)} model(s)")
    print(f"{'=' * 60}")

    for model_type in args.model:
        print(f"\n--- Generating {model_type} plot ---")
        n_fits = create_plot_for_model(model_type, sample_data_list, input_filename,
                                       output_folder, args.xlabel, args.ylabel, y_max_limit)
        print(f"Successful fits: {n_fits}/{len(sample_data_list)}")

    all_results = build_results_table(sample_data_list, args.model)

    results_df = pd.DataFrame(all_results)
    csv_path = os.path.join(output_folder, f'{input_filename}_fit_results.csv')
    results_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    print_summary(all_results, args.model)
    print("\nDone!")


if __name__ == '__main__':
    main()