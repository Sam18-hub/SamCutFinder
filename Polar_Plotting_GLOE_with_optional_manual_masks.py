import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# USER SETTINGS
# =============================================================================

# Use the plot-ready CSV produced directly by the corrected GLOE-seq script.
INPUT_FILE = "raw_nick_position_plot_ready.csv"

# Output name.
OUTPUT_PNG = "gloe_polar_plot.png"
OUTPUT_PDF = "gloe_polar_plot.pdf"

# False = raw reads; True = BPM.
# For a single sample, the spatial shape is identical because BPM is just a
# constant rescaling. For comparisons between libraries, BPM is preferable,
# but without an internal standard it is NOT absolute cleavage efficiency.
USE_BPM = True

# Keep each genomic/plasmid position as one angular bin.
# Increasing this above 1 makes peaks visually wider and can obscure exact
# cleavage positions.
BAR_WIDTH_MULTIPLIER = 3

# Keep None for automatic scaling.
# To compare multiple samples visually, set the SAME numerical RADIAL_MAX
# in every plot.
RADIAL_MAX = None


# -----------------------------------------------------------------
# OPTIONAL PLOT-ONLY MASKS
# -----------------------------------------------------------------
# If a peak is irrelevant for a particular figure and visually hides smaller
# peaks, put its 1-based CSV "Position" here. The original CSV is NOT changed.
#
# Examples:
# TOP_POSITIONS_TO_HIDE_1BASED = [337, 808]
# BOTTOM_POSITIONS_TO_HIDE_1BASED = [137, 171, 172, 173, 174, 999]
#
# Leave a list empty if you do not want to hide anything on that strand.
TOP_POSITIONS_TO_HIDE_1BASED = [272, 337, 808, 1111, 2304, 2705, 2727, 2808, 3216, 3227]
BOTTOM_POSITIONS_TO_HIDE_1BASED = []


# Image resolution.
DPI = 600


# =============================================================================
# LOAD DATA
# =============================================================================

def load_table(path):
    ext = os.path.splitext(path)[1].lower()

    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)

    if ext in [".csv", ".txt", ".tsv"]:
        # sep=None lets pandas detect comma / semicolon / tab automatically.
        return pd.read_csv(path, sep=None, engine="python")

    raise ValueError(
        f"Unsupported input format: {ext}. Use CSV, TSV, TXT, XLSX, or XLS."
    )


df = load_table(INPUT_FILE)

# Support both the current script's underscore column names and the older
# manually prepared Excel column names.
if USE_BPM:
    top_candidates = ["BPM_top", "Top_BPM", "BPM top"]
    bottom_candidates = ["BPM_bottom", "Bottom_BPM", "BPM bottom"]
else:
    top_candidates = ["Reads_top", "Reads top", "Top_Cleavage_Count"]
    bottom_candidates = ["Reads_bottom", "Reads bottom", "Bottom_Cleavage_Count"]


def find_column(candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        "Could not find any of these columns: "
        + ", ".join(candidates)
        + "\nAvailable columns are: "
        + ", ".join(map(str, df.columns))
    )


top_col = find_column(top_candidates)
bottom_col = find_column(bottom_candidates)

y_top = pd.to_numeric(df[top_col], errors="coerce").fillna(0).to_numpy(float)
y_bottom = pd.to_numeric(df[bottom_col], errors="coerce").fillna(0).to_numpy(float)

if len(y_top) != len(y_bottom):
    raise ValueError("Top and bottom arrays do not have the same number of positions.")

N = len(y_top)

if N == 0:
    raise ValueError("Input table contains no positions.")


# =============================================================================
# OPTIONAL PLOT-ONLY MASKING
# =============================================================================

# Work on copies so the values read from the file remain conceptually intact.
y_top = y_top.copy()
y_bottom = y_bottom.copy()

def hide_positions(values, positions_1based, strand_name):
    for position1 in positions_1based:
        i = position1 - 1

        if i < 0 or i >= len(values):
            print(
                f"WARNING: {strand_name} hide position {position1} is outside "
                f"the range 1..{len(values)}"
            )
            continue

        old_value = values[i]
        values[i] = 0.0
        print(
            f"Plot-only mask: {strand_name} position {position1} "
            f"changed from {old_value} to 0 for this figure."
        )

hide_positions(
    y_top,
    TOP_POSITIONS_TO_HIDE_1BASED,
    "TOP"
)
hide_positions(
    y_bottom,
    BOTTOM_POSITIONS_TO_HIDE_1BASED,
    "BOTTOM"
)


# =============================================================================
# POLAR GEOMETRY
# =============================================================================

# One exact angular bin per plasmid position.
theta = np.arange(N) * (2.0 * np.pi / N)
bin_width = (2.0 * np.pi / N) * BAR_WIDTH_MULTIPLIER

# Use a positive radial baseline rather than negative polar radii.
# This makes top extend outward and bottom extend inward without relying on
# Matplotlib's treatment of negative polar coordinates.
data_max = max(float(np.max(y_top)), float(np.max(y_bottom)))

if RADIAL_MAX is None:
    radial_max = data_max
else:
    radial_max = float(RADIAL_MAX)
    if radial_max <= 0:
        raise ValueError("RADIAL_MAX must be positive.")
    if data_max > radial_max:
        print(
            f"WARNING: data maximum ({data_max}) exceeds RADIAL_MAX "
            f"({radial_max}); the tallest bars will be clipped."
        )

# Avoid a zero-radius plot for an all-zero sample.
if radial_max == 0:
    radial_max = 1.0

baseline = radial_max


# =============================================================================
# DRAW
# =============================================================================

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})

# Explicitly define orientation:
# position 1 at 12 o'clock, increasing positions clockwise.
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# TOP: outward from baseline.
ax.bar(
    theta,
    y_top,
    width=bin_width,
    bottom=baseline,
    align="edge",
    color="red",
    linewidth=0,
)

# BOTTOM: inward from the SAME baseline, but using only positive radii.
ax.bar(
    theta,
    y_bottom,
    width=bin_width,
    bottom=baseline - y_bottom,
    align="edge",
    color="black",
    linewidth=0,
)

# Equal linear scale for top and bottom.
ax.set_ylim(0, 2 * radial_max)

# Optional thin circle at the top/bottom dividing baseline.
circle_theta = np.linspace(0, 2 * np.pi, 1000)
ax.plot(circle_theta, np.full_like(circle_theta, baseline),
        color="0.5", linewidth=0.4)

# Clean appearance.
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.spines["polar"].set_visible(False)

plt.tight_layout(pad=0)

# PNG for figures/slides, PDF for vector-quality editing/publication.
plt.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches="tight", pad_inches=0)
plt.savefig(OUTPUT_PDF, bbox_inches="tight", pad_inches=0)

print(f"Input: {INPUT_FILE}")
print(f"Top column: {top_col}")
print(f"Bottom column: {bottom_col}")
print(f"Positions: {N}")
print(f"Maximum plotted value: {radial_max}")
print(f"Saved: {OUTPUT_PNG}")
print(f"Saved: {OUTPUT_PDF}")

plt.show()
