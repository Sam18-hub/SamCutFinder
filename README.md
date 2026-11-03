# SamCutFinder

**Strand-resolved mapping and circular visualization of GLOE-seq cleavage sites on plasmid DNA**

SamCutFinder consists of two Python scripts for analysing plasmid-based GLOE-seq data and visualizing the positions of captured DNA 3′-OH ends on circular DNA substrates. The first script assigns reads to strand-specific plasmid coordinates and produces per-position count tables. The second creates polar plots showing the distribution of those counts on both strands.

> **Scope:** This is a research script collection, not a general-purpose sequencing aligner or a validated end-to-end software package. The current analysis script contains a specific pair of plasmid reference sequences and sample-specific FASTQ paths. Users must configure these for their experiments and validate the mapping/strand conventions before drawing biological conclusions.

## Repository contents

| File | Purpose |
| --- | --- |
| `GLOEseq_Sam_full.py` | Read mapping, strand assignment, cleavage-position counting, QC summaries, and CSV export. |
| `Polar_Plotting_GLOE_with_optional_manual_masks.py` | Circular plotting of the strand-resolved count tables. |
| `requirements.txt` | Third-party packages needed by the plotting script. |
| `CITATION.cff` | Machine-readable instructions for citing the software on GitHub. |
| `LICENSE` | MIT license. |

## How the analysis works

1. Reads are loaded from **gzipped FASTQ** (`.fq.gz` or `.fastq.gz`) files. The analysis script currently uses the read sequence, not its per-base quality scores.
2. The first **75 nt** of each read, starting at `READ_START_OFFSET = 0`, are compared against both user-specified circular plasmid strand references. Exact matches are attempted first.
3. When `ALLOW_ONE_MISMATCH = True`, alignments with one mismatch in the 75-nt anchor can be rescued using candidate 25-nt seeds. With `REQUIRE_UNIQUE_ALIGNMENT = True`, ambiguous best alignments are excluded.
4. The reference strands are aligned to each other using their complementary sequence, allowing the output to put opposite-strand cleavage counts in a common physical coordinate system even when one strand contains an extrahelical loop. Unpaired positions are explicitly flagged.
5. Under the **direct-mode GLOE-seq convention implemented by the script**, a read aligned to the bottom reference is assigned to a 3′-OH on the top strand, and vice versa. The reported cleavage position is the **nucleotide carrying the original 3′-OH**, rather than the following phosphodiester bond.
6. The script reports read counts at each position and **BPM (break reads per million uniquely mapped, coordinate-assigned reads)**:

   `BPM = (reads assigned to the site / total uniquely mapped and coordinate-assigned reads) × 1,000,000`

   BPM normalizes for the number of mapped reads **within a library**. Without an internal standard, it is **not an absolute cleavage fraction or cleavage efficiency**. Reads are counted individually; the script does **not** perform UMI-based deduplication.

## Requirements

- Python 3.10 or newer is recommended.
- The analysis script uses Python's **standard library only**.
- The plotting script uses `numpy`, `pandas`, and `matplotlib`.
- `openpyxl` is included for optional `.xlsx` input and `xlrd` for optional `.xls` input.

These dependency ranges were determined from the imports; exact library versions for the published analysis should be recorded and tested on the authors' computational environment.

## Installation (macOS / Linux)

Open Terminal and run:

```bash
git clone https://github.com/Sam18-hub/SamCutFinder.git
cd SamCutFinder
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\\Scripts\\activate` and use `python` instead of `python3` where appropriate.

**Without Git installed:** On the repository page, click **Code → Download ZIP**, extract the folder, open Terminal in that folder, and run the remaining commands.

## Part 1 — Configure and run the FASTQ analysis

### A. Specify the input FASTQ files

Open `GLOEseq_Sam_full.py` in a text editor. Find the `browse()` function near the beginning of the file. **It currently contains eight absolute paths under `/Users/sam/...`**, which will not exist on another computer. Replace its return value with paths to your own gzipped FASTQ files, for example:

```python
def browse():
    return [
        "/path/to/sample_01.fq.gz",
        "/path/to/sample_02.fq.gz",
    ]
```

You can enter one sample or several. The paths must refer to existing files.

### B. Set the correct plasmid references

The `top_strand` and `bottom_strand` variables are **embedded directly in the script**. The current references have 3,248 and 3,245 nucleotides, respectively, and the code reports three unpaired top-strand positions for that substrate. **These are not universal GLOE-seq references.** Replace both sequences with the correct strand sequences for each substrate and check that their relative orientation, loop position, and nucleotide numbering match your experiment.

The program prints a reference/duplex check and rejects unmapped bottom bases or non-complementary paired positions. This automated check should **not** substitute for an independent validation of the actual plasmid sequence and strand assignment.

### C. Check mapping parameters

Near the `CORRECTED GLOE-SEQ SETTINGS` section, the current defaults are:

```python
MATCH_K = 75
ALLOW_ONE_MISMATCH = True
SEED_K = 25
READ_START_OFFSET = 0
REQUIRE_UNIQUE_ALIGNMENT = True
BPM_SCALE = 1_000_000.0
```

The FASTQ reads must contain the relevant plasmid-derived sequence at the expected offset. This script does not automatically trim adapters, remove UMIs, or filter by sequencing quality. If an experiment uses a different read layout or capture protocol, validate and adapt these settings and the cleavage-coordinate convention.

### D. Run the script

From the repository folder:

```bash
python GLOEseq_Sam_full.py
```

The script asks:

```text
How many reads[0]?
```

Type **`0` and press Return** to process all reads. To test on a subset, type a positive integer such as `10000`. **Do not simply press Return**: in the current version, a blank answer raises a `ValueError`.

For each input FASTQ, the script creates an output folder **next to that FASTQ file**, with a sample-derived name and random numeric suffix. It also writes a combined mapping summary CSV in the directory containing the first processed sample.

### E. Analysis output files

| File | Contents |
| --- | --- |
| `top_cleavage_sites.csv` | Top-strand position, sequence, raw count, BPM. |
| `bottom_cleavage_sites.csv` | Bottom counts placed in the shared top-oriented coordinate, including native bottom positions and unpaired-loop flags. |
| `nick_position.csv` | Combined physically aligned top/bottom positions with counts, BPM and bases. |
| **`raw_nick_position_plot_ready.csv`** | **Direct input to the circular plotting script:** `Position`, `Reads_top`, `Reads_bottom`, `BPM_top`, `BPM_bottom`, and reference annotations. |
| `nick_position_native.csv` | Top and bottom sequences/counts in their respective native 5′→3′ coordinate systems. |
| `results_summary.txt` | Read totals, alignment categories, mapping settings, and normalization convention. |
| `read_alignments.csv` | Read-level mapping, mismatch count and cleavage-coordinate assignments. |
| `orphan_reads.txt`, `ambiguous_reads.txt`, `coordinate_unmappable_reads.txt`, `short_reads.txt`, `invalid_base_reads.txt` | Reads excluded for the indicated reasons. |

**Read `results_summary.txt` before interpreting any plots**, especially the ambiguous, unmapped and coordinate-unmappable fractions.

## Part 2 — Generate a circular cleavage-site plot

1. Find the sample's `raw_nick_position_plot_ready.csv` from Part 1.
2. Open `Polar_Plotting_GLOE_with_optional_manual_masks.py` in a text editor.
3. Edit `INPUT_FILE` to point to that CSV. Use an absolute path if it is outside the repository folder:

   ```python
   INPUT_FILE = "/path/to/sample_output/raw_nick_position_plot_ready.csv"
   ```

4. Choose the quantity to plot:

   ```python
   USE_BPM = True  # False = raw read counts; True = BPM
   ```

5. **Important: review the position masks before generating scientific figures.** The uploaded plotting script currently contains a list of **10 top-strand positions to hide**. These signals are set to zero **in the plotted arrays only**, while the original CSV is unchanged. For an unmasked figure, replace both lists with:

   ```python
   TOP_POSITIONS_TO_HIDE_1BASED = []
   BOTTOM_POSITIONS_TO_HIDE_1BASED = []
   ```

   If you deliberately hide any positions in a publication figure, disclose the exact exclusion criteria, affected positions, and rationale, and preserve an unmasked visualization for transparency.

6. Optional graphical settings:

   ```python
   BAR_WIDTH_MULTIPLIER = 3  # visual width; >1 overlaps neighbouring angular bins
   RADIAL_MAX = None          # auto scaling; set a shared value across samples for comparisons
   DPI = 600
   ```

7. Run:

   ```bash
   python Polar_Plotting_GLOE_with_optional_manual_masks.py
   ```

By default, the script saves `gloe_polar_plot.png` and `gloe_polar_plot.pdf` in the **current working directory**, and opens a Matplotlib plot window if a graphical display is available.

### Reading the circular plot

- **Clockwise positions:** position 1 starts at 12 o'clock, with positions increasing clockwise.
- **Red bars:** top-strand signal extending outward from the central dividing circle.
- **Black bars:** bottom-strand signal extending inward from the same dividing circle.
- **Scale:** top and bottom share the same linear radial scale within each plot. For across-sample comparisons, use the same `RADIAL_MAX` and the same signal unit (`USE_BPM`).
- `BAR_WIDTH_MULTIPLIER = 3` makes bars wider than their original position spacing. This is a **visual display choice**, not pooling or smoothing of adjacent counts; individual positions may visually overlap.

The plotting script also accepts compatible CSV/TSV/TXT and Excel (`.xlsx`/`.xls`) tables with the expected column names. Its output is not interactive: it creates static PNG and PDF figures.

## Scientific assumptions and limitations

- This workflow targets the **specific circular plasmid substrates** supplied as reference sequences, rather than performing general genomic alignment.
- Mapping relies on the first 75 nt of the read (by default); confirm the read orientation, adapter structure and strand-capture chemistry for your particular GLOE-seq library.
- Reads with ambiguous best matches are excluded; repeat-rich regions may therefore have reduced mappability.
- The mapping accepts at most one mismatch in the anchor when rescue is enabled; it does not implement general gapped alignment or indel correction.
- The analysis does not use read-quality scores, deduplicate PCR molecules/UMIs, or estimate absolute cleavage efficiency.
- The reference-duplex correspondence is inferred using sequence matching; inspect the coordinate map carefully for alternative plasmids or more complex loop/indel architectures.
- Entire FASTQ read sequences and per-read output records are retained in memory, which can be demanding for very large libraries.
- The plotting masks are **visualization-only** and must not be confused with exclusions applied during read mapping or quantification.

## Reproducibility and citation

The main repository is available at <https://github.com/Sam18-hub/SamCutFinder>.

For a manuscript, please cite a **specific versioned release** (ideally archived with a Zenodo DOI), rather than only the moving `main` branch. Record the code version, Python/dependency versions, mapping settings, reference sequences, any display masks, and the relevant figure-generating CSVs. A `CITATION.cff` file is provided for GitHub's **“Cite this repository”** feature.

## Author and license

**Issam Senoussi** — Institute for Research in Biomedicine (IRB), Bellinzona, Switzerland.

Distributed under the **MIT License**; see [`LICENSE`](LICENSE).
