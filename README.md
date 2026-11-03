# SamCutFinder

**Strand-resolved mapping and circular visualization of GLOE-seq cleavage sites on plasmid DNA**

SamCutFinder consists of two Python scripts for analysing plasmid-based GLOE-seq data and visualizing the positions of captured DNA 3′-OH ends on circular DNA substrates. The first script assigns reads to strand-specific plasmid coordinates and produces per-position count tables. The second creates polar plots showing the distribution of those counts on both strands.

> **Scope:** This is a research script collection, not a general-purpose sequencing aligner or a validated end-to-end software package. The analysis script contains a specific pair of plasmid reference sequences. Users must verify or replace these sequences for their experiments and validate the mapping/strand conventions before drawing biological conclusions.

## Repository contents

| File | Purpose |
| --- | --- |
| `GLOEseq_Sam_full.py` | Read mapping, strand assignment, cleavage-position counting, QC summaries, and CSV export. |
| `Polar_Plotting_GLOE_with_optional_manual_masks(1).py` | Circular plotting of the strand-resolved count tables. |
| `requirements.txt` | Third-party packages needed by the plotting script. |
| `CITATION.cff` | Machine-readable instructions for citing the software on GitHub. |
| `LICENSE` | MIT license. |

## How the analysis works

1. Reads are loaded from **gzipped FASTQ** (`.fq.gz` or `.fastq.gz`) files. The analysis uses read sequences but not per-base quality scores; FASTQ record structure and sequence/quality lengths are checked.
2. The first **75 nt** of each read, starting at `READ_START_OFFSET = 0`, are compared against both user-specified circular plasmid strand references. Exact matches are attempted first.
3. When `ALLOW_ONE_MISMATCH = True`, alignments with one mismatch in the 75-nt anchor can be rescued using candidate 25-nt seeds. With `REQUIRE_UNIQUE_ALIGNMENT = True`, ambiguous best alignments are excluded.
4. The two reference strands are aligned using their complementary sequences so that strand-specific counts share a physical coordinate system. Fully complementary circular plasmids without an extrahelical loop are supported. Where extrahelical positions occur on the top reference strand, they are flagged explicitly.
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

### A. Select the input FASTQ files (no Mac-specific code paths)

The updated analysis script accepts a directory containing gzipped FASTQ files, or explicit file paths. It does **not** require editing the `browse()` function.

Process every `.fq.gz` or `.fastq.gz` file directly inside a folder (not recursively):

```bash
python GLOEseq_Sam_full.py --input-dir "/path/to/your/FASTQ_folder"
```

Alternatively, process one file:

```bash
python GLOEseq_Sam_full.py "/path/to/sample_01.fq.gz"
```

Or multiple named files:

```bash
python GLOEseq_Sam_full.py "/path/to/sample_01.fq.gz" "/path/to/sample_02.fq.gz"
```

For a short test run of 10,000 reads per sample, add `--max-reads 10000`. The default (`--max-reads 0`) processes all reads. Use `python GLOEseq_Sam_full.py --help` for the full command-line help.

### B. Set the correct plasmid references

The `top_strand` and `bottom_strand` variables are **embedded directly in the script**. The current references have 3,248 and 3,245 nucleotides, respectively, and the code reports three unpaired top-strand positions for that substrate. **These are example-specific GLOE-seq reference sequences, not a required plasmid architecture.** Replace both sequences with the correct strand sequences for each substrate and check their 5′→3′ orientation and nucleotide numbering. A plasmid with no extrahelical loop can be analysed with fully complementary references. For substrates with extrahelical positions, ensure the strand-pairing map corresponds to the DNA construct.

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

The FASTQ reads must contain the relevant plasmid-derived sequence at the expected offset. This script does not automatically trim adapters, remove UMIs, or filter by sequencing quality. It rejects malformed FASTQ records rather than silently processing partial inputs. If an experiment uses a different read layout or capture protocol, validate and adapt these settings and the cleavage-coordinate convention.

### D. Run the script

Open Terminal in the repository folder and execute the command in Part A, substituting your actual FASTQ folder/path. For example:

```bash
python GLOEseq_Sam_full.py --input-dir "/Users/yourname/Documents/GLOEseq_FASTQ"
```

The updated version no longer prompts for a read count. All reads are processed by default. For each FASTQ, the script creates an output folder **next to the FASTQ**, with a sample-derived name and random numeric suffix. A combined mapping summary CSV is also written in the directory containing the first processed sample.

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

1. Locate `raw_nick_position_plot_ready.csv` from the corresponding sample's analysis output folder.
2. From the repository folder, run this command (replace the path with your actual CSV):

   ```bash
   python "Polar_Plotting_GLOE_with_optional_manual_masks(1).py" "/path/to/raw_nick_position_plot_ready.csv"
   ```

3. By default the script plots **BPM**, and generates **`gloe_polar_plot.png`** and **`gloe_polar_plot.pdf`** in the current working directory.
4. To plot raw read counts rather than BPM, add `--raw`.
5. For the same radial scale in compared plots, add `--radial-max 1000` (replace 1000 with an appropriate fixed value). To choose an output filename, use `--output-prefix "/path/to/sample1_circle"`. If `--radial-max` is less than the largest peak, bars exceeding that fixed value are clipped in the image (with a warning); the input data remain unchanged.

**No cleavage positions are hidden by default** in this updated script:

```python
TOP_POSITIONS_TO_HIDE_1BASED = []
BOTTOM_POSITIONS_TO_HIDE_1BASED = []
```

Optional display-only masking can be requested explicitly using `--hide-top 272 337` or `--hide-bottom 137 171`, but **do not use masked representations as an undisclosed depiction of the complete cleavage profile**. The input CSV is not altered. If masking is ever scientifically justified, document the positions, rationale, and unmasked result.

**Visual resolution:** `BAR_WIDTH_MULTIPLIER = 3` in the original plot style makes bars three times the angular spacing of adjacent positions; it may visually merge adjacent peaks. Set this constant to `1` when an exact one-position-wide representation is required. Maintain the same scale and graphical parameters across compared samples.

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
- Optional plotting masks are **visualization-only** and must not be confused with exclusions applied during read mapping or quantification.

## Reproducibility and citation

The main repository is available at <https://github.com/Sam18-hub/SamCutFinder>.

For a manuscript, please cite a **specific versioned release** (ideally archived with a Zenodo DOI), rather than only the moving `main` branch. Record the code version, Python/dependency versions, mapping settings, reference sequences, any display masks, and the relevant figure-generating CSVs. A `CITATION.cff` file is provided for GitHub's **“Cite this repository”** feature.

## Author and license

**Issam Senoussi** — Institute for Research in Biomedicine (IRB), Bellinzona, Switzerland.

Distributed under the **MIT License**; see [`LICENSE`](LICENSE).
