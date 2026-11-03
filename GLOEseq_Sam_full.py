import os
import sys
import gzip
import random
import csv
import difflib
from collections import defaultdict

def import_reads(file_name, number_of_reads):
    """Reads gzipped FASTQ files and optionally limits the number of reads."""
    reads = []
    try:
        with gzip.open(file_name, 'rt') as file:
            read_counter = 0
            while True:
                header = file.readline()
                if not header:
                    break
                sequence = file.readline().strip().upper()
                plus = file.readline()
                quality = file.readline()
                reads.append(sequence)
                read_counter += 1
                if number_of_reads != 0 and read_counter >= number_of_reads:
                    break
    except IOError as e:
        print(f"Error opening file: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
    return reads

def browse():
    """Returns the file names to be processed."""
    return ["/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_17_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_19_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_15_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_16_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_2_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_3_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_4_GLOESEQ_L00.fq.gz",
            "/Users/sam/Documents/Experiments/Experiments/python files/FastQ datas/GLOEseq08_09_26/202609041656_AE01-250401014_4P260523176US293232S2_B_runG58_6_GLOESEQ_L00.fq.gz",
            ]

### Reference Sequences Input

top_strand = (
    "cagtcctcagctaccctcagcctgacctcagcgtcacctcagcagcttgagtattctatagtgtcacctaaatagcttggcgtaatcatggtcatagctgtttcctgtgtgaaattgttatccgctcacaattccacacaacatacgagccggaagcataaagtgtaaagcctggggtgcctaatgagtgagctaactcacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatctgccaacgcgcggggagaggcggtttgcgtattgggcgctcttccgcttcctcgctcactgagtcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgaatccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttctgaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctggctccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgacgtctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtctcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtacttagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggaaattgtaagcgttaatattttgttaaaattcgcgttaaatttttgttaaatcagctcattttttaaccaataggccgaaatcggcaaaatcccttataaatcaaaagaatagaccgagatagggttgagtgttgttccagtttggaacaagaatccactattaaagaacgtggaatccaacgtcaaagggcgaaaaaccgtctatcagggcgatggcccactacgtgaaccatcaccctaatcaagttttttggggtcgaggtgccgtaaagcactaaatcggaaccctaaagggagcccccgatttagagcttgacggggaaagccggcgaacgtggcgagaaaggaagggaagaaagcgaaagtagcgggcgctagggcgctggcaagtgtagcggtcacgctgcgcgtaaccaccacacccgccgcgcttaatgcgccgctacagggcgcgtccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcgggcctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcgattaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgacggccagtgaattgtaatacgaacactatagggcgaattggcggccgcgatctgatcagatccacctcagcacttcctcagcgagcctcagcgga"
).upper()

bottom_strand = (
    "tccgctgaggctcgctgaggaagtgctgaggtggatctgatcagatcgcggccgccaattcgccctatagtgttcgtattacaattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctactttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggattccacgttctttaatagtggattcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaatttcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctaagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagccaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttcagaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggattcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgactcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggcagattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagctatttaggtgacactatagaatactcaagctgctgaggtgacgctgaggtcaggctgagggtagctgagga"
).upper()

# Generate circular mapping support
# =============================================================================
# CORRECTED GLOE-SEQ SETTINGS
# =============================================================================

# Use a longer 5' anchor than the original 30 nt.
# This is especially important for the repeat-rich synthetic cassette.
MATCH_K = 75

# Efficiently rescue a read with one sequencing mismatch, but only if the
# resulting mapping is unique.
ALLOW_ONE_MISMATCH = True
SEED_K = 25

# Normally 0. Change only if R1 contains known technical bases before the
# true plasmid-derived sequence.
READ_START_OFFSET = 0

# Ambiguous/repetitive reads are excluded instead of being silently assigned
# to the first occurrence.
REQUIRE_UNIQUE_ALIGNMENT = True

# Because these 20 samples do not contain an internal reference, also report
# break reads per million uniquely mapped / coordinate-assigned reads.
BPM_SCALE = 1_000_000.0


# =============================================================================
# BASIC SEQUENCE UTILITIES
# =============================================================================

def reverse_complement(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()


def circular_slice(seq, start0, length):
    """Return `length` bases starting at `start0` on a circular sequence."""
    n = len(seq)
    return "".join(seq[(start0 + j) % n] for j in range(length))


def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))


# =============================================================================
# BUILD A TRUE TOP <-> BOTTOM PHYSICAL COORDINATE MAP
# =============================================================================

def build_strand_coordinate_maps():
    """
    Align top_strand to reverse-complement(bottom_strand).

    This allows the two references to have different lengths because of an
    extrahelical loop.

    Returns
    -------
    top_to_bottom_native:
        top index -> physically paired bottom index in bottom native 5'->3'
        coordinates; None for an unpaired top-loop base.

    bottom_native_to_top:
        bottom native index -> physically paired top index.
    """
    bottom_rc = reverse_complement(bottom_strand)

    top_to_bottom_native = [None] * len(top_strand)
    bottom_native_to_top = [None] * len(bottom_strand)

    matcher = difflib.SequenceMatcher(
        None, top_strand, bottom_rc, autojunk=False
    )

    for tag, top_start, top_end, rc_start, rc_end in matcher.get_opcodes():
        if tag != "equal":
            continue

        for offset in range(top_end - top_start):
            top_i = top_start + offset
            bottom_rc_i = rc_start + offset
            bottom_native_i = len(bottom_strand) - 1 - bottom_rc_i

            top_to_bottom_native[top_i] = bottom_native_i
            bottom_native_to_top[bottom_native_i] = top_i

    return top_to_bottom_native, bottom_native_to_top


TOP_TO_BOTTOM_NATIVE, BOTTOM_NATIVE_TO_TOP = build_strand_coordinate_maps()


def validate_strand_map():
    """
    Check that paired regions are complementary and report unpaired loop bases.
    """
    unpaired_top = [
        i + 1 for i, x in enumerate(TOP_TO_BOTTOM_NATIVE) if x is None
    ]
    unmapped_bottom = [
        i + 1 for i, x in enumerate(BOTTOM_NATIVE_TO_TOP) if x is None
    ]

    mismatches = []

    for top_i, bottom_i in enumerate(TOP_TO_BOTTOM_NATIVE):
        if bottom_i is None:
            continue

        expected_top_base = reverse_complement(bottom_strand[bottom_i])
        if top_strand[top_i] != expected_top_base:
            mismatches.append((top_i + 1, bottom_i + 1))

    print("\nReference / duplex check")
    print("------------------------")
    print(f"Top strand length:     {len(top_strand)}")
    print(f"Bottom strand length:  {len(bottom_strand)}")
    print(f"Unpaired top positions:{unpaired_top}")
    print(f"Unmapped bottom bases: {unmapped_bottom}")
    print(f"Paired-base mismatches:{len(mismatches)}")

    if unmapped_bottom:
        raise ValueError(
            "Some bottom-strand bases cannot be mapped to the top strand. "
            "Please re-check the reference sequences."
        )

    if mismatches:
        raise ValueError(
            "Paired parts of the top/bottom references are not complementary. "
            "Please re-check the reference sequences."
        )

    return unpaired_top


# =============================================================================
# UNIQUE CIRCULAR MAPPING OF THE 5' GLOE READ ANCHOR
# =============================================================================

REFS = {
    "top_ref": top_strand,
    "bottom_ref": bottom_strand,
}


def build_full_kmer_index(k):
    index = defaultdict(list)

    for ref_name, seq in REFS.items():
        for start0 in range(len(seq)):
            kmer = circular_slice(seq, start0, k)
            index[kmer].append((ref_name, start0))

    return dict(index)


def build_seed_index(seed_k):
    index = defaultdict(list)

    for ref_name, seq in REFS.items():
        for start0 in range(len(seq)):
            seed = circular_slice(seq, start0, seed_k)
            index[seed].append((ref_name, start0))

    return dict(index)


print("Building circular unique-mapping indices...")
FULL_INDEX = build_full_kmer_index(MATCH_K)
SEED_INDEX = build_seed_index(SEED_K)


def find_unique_alignment(read_fragment):
    """
    Map the first MATCH_K bases of a read.

    Returns
    -------
    ref_name, start0, mismatches, status

    status is one of:
      exact
      one_mismatch
      ambiguous
      orphan
      short_read
      invalid_base
    """
    if len(read_fragment) < MATCH_K:
        return None, None, None, "short_read"

    read_fragment = read_fragment[:MATCH_K]

    if any(base not in "ACGT" for base in read_fragment):
        return None, None, None, "invalid_base"

    # First try an exact full-length anchor.
    exact_hits = sorted(set(FULL_INDEX.get(read_fragment, [])))

    if exact_hits:
        if REQUIRE_UNIQUE_ALIGNMENT and len(exact_hits) != 1:
            return None, None, 0, "ambiguous"

        ref_name, start0 = exact_hits[0]
        return ref_name, start0, 0, "exact"

    if not ALLOW_ONE_MISMATCH:
        return None, None, None, "orphan"

    # Efficient one-mismatch rescue:
    # nominate positions from exact short seeds, then verify the full MATCH_K
    # anchor with Hamming distance.
    candidate_hits = set()

    for offset in range(0, MATCH_K, SEED_K):
        seed = read_fragment[offset:offset + SEED_K]

        if len(seed) != SEED_K:
            continue

        for ref_name, seed_start0 in SEED_INDEX.get(seed, []):
            seq_len = len(REFS[ref_name])
            candidate_start0 = (seed_start0 - offset) % seq_len
            candidate_hits.add((ref_name, candidate_start0))

    verified = []

    for ref_name, start0 in candidate_hits:
        ref_anchor = circular_slice(REFS[ref_name], start0, MATCH_K)
        distance = hamming_distance(read_fragment, ref_anchor)

        if distance <= 1:
            verified.append((ref_name, start0, distance))

    if not verified:
        return None, None, None, "orphan"

    best_distance = min(x[2] for x in verified)
    best_hits = sorted(
        set((ref_name, start0, distance)
            for ref_name, start0, distance in verified
            if distance == best_distance)
    )

    if REQUIRE_UNIQUE_ALIGNMENT and len(best_hits) != 1:
        return None, None, best_distance, "ambiguous"

    ref_name, start0, distance = best_hits[0]
    return ref_name, start0, distance, "one_mismatch"


# =============================================================================
# GLOE DIRECT-MODE CLEAVAGE POSITION
# =============================================================================

def cleavage_from_alignment(ref_name, mapped_start0):
    """
    Convert the mapped 5' end of the GLOE read into the original 3'-OH.

    GLOE reads are reverse complements of the captured DNA fragment.

      read maps bottom_ref -> original 3'-OH on TOP
      read maps top_ref    -> original 3'-OH on BOTTOM

    The returned nucleotide is the original 3'-terminal nucleotide itself
    (GLOE direct mode). The physical phosphodiester break is immediately 3'
    of that nucleotide on the cleaved strand.
    """
    if ref_name == "bottom_ref":
        top_i = BOTTOM_NATIVE_TO_TOP[mapped_start0]

        if top_i is None:
            return None

        return {
            "strand": "top",
            "native_index0": top_i,
            "plot_index0": top_i,
        }

    if ref_name == "top_ref":
        bottom_native_i = TOP_TO_BOTTOM_NATIVE[mapped_start0]

        # If the read begins exactly on an unpaired top-loop base, there is no
        # physically paired bottom nucleotide to assign in direct mode.
        if bottom_native_i is None:
            return None

        plot_i = BOTTOM_NATIVE_TO_TOP[bottom_native_i]

        return {
            "strand": "bottom",
            "native_index0": bottom_native_i,
            "plot_index0": plot_i,
        }

    raise ValueError(f"Unknown reference: {ref_name}")


# =============================================================================
# READ ANALYSIS
# =============================================================================

def analyse_reads(reads):
    top_counts = [0] * len(top_strand)
    bottom_counts = [0] * len(bottom_strand)

    orphan_reads = []
    ambiguous_reads = []
    coordinate_unmappable_reads = []
    short_reads = []
    invalid_base_reads = []
    read_alignments = []

    exact_aligned = 0
    one_mismatch_aligned = 0

    for read in reads:
        read_fragment = read[
            READ_START_OFFSET:READ_START_OFFSET + MATCH_K
        ]

        ref_name, mapped_start0, mismatches, status = (
            find_unique_alignment(read_fragment)
        )

        if status == "short_read":
            short_reads.append(read)
            continue

        if status == "invalid_base":
            invalid_base_reads.append(read)
            continue

        if status == "ambiguous":
            ambiguous_reads.append(read)
            continue

        if status == "orphan":
            orphan_reads.append(read)
            continue

        cleavage = cleavage_from_alignment(ref_name, mapped_start0)

        if cleavage is None:
            coordinate_unmappable_reads.append(read)
            continue

        strand = cleavage["strand"]
        native_i = cleavage["native_index0"]
        plot_i = cleavage["plot_index0"]

        if strand == "top":
            top_counts[native_i] += 1
        else:
            bottom_counts[native_i] += 1

        if status == "exact":
            exact_aligned += 1
        elif status == "one_mismatch":
            one_mismatch_aligned += 1

        read_alignments.append({
            "read_sequence": read,
            "mapped_reference": ref_name,
            "mapped_start_native_1based": mapped_start0 + 1,
            "alignment_status": status,
            "mismatches": mismatches,
            "cleavage_strand": strand,
            "cleavage_position_native_1based": native_i + 1,
            "cleavage_position_plot_1based": plot_i + 1,
        })

    qc = {
        "exact_aligned": exact_aligned,
        "one_mismatch_aligned": one_mismatch_aligned,
    }

    return (
        top_counts,
        bottom_counts,
        orphan_reads,
        ambiguous_reads,
        coordinate_unmappable_reads,
        short_reads,
        invalid_base_reads,
        read_alignments,
        qc,
    )


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def bpm(count, mapped_total):
    if mapped_total <= 0:
        return 0.0

    return (count / mapped_total) * BPM_SCALE


def bottom_at_plot_position(plot_i, bottom_counts):
    """
    Return bottom information physically aligned with a top/clockwise coordinate.
    """
    bottom_native_i = TOP_TO_BOTTOM_NATIVE[plot_i]

    if bottom_native_i is None:
        return {
            "base": "-",
            "count": 0,
            "native_position": "",
            "is_unpaired": True,
        }

    return {
        "base": bottom_strand[bottom_native_i],
        "count": bottom_counts[bottom_native_i],
        "native_position": bottom_native_i + 1,
        "is_unpaired": False,
    }


def write_outputs(
    folder_name,
    reads,
    top_counts,
    bottom_counts,
    orphan_reads,
    ambiguous_reads,
    coordinate_unmappable_reads,
    short_reads,
    invalid_base_reads,
    read_alignments,
    qc,
    unpaired_top_positions,
):
    total_reads = len(reads)
    top_aligned_reads = sum(top_counts)
    bottom_aligned_reads = sum(bottom_counts)
    mapped_total = top_aligned_reads + bottom_aligned_reads

    # -------------------------------------------------------------------------
    # top_cleavage_sites.csv
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "top_cleavage_sites.csv"),
        "w",
        newline=""
    ) as f_top:
        writer = csv.writer(f_top)
        writer.writerow(["Position", "Sequence", "Read_Count", "BPM"])

        for i in range(len(top_strand)):
            writer.writerow([
                i + 1,
                top_strand[i],
                top_counts[i],
                f"{bpm(top_counts[i], mapped_total):.6f}",
            ])

    # -------------------------------------------------------------------------
    # bottom_cleavage_sites.csv
    #
    # This file is already physically oriented to the same top/clockwise
    # coordinate. Do NOT reverse it again.
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "bottom_cleavage_sites.csv"),
        "w",
        newline=""
    ) as f_bottom:
        writer = csv.writer(f_bottom)
        writer.writerow([
            "Position",
            "Sequence",
            "Read_Count",
            "BPM",
            "Bottom_Native_Position_1based",
            "Is_Unpaired_Top_Loop_Position",
        ])

        for plot_i in range(len(top_strand)):
            info = bottom_at_plot_position(plot_i, bottom_counts)

            writer.writerow([
                plot_i + 1,
                info["base"],
                info["count"],
                f"{bpm(info['count'], mapped_total):.6f}",
                info["native_position"],
                info["is_unpaired"],
            ])

    # -------------------------------------------------------------------------
    # nick_position.csv
    #
    # This is the corrected physical coordinate table.
    # Both strands are aligned to the same top/clockwise plasmid coordinate.
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "nick_position.csv"),
        "w",
        newline=""
    ) as f_nick:
        writer = csv.writer(f_nick)
        writer.writerow([
            "Position",
            "Top_Sequence",
            "Top_Cleavage_Count",
            "Bottom_Sequence_3to5",
            "Bottom_Cleavage_Count",
            "Top_BPM",
            "Bottom_BPM",
            "Bottom_Native_Position_1based",
            "Is_Unpaired_Top_Loop_Position",
        ])

        for plot_i in range(len(top_strand)):
            info = bottom_at_plot_position(plot_i, bottom_counts)

            writer.writerow([
                plot_i + 1,
                top_strand[plot_i],
                top_counts[plot_i],
                info["base"],
                info["count"],
                f"{bpm(top_counts[plot_i], mapped_total):.6f}",
                f"{bpm(info['count'], mapped_total):.6f}",
                info["native_position"],
                info["is_unpaired"],
            ])

    # -------------------------------------------------------------------------
    # raw_nick_position_plot_ready.csv
    #
    # Simple file for your polar plot scripts.
    # Bottom is ALREADY in the correct physical orientation.
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "raw_nick_position_plot_ready.csv"),
        "w",
        newline=""
    ) as f_plot:
        writer = csv.writer(f_plot)
        writer.writerow([
            "Position",
            "Reads_top",
            "Reads_bottom",
            "Top_base_5to3",
            "Bottom_base_3to5",
            "BPM_top",
            "BPM_bottom",
            "Bottom_native_position_1based",
            "Is_unpaired_top_loop_position",
        ])

        for plot_i in range(len(top_strand)):
            info = bottom_at_plot_position(plot_i, bottom_counts)

            writer.writerow([
                plot_i + 1,
                top_counts[plot_i],
                info["count"],
                top_strand[plot_i],
                info["base"],
                f"{bpm(top_counts[plot_i], mapped_total):.6f}",
                f"{bpm(info['count'], mapped_total):.6f}",
                info["native_position"],
                info["is_unpaired"],
            ])

    # -------------------------------------------------------------------------
    # nick_position_native.csv
    #
    # Keeps each strand in its own native 5'->3' coordinate for debugging.
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "nick_position_native.csv"),
        "w",
        newline=""
    ) as f_native:
        writer = csv.writer(f_native)
        writer.writerow([
            "Native_Position_1based",
            "Top_Base_5to3",
            "Top_Count",
            "Bottom_Base_5to3_Native",
            "Bottom_Count_Native",
        ])

        for i in range(max(len(top_strand), len(bottom_strand))):
            writer.writerow([
                i + 1,
                top_strand[i] if i < len(top_strand) else "",
                top_counts[i] if i < len(top_counts) else "",
                bottom_strand[i] if i < len(bottom_strand) else "",
                bottom_counts[i] if i < len(bottom_counts) else "",
            ])

    # -------------------------------------------------------------------------
    # results_summary.txt
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "results_summary.txt"),
        "w"
    ) as f_summary:

        def pct(value):
            if total_reads == 0:
                return 0.0
            return value * 100.0 / total_reads

        f_summary.write(f"Number of reads: {total_reads}\n")
        f_summary.write(
            f"Top strand cleavage events: "
            f"{top_aligned_reads} ({pct(top_aligned_reads):.2f}%)\n"
        )
        f_summary.write(
            f"Bottom strand cleavage events: "
            f"{bottom_aligned_reads} ({pct(bottom_aligned_reads):.2f}%)\n"
        )
        f_summary.write(
            f"Unique coordinate-assigned mapped reads: "
            f"{mapped_total} ({pct(mapped_total):.2f}%)\n"
        )
        f_summary.write(
            f"Exact aligned: {qc['exact_aligned']} "
            f"({pct(qc['exact_aligned']):.2f}%)\n"
        )
        f_summary.write(
            f"One-mismatch aligned: {qc['one_mismatch_aligned']} "
            f"({pct(qc['one_mismatch_aligned']):.2f}%)\n"
        )
        f_summary.write(
            f"Ambiguous reads: {len(ambiguous_reads)} "
            f"({pct(len(ambiguous_reads)):.2f}%)\n"
        )
        f_summary.write(
            f"Orphan reads: {len(orphan_reads)} "
            f"({pct(len(orphan_reads)):.2f}%)\n"
        )
        f_summary.write(
            f"Coordinate-unmappable loop-start reads: "
            f"{len(coordinate_unmappable_reads)} "
            f"({pct(len(coordinate_unmappable_reads)):.2f}%)\n"
        )
        f_summary.write(
            f"Short reads: {len(short_reads)} "
            f"({pct(len(short_reads)):.2f}%)\n"
        )
        f_summary.write(
            f"Reads with N/invalid bases in anchor: "
            f"{len(invalid_base_reads)} "
            f"({pct(len(invalid_base_reads)):.2f}%)\n"
        )

        f_summary.write("\nReference information\n")
        f_summary.write("---------------------\n")
        f_summary.write(f"Top length: {len(top_strand)}\n")
        f_summary.write(f"Bottom length: {len(bottom_strand)}\n")
        f_summary.write(
            f"Unpaired top-loop positions: {unpaired_top_positions}\n"
        )

        f_summary.write("\nMapping settings\n")
        f_summary.write("----------------\n")
        f_summary.write(f"MATCH_K: {MATCH_K}\n")
        f_summary.write(
            f"ALLOW_ONE_MISMATCH: {ALLOW_ONE_MISMATCH}\n"
        )
        f_summary.write(
            f"READ_START_OFFSET: {READ_START_OFFSET}\n"
        )

        f_summary.write("\nGLOE direct-mode convention\n")
        f_summary.write("---------------------------\n")
        f_summary.write(
            "Read maps to bottom -> original 3'-OH assigned to TOP.\n"
        )
        f_summary.write(
            "Read maps to top -> original 3'-OH assigned to BOTTOM.\n"
        )
        f_summary.write(
            "Reported position = nucleotide carrying the original 3'-OH.\n"
        )
        f_summary.write(
            "Physical phosphodiester break = immediately 3' of that nucleotide.\n"
        )

        f_summary.write("\nNormalisation\n")
        f_summary.write("-------------\n")
        f_summary.write(
            "BPM = count / uniquely mapped coordinate-assigned reads * 1,000,000.\n"
        )
        f_summary.write(
            "These samples have no internal standard, so BPM corrects sequencing "
            "depth but is NOT absolute cleavage quantification.\n"
        )

    # -------------------------------------------------------------------------
    # read_alignments.csv
    # -------------------------------------------------------------------------
    with open(
        os.path.join(folder_name, "read_alignments.csv"),
        "w",
        newline=""
    ) as f_alignments:
        writer = csv.writer(f_alignments)
        writer.writerow([
            "Read_Sequence",
            "Mapped_Reference",
            "Mapped_Start_Native_1based",
            "Alignment_Status",
            "Mismatches",
            "Cleavage_Strand",
            "Cleavage_Position_Native_1based",
            "Cleavage_Position_Plot_1based",
        ])

        for a in read_alignments:
            writer.writerow([
                a["read_sequence"],
                a["mapped_reference"],
                a["mapped_start_native_1based"],
                a["alignment_status"],
                a["mismatches"],
                a["cleavage_strand"],
                a["cleavage_position_native_1based"],
                a["cleavage_position_plot_1based"],
            ])

    # -------------------------------------------------------------------------
    # excluded reads
    # -------------------------------------------------------------------------
    excluded_files = [
        ("orphan_reads.txt", orphan_reads),
        ("ambiguous_reads.txt", ambiguous_reads),
        (
            "coordinate_unmappable_reads.txt",
            coordinate_unmappable_reads
        ),
        ("short_reads.txt", short_reads),
        ("invalid_base_reads.txt", invalid_base_reads),
    ]

    for filename, excluded_reads in excluded_files:
        with open(os.path.join(folder_name, filename), "w") as fh:
            for read in excluded_reads:
                fh.write(read + "\n")

    return {
        "total_reads": total_reads,
        "top_aligned": top_aligned_reads,
        "bottom_aligned": bottom_aligned_reads,
        "mapped_total": mapped_total,
        "exact": qc["exact_aligned"],
        "one_mismatch": qc["one_mismatch_aligned"],
        "ambiguous": len(ambiguous_reads),
        "orphan": len(orphan_reads),
        "coordinate_unmappable": len(coordinate_unmappable_reads),
        "short_reads": len(short_reads),
        "invalid_base_reads": len(invalid_base_reads),
    }


# =============================================================================
# PROGRAM STARTS HERE
# =============================================================================

def main():
    unpaired_top_positions = validate_strand_map()

    file_names = browse()

    if not file_names:
        raise SystemExit("No FASTQ files are listed in browse().")

    number_of_reads = int(input("How many reads[0]? "))

    combined_summary = []

    for file_name in file_names:
        if not os.path.exists(file_name):
            print(f"WARNING: file not found, skipping: {file_name}")
            continue

        folder_name = os.path.join(
            os.path.dirname(file_name),
            f"{os.path.basename(file_name)[:-6]}_"
            f"{random.randint(0, 100000)}"
        )
        os.makedirs(folder_name, exist_ok=True)

        print("\nAnalyzing: " + os.path.basename(file_name))

        reads = import_reads(file_name, number_of_reads)
        print("Number of reads: " + str(len(reads)))

        (
            top_counts,
            bottom_counts,
            orphan_reads,
            ambiguous_reads,
            coordinate_unmappable_reads,
            short_reads,
            invalid_base_reads,
            read_alignments,
            qc,
        ) = analyse_reads(reads)

        summary = write_outputs(
            folder_name,
            reads,
            top_counts,
            bottom_counts,
            orphan_reads,
            ambiguous_reads,
            coordinate_unmappable_reads,
            short_reads,
            invalid_base_reads,
            read_alignments,
            qc,
            unpaired_top_positions,
        )

        print(
            f"Top strand cleavage events: "
            f"{summary['top_aligned']}"
        )
        print(
            f"Bottom strand cleavage events: "
            f"{summary['bottom_aligned']}"
        )
        print(
            f"Unique mapped / coordinate-assigned reads: "
            f"{summary['mapped_total']}"
        )
        print(f"Ambiguous reads: {summary['ambiguous']}")
        print(f"Orphan reads: {summary['orphan']}")
        print("Output folder:", folder_name)

        combined_summary.append({
            "FASTQ": file_name,
            "Output_Folder": folder_name,
            **summary,
        })

    # One overview CSV for all samples processed in the same run.
    if combined_summary:
        parent = os.path.dirname(combined_summary[0]["FASTQ"])
        combined_path = os.path.join(
            parent,
            f"GLOEseq_all_samples_mapping_summary_"
            f"{random.randint(0, 100000)}.csv"
        )

        fields = [
            "FASTQ",
            "Output_Folder",
            "total_reads",
            "top_aligned",
            "bottom_aligned",
            "mapped_total",
            "exact",
            "one_mismatch",
            "ambiguous",
            "orphan",
            "coordinate_unmappable",
            "short_reads",
            "invalid_base_reads",
        ]

        with open(combined_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()

            for row in combined_summary:
                writer.writerow(row)

        print("\nCombined mapping summary:", combined_path)


if __name__ == "__main__":
    main()
