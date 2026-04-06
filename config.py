# src/config.py

"""
Configuration file for the WCSA DNA Sequencing Pipeline.

This module centralizes all tunable constants used across the repository,
including simulation settings, signal generation parameters, quality-control
thresholds, WCSA scoring constants, and plotting defaults.

Designed for the NRT Researcher-a-thon repository structure.
"""

from __future__ import annotations

import numpy as np

# ============================================================
# 1. GLOBAL REPRODUCIBILITY
# ============================================================

RANDOM_SEED: int = 42

# ============================================================
# 2. REFERENCE / READ GENERATION
# ============================================================

REF_LEN: int = 800
READ_LEN: int = 60
NUM_READS: int = 60

DNA_BASES: tuple[str, ...] = ("A", "C", "G", "T")

# ============================================================
# 3. RAW SIGNAL SIMULATION
# ============================================================

DWELL_LEN: int = 5
NOISE_STD: float = 8.0

# Nanopore-like current levels used in synthetic signal simulation
LEVELS: dict[str, int] = {
    "A": 60,
    "C": 75,
    "G": 90,
    "T": 105,
}

# Backward-compatible alias if other modules use this name
BASE_CURRENT = LEVELS

# ============================================================
# 4. QUALITY CONTROL
# ============================================================

QC_THRESHOLD: float = 0.15

# ============================================================
# 5. WCSA ENCODING
# ============================================================

ALPHA: int = 5

BASE_TO_DIGIT: dict[str, int] = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "-": 0,
}

DIGIT_TO_BASE: dict[int, str] = {value: key for key, value in BASE_TO_DIGIT.items()}

# ============================================================
# 6. SCORING CONSTANTS
# ============================================================

MATCH_SCORE: int = 2
MISMATCH_SCORE: int = -1
GAP_SCORE: int = -2

# ============================================================
# 7. CANDIDATE WINDOW GENERATION
# ============================================================

WINDOW_STRIDE: int = 5

# ============================================================
# 8. INDEL / ERROR MODEL SETTINGS
# ============================================================

DEFAULT_INDEL_PROB: float = 0.01
DEFAULT_MAX_INDEL_LEN: int = 2

# ============================================================
# 9. WCSA LOOKUP TABLE
# ============================================================
# Maximum possible S value:
# max_ref_digit + ALPHA * max_read_digit = 4 + 5*4 = 24
# Therefore LUT size = 25 to index from 0 to 24.

WCSA_LUT_SIZE: int = 25
WCSA_SCORE_LUT: np.ndarray = np.zeros(WCSA_LUT_SIZE, dtype=np.int8)

# Specific mismatch penalties from the notebook logic
SPECIFIC_MISMATCH_PENALTIES: dict[int, int] = {
    7: -1,   # C vs A
    8: -2,   # G vs A
    9: -3,   # T vs A
    11: -1,  # A vs C
    13: -1,  # G vs C
    14: -2,  # T vs C
    16: -2,  # A vs G
    17: -1,  # C vs G
    19: -1,  # T vs G
    21: -3,  # A vs T
    22: -2,  # C vs T
    23: -1,  # G vs T
}


def build_wcsa_score_lut() -> np.ndarray:
    """
    Build a lookup table for WCSA score decoding based on the encoded
    column-sum value S = ref + ALPHA * read.

    Match positions are assigned MATCH_SCORE.
    Mismatch positions use SPECIFIC_MISMATCH_PENALTIES when available,
    otherwise fall back to the default MISMATCH_SCORE.
    """
    lut = np.zeros(WCSA_LUT_SIZE, dtype=np.int8)

    for ref_base, ref_digit in BASE_TO_DIGIT.items():
        for read_base, read_digit in BASE_TO_DIGIT.items():
            s_value = ref_digit + ALPHA * read_digit

            if s_value >= WCSA_LUT_SIZE:
                continue

            if ref_base == "-" or read_base == "-":
                lut[s_value] = GAP_SCORE
            elif ref_base == read_base:
                lut[s_value] = MATCH_SCORE
            else:
                lut[s_value] = SPECIFIC_MISMATCH_PENALTIES.get(s_value, MISMATCH_SCORE)

    return lut


# Initialize LUT once at import time
WCSA_SCORE_LUT = build_wcsa_score_lut()

# ============================================================
# 10. PLOTTING / OUTPUT DEFAULTS
# ============================================================

FIG_DPI: int = 300
HIST_BINS: int = 30
RESULTS_FIG_DIR: str = "results/figures"
RESULTS_TABLE_DIR: str = "results/tables"
RESULTS_OUTPUT_DIR: str = "results/mapped_output"
