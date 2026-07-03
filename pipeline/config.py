from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent / "analysis"

FASTA_FILE = PROJECT_DIR / "organism.fasta"

LCRBYTOOLS_DIR = PROJECT_DIR / "lcrbytools"

BED_DIR = PROJECT_DIR / "bed"

DATA_DIR = PROJECT_DIR / "data"

FIG1_DATA_DIR = DATA_DIR / "Fig1"

FIG2_DATA_DIR = DATA_DIR / "Fig2"

EXTRACTED_SEQUENCES_DIR = FIG1_DATA_DIR / "extracted_sequences"

PROTEIN_LENGTHS_FILE = FIG1_DATA_DIR / "protein_lengths.tsv"

BED_GLOB_PATTERN = "*.bed"

FIGURES_DIR = PROJECT_DIR / "figures"

FIG1_FIGURES_DIR = FIGURES_DIR / "Fig1"

FIG2_FIGURES_DIR = FIGURES_DIR / "Fig2"

REPORT_DIR = PROJECT_DIR / "report"

LOG_DIR = PROJECT_DIR / "logs"