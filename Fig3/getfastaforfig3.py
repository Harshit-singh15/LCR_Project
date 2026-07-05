from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from Fig2.getfasta import extract_fasta_from_beds
import sys

extract_fasta_from_beds(
    bed_dir=sys.argv[1],
    fasta_file= sys.argv[2],
    output_dir=sys.argv[3],
)