from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config

SOURCE_DIR = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path("ecoli")

print("[INFO] Preparing analysis workspace...")

source_fasta = None
for fasta_path in SOURCE_DIR.glob("*.fasta"):
    if source_fasta is None:
        source_fasta = fasta_path
    else:
        raise ValueError(f"Expected one FASTA file in {SOURCE_DIR}, found multiple")

if source_fasta is None:
    raise FileNotFoundError(f"No FASTA file found in {SOURCE_DIR}")

config.PROJECT_DIR.mkdir(parents=True, exist_ok=True)

print("[INFO] Copying FASTA...")
shutil.copy2(source_fasta, config.FASTA_FILE)

lcrbytools_source = SOURCE_DIR / "lcrbytools"
analysis_lcrbytools = config.LCRBYTOOLS_DIR

if analysis_lcrbytools.exists():
    shutil.rmtree(analysis_lcrbytools)

print("[INFO] Copying tool outputs...")
shutil.copytree(lcrbytools_source, analysis_lcrbytools)

print("[SUCCESS] Analysis workspace ready.")
