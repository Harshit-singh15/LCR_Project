from pathlib import Path
import subprocess
import time
from multiinter_01 import windows_to_wsl 

# ==========================================================
# Figure 2 - Step 2
# Extract FASTA sequences from consensus BED files
# ==========================================================

def run_getfasta(fasta_file):

    """
    Equivalent bash command

    bedtools getfasta \
        -fi fasta_file \
        -bed consensus_x.bed \
        -name \
        -fo consensus_x.fa
    """

    start = time.time()

    print("\n" + "=" * 60)
    print("Figure 2 : BEDTools getfasta")
    print("=" * 60)

    project_dir = Path.cwd()

    fasta_file = Path(fasta_file)

    if not fasta_file.is_absolute():
        fasta_file = project_dir / fasta_file

    if not fasta_file.exists():
        raise FileNotFoundError(
            f"FASTA file not found:\n{fasta_file}"
        )

    bed_dir = project_dir / "ecoli" / "dataforFig2" / "02_consensus_beds"

    out_dir = project_dir / "ecoli" / "dataforFig2" / "03_consensus_fastas"

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    bed_files = sorted(
        bed_dir.glob("consensus_*.bed")
    )

    if len(bed_files) == 0:
        raise RuntimeError(
            "No consensus BED files found."
        )

    fasta_wsl = windows_to_wsl(fasta_file)

    count = 0

    for bed in bed_files:

        if bed.stat().st_size == 0:
            print(f"Skipping empty file : {bed.name}")
            continue

        output = out_dir / f"{bed.stem}.fa"

        print(f"Processing {bed.stem}")

        command = (
            f"bedtools getfasta "
            f"-fi {fasta_wsl} "
            f"-bed {windows_to_wsl(bed)} "
            f"-name "
            f"-fo {windows_to_wsl(output)}"
        )

        result = subprocess.run(
            [
                "wsl",
                "bash",
                "-c",
                command
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:

            print(result.stderr)

            print(f"Skipped : {bed.name}")

            continue

        count += 1

    elapsed = time.time() - start

    print("\n" + "=" * 60)
    print("Finished.")
    print(f"Consensus FASTA files created : {count}")
    print(f"Output folder : {out_dir}")
    print(f"Time : {elapsed:.2f} sec")
    print("=" * 60)

    return out_dir

run_getfasta("ecoli\ecoli_cleaned.fasta")