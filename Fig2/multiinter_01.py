"""
Figure 2 - Step 1
Generate consensus LCR regions using BEDTools multiinter.

Equivalent command:

bedtools multiinter \
-i bed_bedtools/*.bed \
> dataforFig2/multiinter.tsv
"""

from pathlib import Path
import subprocess
import time


# ==========================================================
# Helper
# ==========================================================

def windows_to_wsl(path: Path) -> str:
    """
    Convert a Windows path to a WSL path.

    Example
    -------
    D:\\Project\\bed.bed

    becomes

    /mnt/d/Project/bed.bed
    """

    path = path.resolve()

    drive = path.drive[0].lower()

    rest = path.as_posix().split(":", 1)[1]

    return f"/mnt/{drive}{rest}"


# ==========================================================
# BEDTools Multiinter
# ==========================================================

def run_multiinter():

    start = time.time()

    print("\n" + "=" * 60)
    print("Figure 2 : BEDTools multiinter")
    print("=" * 60)

    project_dir = Path.cwd()

    input_dir = project_dir / "ecoli" / "bed_bedtools_Ecoli"

    output_dir = project_dir / "ecoli" / "dataforFig2"

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = output_dir / "multiinter.tsv"

    if not input_dir.exists():
        raise FileNotFoundError(
            f"Input folder not found:\n{input_dir}"
        )

    bed_files = sorted(input_dir.glob("*.bed"))

    if len(bed_files) == 0:
        raise RuntimeError(
            "No BED files found inside bed_bedtools_Ecoli."
        )

    print(f"\nFound {len(bed_files)} BED files\n")

    wsl_beds = [
        windows_to_wsl(f)
        for f in bed_files
        if f.stat().st_size > 0
    ]

    if len(wsl_beds) == 0:
        raise RuntimeError(
            "All BED files are empty."
        )

    output_wsl = windows_to_wsl(output_file)

    command = (
        f"bedtools multiinter "
        f"-i {' '.join(wsl_beds)} "
        f"> {output_wsl}"
    )

    print("Running BEDTools...\n")

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

        raise RuntimeError(
            "BEDTools multiinter failed."
        )

    elapsed = time.time() - start

    print("Done.\n")

    print(f"Output : {output_file}")

    print(f"Time   : {elapsed:.2f} sec")

    print("=" * 60)

    return output_file


# ==========================================================
# Test
# ==========================================================

if __name__ == "__main__":

    run_multiinter()