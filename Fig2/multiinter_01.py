from pathlib import Path
import subprocess
import sys
import time


# ==========================================================
# Windows -> WSL
# ==========================================================

def windows_to_wsl(path: Path):

    path = path.resolve()

    drive = path.drive[0].lower()

    rest = path.as_posix().split(":", 1)[1]

    return f"/mnt/{drive}{rest}"


# ==========================================================
# Figure 2
# BEDTools Multiinter
# ==========================================================

def run_multiinter(input_dir, output_dir):

    overall_start = time.perf_counter()

    print("\n" + "=" * 60)
    print("Figure 2 : BEDTools Multiinter")
    print("=" * 60)

    project = Path.cwd()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = output_dir 

    bed_files = sorted(
        input_dir.glob("*.bed")
    )

    if len(bed_files) == 0:
        raise RuntimeError(
            "No BED files found."
        )

    print(f"\nFound {len(bed_files)} BED files.\n")

    command = [
        "wsl",
        "bedtools",
        "multiinter",
        "-i"
    ]

    for bed in bed_files:

        if bed.stat().st_size == 0:

            print(f"Skipping empty file : {bed.name}")

            continue

        command.append(
            windows_to_wsl(bed)
        )

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:

        print(result.stderr)

        raise RuntimeError(
            "BEDTools multiinter failed."
        )

    output_file.write_text(result.stdout)

    elapsed = time.perf_counter() - overall_start

    print("Finished.\n")

    print(f"Output : {output_file}")

    print(f"Time   : {elapsed:.2f} sec")

    print("=" * 60)

    return output_file


# ==========================================================
# Test
# ==========================================================

if __name__ == "__main__":

    run_multiinter(sys.argv[1], sys.argv[2])