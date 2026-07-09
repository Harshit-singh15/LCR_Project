from pathlib import Path
import subprocess
import sys
import time
import os


# ==========================================================
# Windows -> WSL
# ==========================================================

def windows_to_wsl(path: Path):
    path = path.resolve()
    
    # FIX: If there's no drive letter, we are already on Linux/Render
    if not path.drive:
        return path.as_posix()

    # This part runs safely if a Windows drive letter (like C:) is found
    drive = path.drive[0].lower()
    rest = path.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


# ==========================================================
# Figure 2
# BEDTools Multiinter
# ==========================================================

import os
import subprocess
import time
from pathlib import Path
import pandas as pd

def run_multiinter(input_dir, output_dir):
    overall_start = time.perf_counter()

    print("\n" + "=" * 60)
    print("Figure 2 : BEDTools Multiinter")
    print("=" * 60)

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

    # FIX: Only include "wsl" if running on Windows
    if os.name == 'nt':
        command = ["wsl", "bedtools", "multiinter", "-i"]
    else:
        command = ["bedtools", "multiinter", "-i"]

    for bed in bed_files:
        if bed.stat().st_size == 0:
            print(f"Skipping empty file : {bed.name}")
            continue

        # Convert path or leave it native depending on OS
        command.append(
            str(windows_to_wsl(bed)) if os.name == 'nt' else str(bed.resolve())
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