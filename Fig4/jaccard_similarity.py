from pathlib import Path
import subprocess
import pandas as pd
import time
import sys
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
# Compute Jaccard Matrix
# ==========================================================



def run_jaccard_matrix(input_dir, output_dir):
    overall_start = time.perf_counter()

    print("\n" + "=" * 60)
    print("Figure 4 : BEDTools Jaccard Matrix")
    print("=" * 60)

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = output_dir 

    bed_files = sorted(input_dir.glob("*.bed"))

    if len(bed_files) < 2:
        raise RuntimeError(
            "At least two BED files are required."
        )

    names = [f.stem for f in bed_files]
    n = len(bed_files)

    matrix = pd.DataFrame(
        index=names,
        columns=names,
        dtype=float
    )

    total_pairs = n * (n + 1) // 2
    current = 0

    for i in range(n):
        matrix.iloc[i, i] = 1.0

        for j in range(i + 1, n):
            current += 1
            print(
                f"[{current}/{total_pairs}] "
                f"{names[i]}  vs  {names[j]}"
            )

            # FIX: Only format paths for WSL if we're on Windows
            if os.name == 'nt':
                file1 = str(windows_to_wsl(bed_files[i]))
                file2 = str(windows_to_wsl(bed_files[j]))
                base_cmd = ["wsl", "bedtools", "jaccard"]
            else:
                file1 = str(bed_files[i].resolve())
                file2 = str(bed_files[j].resolve())
                base_cmd = ["bedtools", "jaccard"]

            result = subprocess.run(
                base_cmd + ["-a", file1, "-b", file2],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(result.stderr)
                value = 0
            else:
                lines = result.stdout.strip().splitlines()
                if len(lines) < 2:
                    value = 0
                else:
                    cols = lines[-1].split()
                    value = float(cols[2])

            matrix.iloc[i, j] = value
            matrix.iloc[j, i] = value

    matrix.index.name = "Method"

    matrix.to_csv(
        output_file,
        sep="\t",
        float_format="%.6f"
    )

    elapsed = time.perf_counter() - overall_start

    print("\n" + "=" * 60)
    print("Finished")
    print("=" * 60)
    print(f"Methods      : {n}")
    print(f"Comparisons  : {n*(n-1)//2}")
    print(f"Saved        : {output_file}")
    print(f"Time         : {elapsed:.2f} sec")
    print("=" * 60)

    return output_file


# ==========================================================
# Test
# ==========================================================

if __name__ == "__main__":

    run_jaccard_matrix(sys.argv[1], sys.argv[2])