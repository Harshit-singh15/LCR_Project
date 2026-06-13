from pathlib import Path
import pandas as pd

INPUT_DIR = Path("outputs_prerequisite\\LCR_Count")

print("Directory exists:", INPUT_DIR.exists())
print("Files found:")

for f in INPUT_DIR.glob("*"):
    print(f)

TOTAL_PROTEINS = 21853


print("\n=== COUNT QC ===\n")

for file in INPUT_DIR.glob("*.tsv"):

    df = pd.read_csv(file, sep="\t")

    total = df["Count"].sum()

    status = "PASS"

    if total != TOTAL_PROTEINS:
        status = "FAIL"

    print(
        f"{file.stem:40s} "
        f"{total:6d} "
        f"{status}"
    )