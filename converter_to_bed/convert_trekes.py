from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config

print("[INFO] Running convert_trekes")

input_folder = config.LCRBYTOOLS_DIR
output_folder = config.BED_DIR
output_folder.mkdir(parents=True, exist_ok=True)

for input_file in sorted(input_folder.iterdir()):
    if not input_file.is_file() or input_file.suffix != ".tsv":
        continue

    print(f"[INFO] Processing {input_file.name}")

    df = pd.read_csv(input_file, sep="\t")

    # Extract UniProt accession
    def get_accession(seqid):

        parts = str(seqid).split("|")

        if len(parts) >= 3:
            return parts[1]

        return seqid

    df["Protein_ID"] = df["seqid"].apply(get_accession)

    output_bed = output_folder / f"{input_file.stem}.bed"

    df[["Protein_ID", "start", "end"]].rename(
        columns={
            "start": "Start",
            "end": "End"
        }
    ).to_csv(
        output_bed,
        sep="\t",
        index=False,
        header=False
    )

    print(f"[SUCCESS] {input_file.name} converted")
    print("T-REKS BED file written:", output_bed)
    print("Regions:", len(df))

print("[SUCCESS] Converter completed")