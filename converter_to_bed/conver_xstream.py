from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config

print("[INFO] Running conver_xstream")

input_folder = config.LCRBYTOOLS_DIR
output_folder = config.BED_DIR
output_folder.mkdir(parents=True, exist_ok=True)

for input_file in sorted(input_folder.iterdir()):
    if not input_file.is_file() or input_file.suffix != ".html":
        continue

    if "xstream" not in input_file.stem.lower():
        continue

    print(f"[INFO] Processing {input_file.name}")

    with open(
        input_file,
        "r",
        encoding="utf-8",
        errors="ignore"
    ) as f:
        html = f.read()

    protein_pattern = re.compile(
        r'(?:tr|sp)\|([^|]+)\|[^<]+'
    )

    position_pattern = re.compile(
        r'>(\d+)-(\d+)<'
    )

    lines = []

    current_protein = None

    for block in re.split(
        r'<A NAME="\d+"></A>',
        html
    ):

        prot = protein_pattern.search(block)

        if prot:
            current_protein = prot.group(1)

        if current_protein:

            for start, end in position_pattern.findall(block):

                lines.append(
                    (
                        current_protein,
                        int(start),
                        int(end)
                    )
                )

    output_bed = output_folder / f"{input_file.stem}.bed"

    with output_bed.open("w", encoding="utf-8") as out:

        for prot, start, end in lines:

            out.write(
                f"{prot}\t{start}\t{end}\n"
            )

    print(f"[SUCCESS] {input_file.name} converted")
    print(f"Extracted {len(lines)} regions")
    print(f"Output: {output_bed}")

print("[SUCCESS] Converter completed")