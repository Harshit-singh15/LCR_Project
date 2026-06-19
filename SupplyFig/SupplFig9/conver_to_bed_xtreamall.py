import re
from pathlib import Path

# ==========================
input_dir = Path(r"SupplyFig\SupplFig9\xtreamwithmulti_m")
output_dir = Path(r"SupplyFig\SupplFig9")
# ==========================

output_dir.mkdir(parents=True, exist_ok=True)

# Protein IDs
protein_pattern = re.compile(
    r'(?:tr|sp)\|[^|]+\|[^<]+'
)

# Position ranges
position_pattern = re.compile(
    r'>(\d+)-(\d+)<'
)

html_files = list(input_dir.glob("*.html"))

for input_html in html_files:

    with open(input_html, "r", encoding="utf-8", errors="ignore") as f:
        html = f.read()

    lines = []

    current_protein = None

    for block in re.split(r'<A NAME="\d+"></A>', html):

        prot = protein_pattern.search(block)

        if prot:
            current_protein = prot.group()

        if current_protein:

            positions = position_pattern.findall(block)

            if positions:

                start, end = positions[0]

                lines.append(
                    (current_protein, start, end)
                )

    output_bed = output_dir / f"{input_html.stem}.bed"

    with open(output_bed, "w") as out:

        out.write("Protein_ID\tStart\tEnd\n")

        for prot, start, end in lines:

            out.write(
                f"{prot}\t{start}\t{end}\n"
            )

    print(
        f"{input_html.name} -> {output_bed.name} "
        f"({len(lines)} regions)"
    )

print(f"\nProcessed {len(html_files)} files.")