import re

# ==========================
input_html = r"outputs\output_m5XSTREAM__i0.7_g3_m5_e2.0_out_2.html"
output_bed = r"SupplyFig\SupplFig9\xstream_m1.bed"
# ==========================

with open(input_html, "r", encoding="utf-8", errors="ignore") as f:
    html = f.read()

# Protein IDs
protein_pattern = re.compile(
    r'(?:tr|sp)\|[^|]+\|[^<]+'
)

# Position ranges
position_pattern = re.compile(
    r'>(\d+)-(\d+)<'
)

proteins = protein_pattern.findall(html)

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

with open(output_bed, "w") as out:

    out.write("Protein_ID\tStart\tEnd\n")

    for prot, start, end in lines:

        out.write(
            f"{prot}\t{start}\t{end}\n"
        )

print(f"Extracted {len(lines)} regions")