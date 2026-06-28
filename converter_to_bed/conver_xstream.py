import os
import re

# ==========================
input_html = r"zebrafish\lcrbytools_zebrafish\xstream_m1_zebrafish.html"
output_bed = r"zebrafish\bed_zebrafish\xstream_m1_zebrafish.bed"
# ==========================

os.makedirs(
    os.path.dirname(output_bed),
    exist_ok=True
)

with open(
    input_html,
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

with open(output_bed, "w") as out:

    for prot, start, end in lines:

        out.write(
            f"{prot}\t{start}\t{end}\n"
        )

print(f"Extracted {len(lines)} regions")
print(f"Output: {output_bed}")