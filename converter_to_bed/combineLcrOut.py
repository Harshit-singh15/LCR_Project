import os
import re

# load mapping
mapping = {}

with open("C:\\Users\\91892\\Documents\\MATLAB\\LcrFinder\\combined\\header_map_arabidopsis.tsv") as f:
    for line in f:
        idx, protein = line.rstrip().split("\t",1)
        mapping[idx] = protein

out = open("C:\\Users\\91892\\Documents\\MATLAB\\LcrFinder\\combined\\LCRFinder_combined_arabidopsis.tsv","w")

for fname in sorted(os.listdir("C:\\Users\\91892\\Documents\\MATLAB\\LcrFinder\\outputLcr")):

    if not fname.endswith("_"):
        continue

    idx = fname.rstrip("_")

    if idx not in mapping:
        continue

    protein = mapping[idx]

    with open(os.path.join("C:\\Users\\91892\\Documents\\MATLAB\\LcrFinder\\outputLcr",fname)) as f:

        for line in f:

            if line.startswith(">"):
                continue

            m = re.search(r'\((\d+)\).*?\((\d+)\)', line)

            if m:
                start = m.group(1)
                end = m.group(2)

                out.write(
                    f"{protein}\t{start}\t{end}\n"
                )

out.close()