from importlib.resources import path
from pathlib import Path
import re
import sys

print("[INFO] Running convert_xstream")

# Define your exact input file and desired output file paths
input_file = Path(sys.argv[1])  # Get the input file path from command line argument
output_file = Path(sys.argv[2])  # Get the output file path from command line argument

# Ensure the output directory exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Validate that the input file exists before proceeding
if not input_file.is_file():
    print(f"[ERROR] Input file not found: {input_file}")
    exit(1)

print(f"[INFO] Processing {input_file.name}")

# Read the HTML content
with open(input_file, "r", encoding="utf-8", errors="ignore") as f:
    html = f.read()

# Regex patterns
protein_pattern = re.compile(r'(?:tr|sp)\|([^|]+)\|[^<]+')
position_pattern = re.compile(r'>(\d+)-(\d+)<')

lines = []
current_protein = None

# Split and parse the blocks
for block in re.split(r'<A NAME="\d+"></A>', html):
    prot = protein_pattern.search(block)
    if prot:
        current_protein = prot.group(1)

    if current_protein:
        for start, end in position_pattern.findall(block):
            lines.append((current_protein, int(start), int(end)))

# Write directly to the specified output file
with output_file.open("w", encoding="utf-8") as out:
    for prot, start, end in lines:
        out.write(f"{prot}\t{start}\t{end}\n")

print(f"[SUCCESS] {input_file.name} converted")
print(f"Extracted {len(lines)} regions")
print(f"Output saved to: {output_file}")
print("[SUCCESS] Converter completed")