input_file = r"outputs\flps2_strict.out"
output_bed = r"bed_files\flps2_strict_mouse.bed"

with open(input_file) as inp, open(output_bed, "w") as out:

    out.write("Protein_ID\tStart\tEnd\n")

    for line in inp:

        line = line.strip()

        if not line:
            continue

        fields = line.split()

        protein_id = fields[0]
        start = fields[4]
        end = fields[5]

        out.write(f"{protein_id}\t{start}\t{end}\n")

print("Conversion complete.")