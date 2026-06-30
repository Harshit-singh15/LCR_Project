input_file = "outputs\\flps_strict_again_mouse.txt"
output_bed = "bed_files\\flps_strict_again_mouse.bed"

with open(input_file) as fin, open(output_bed, "w") as fout:
    fout.write("Protein_ID\tStart\tEnd\n")

    for line in fin:
        line = line.strip()

        if not line:
            continue

        cols = line.split()

        if len(cols) < 5:
            continue

        try:
            protein_id = cols[0]
            start = int(cols[3])
            end   = int(cols[4])

            fout.write(f"{protein_id}\t{start}\t{end}\n")

        except ValueError:
            pass

print("Done")