# inspect_reference_bed.py

bed_file = "bed_files\\use_reference_lcr.bed"

rows = 0
proteins = set()
lcr_ids = set()

with open(bed_file) as f:
    next(f)

    for line in f:
        cols = line.rstrip().split("\t")

        rows += 1
        proteins.add(cols[0])

        if len(cols) >= 4:
            lcr_ids.add(cols[3])

print("Rows:", rows)
print("Unique proteins:", len(proteins))
print("Unique LCR_IDs:", len(lcr_ids))