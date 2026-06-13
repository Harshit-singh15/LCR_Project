from Bio import SeqIO
import csv
import sys

input_fasta = sys.argv[1]
output_bed = sys.argv[2]

with open(output_bed, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Protein_ID", "Start", "End"])

    for record in SeqIO.parse(input_fasta, "fasta"):

        seq = str(record.seq)
        protein_id = record.id

        in_lcr = False
        start = None

        for i, aa in enumerate(seq, start=1):

            if aa.islower():

                if not in_lcr:
                    start = i
                    in_lcr = True

            else:

                if in_lcr:
                    writer.writerow([protein_id, start, i - 1])
                    in_lcr = False

        if in_lcr:
            writer.writerow([protein_id, start, len(seq)])