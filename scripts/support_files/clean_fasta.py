from Bio import SeqIO

input_fasta = "outputs\\mouse_proteome.fasta"
output_fasta = "Fig2\\mouse_proteome_clean.fasta"

with open(output_fasta, "w") as out:

    for record in SeqIO.parse(input_fasta, "fasta"):

        accession = record.id.split("|")[1]

        record.id = accession
        record.name = accession
        record.description = accession

        SeqIO.write(record, out, "fasta")

print("Done")