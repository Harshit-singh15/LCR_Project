#!/bin/bash

JAR=~/Lcr_Project/tools/treks-hpc/target/T-ReksHPC_0.1-SNAPSHOT.jar

OUTDIR=treks_single
mkdir -p "$OUTDIR"

touch "$OUTDIR/failed_sequences.txt"

HEADER_WRITTEN=0

for fasta in single_seqs/*.fasta
do
    base=$(basename "$fasta" .fasta)

    echo "Processing $base"

    timeout 120 \
    java -Xmx1G \
    -jar "$JAR" \
    -f "$fasta" \
    -t "$OUTDIR/${base}.tsv" \
    -a "$OUTDIR/${base}.aln" \
    -c /usr/bin/clustalw \
    > "$OUTDIR/${base}.log" 2>&1

    status=$?

    if [ $status -eq 124 ]; then
        echo "$base    TIMEOUT" >> "$OUTDIR/failed_sequences.txt"
        continue
    fi

    if [ $status -ne 0 ]; then
        echo "$base    ERROR" >> "$OUTDIR/failed_sequences.txt"
        continue
    fi

    if [ -s "$OUTDIR/${base}.tsv" ]; then

        if [ $HEADER_WRITTEN -eq 0 ]; then
            cat "$OUTDIR/${base}.tsv" > "$OUTDIR/treks_combined.tsv"
            HEADER_WRITTEN=1
        else
            tail -n +2 "$OUTDIR/${base}.tsv" >> "$OUTDIR/treks_combined.tsv"
        fi

    fi

done
