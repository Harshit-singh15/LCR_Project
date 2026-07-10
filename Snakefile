"""
Snakemake workflow: LCR benchmarking pipeline, BED files -> Fig 7.
Fig 8 is intentionally excluded (depends on external DisProt/SIFTS data — separate workflow).

Run:
    snakemake -n                       # dry run — shows the plan without executing
    snakemake --dag | dot -Tpng > dag.png   # visualize the DAG (requires graphviz)
    snakemake --cores 4                # actually run, using up to 4 cores in parallel
    snakemake --cores 4 fig7           # run only what's needed to build the fig7 target
"""
from pathlib import Path

PROJECT_ROOT = workflow.basedir
configfile: "config.yaml"

ORG = config["organism"]
BASE = config["base_dir"]
S = config["scripts_dir"]
COMBINED = config.get("combined_figures", True)

# ---------------------------------------------------------------------------
# Core inputs / shared preprocessing outputs
# ---------------------------------------------------------------------------
FASTA          = f"{BASE}/{ORG}.fasta"
BED_DIR        = f"{BASE}/bed_{ORG}"                 # you said this already exists
CLEANED_FASTA  = f"{BASE}/{ORG}_cleaned.fasta"
BED_BT_DIR     = f"{BASE}/bed_bedtools_{ORG}"


# ---------------------------------------------------------------------------
# Final target — everything Fig1-Fig7 needs to produce
# ---------------------------------------------------------------------------
rule all:
    input:
        f"{BASE}/Fig_outputs/Fig1/Figure1A_Length_Heatmap.png",
        f"{BASE}/Fig_outputs/Fig1/Figure1B_Coverage_Heatmap.png",
        f"{BASE}/Fig_outputs/Fig1/Figure1C_CountDistribution.png",
        f"{BASE}/Fig_outputs/Fig1/Figure1D_AminoAcidComposition.png",
        f"{BASE}/Fig_outputs/Fig1/Figure1E_Entropy_Boxplot.png",
        (f"{BASE}/Fig_outputs/Fig1/Figure1_Combined.png" if COMBINED else []),
        f"{BASE}/Fig_outputs/Fig2/Fig2A_PeptideMotifs.png",
        f"{BASE}/Fig_outputs/Fig2/Fig2B_Entropy.png",
        f"{BASE}/Fig_outputs/Fig2/Fig2C_Purity.png",
        f"{BASE}/Fig_outputs/Fig3/Fig3_Purity.png",
        f"{BASE}/Fig_outputs/Fig4/Fig4_Jaccard.png",
        (f"{BASE}/Fig_outputs/Fig5/Fig5_heatmaps_done.flag"),
        (f"{BASE}/Fig_outputs/Fig5/Fig5_Combined.png" if COMBINED else []),
        f"{BASE}/Fig_outputs/Fig6/Fig6A_Genelength.png",
        f"{BASE}/Fig_outputs/Fig6/Fig6B_LCRCount.png",
        f"{BASE}/Fig_outputs/Fig6/Fig6C_Coverage.png",
        f"{BASE}/Fig_outputs/Fig6/Fig6D_EntropyRatio.png",
        (f"{BASE}/Fig_outputs/Fig6/Fig6_Combined.png" if COMBINED else []),
        f"{BASE}/Fig_outputs/Fig7/Fig7_summary.png",
        f"{BASE}/reports/Benchmark_Report.pdf"


# ---------------------------------------------------------------------------
# Shared preprocessing (used by Fig2, Fig3, Fig4, Fig5, Fig6)
# ---------------------------------------------------------------------------

# 01_cleanfasta_headers.py
# argv order: fasta_in, fasta_out
rule clean_fasta:
    input:
        fasta=FASTA
    output:
        cleaned=CLEANED_FASTA
    shell:
        "python {S}/Fig2/01_cleanfasta_headers.py {input.fasta} {output.cleaned}"

# 00_prepare_bedfiles.py
# argv order: bed_in_dir, bed_out_dir
rule prepare_bedtools_beds:
    input:
        bed_dir=BED_DIR
    output:
        directory(BED_BT_DIR)
    shell:
        "python {S}/Fig2/00_prepare_bedfiles.py {input.bed_dir} {output}"


# ---------------------------------------------------------------------------
# Fig 1 — length, coverage, count, AA composition, entropy
# ---------------------------------------------------------------------------

# Run_pipeline_1.py
# argv order: fasta, bed_dir, output_base_dir (dataforFig1/)
rule fig1_extract_sequences:
    input:
        fasta=FASTA,
        bed_dir=BED_DIR
    output:
        extracted=directory(f"{BASE}/dataforFig1/extracted_sequences"),
        lengths=f"{BASE}/dataforFig1/protein_lengths.tsv"
    shell:
        "python {S}/Fig1/run_pipeline_1.py {input.fasta} {input.bed_dir} {BASE}/dataforFig1"

# gen_len_counts_2.py
# argv order: extracted_sequences_dir, length_out_dir, count_out_dir
# gen_len_counts_2.py
# argv order:
# extracted_sequences_dir, protein_lengths.tsv, length_out_dir, count_out_dir

rule fig1_length_counts:
    input:
        extracted = f"{BASE}/dataforFig1/extracted_sequences",
        protein_lengths = f"{BASE}/dataforFig1/protein_lengths.tsv"
    output:
        length_dir = directory(f"{BASE}/dataforFig1/LCR_Length"),
        count_dir = directory(f"{BASE}/dataforFig1/LCR_Count")
    shell:
        """
        python {S}/Fig1/gen_len_counts_2.py \
            {input.extracted} \
            {input.protein_lengths} \
            {output.length_dir} \
            {output.count_dir}
        """
# Gen_diversity_3.py
# argv order: extracted_sequences_dir, entropy_out_dir
rule fig1_entropy:
    input:
        extracted=f"{BASE}/dataforFig1/extracted_sequences"
    output:
        directory(f"{BASE}/dataforFig1/ShanonEntropy")
    shell:
        "python {S}/Fig1/gen_diversity_3.py {input.extracted} {output}"

# Gen_coverage_4.py
# argv order: extracted_sequences_dir, protein_lengths_tsv, out_dir, fasta
rule fig1_coverage:
    input:
        extracted=f"{BASE}/dataforFig1/extracted_sequences",
        lengths=f"{BASE}/dataforFig1/protein_lengths.tsv",
        fasta=FASTA
    output:
        directory(f"{BASE}/dataforFig1/LCR_Coverage")
    shell:
        "python {S}/Fig1/Gen_coverage_4.py {input.extracted} {input.lengths} {output} {input.fasta}"

# Aa_composition_5.py
# argv order: extracted_sequences_dir, out_dir
rule fig1_aa_composition:
    input:
        extracted=f"{BASE}/dataforFig1/extracted_sequences"
    output:
        directory(f"{BASE}/dataforFig1/Amino_acid")
    shell:
        "python {S}/Fig1/aa_composition_5.py {input.extracted} {output}"

# fig1A_length_hm.py  — argv order: length_dir, out_dir
rule fig1A_plot:
    input:
        f"{BASE}/dataforFig1/LCR_Length"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1A_Length_Heatmap.png"
    shell:
        "python {S}/Fig1/fig1A_length_hm.py {input} {BASE}/Fig_outputs/Fig1"

# fig1B_coverage_hm.py — argv order: coverage_dir, out_dir
rule fig1B_plot:
    input:
        f"{BASE}/dataforFig1/LCR_Coverage"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1B_Coverage_Heatmap.png"
    shell:
        "python {S}/Fig1/fig1B_coverage_hm.py {input} {BASE}/Fig_outputs/Fig1"

# fig1C_count_dist.py — argv order: count_dir, out_dir
rule fig1C_plot:
    input:
        f"{BASE}/dataforFig1/LCR_Count"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1C_CountDistribution.png"
    shell:
        "python {S}/Fig1/fig1C_count_dist.py {input} {BASE}/Fig_outputs/Fig1"

# fig1D_aa_comp.py — argv order: amino_acid_dir, out_dir
rule fig1D_plot:
    input:
        f"{BASE}/dataforFig1/Amino_acid"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1D_AminoAcidComposition.png"
    shell:
        "python {S}/Fig1/fig1D_aa_comp.py {input} {BASE}/Fig_outputs/Fig1"

# fig1E_entropybox.py — argv order: entropy_dir, out_dir
rule fig1E_plot:
    input:
        f"{BASE}/dataforFig1/ShanonEntropy"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1E_Entropy_Boxplot.png"
    shell:
        "python {S}/Fig1/fig1E_entropybox.py {input} {BASE}/Fig_outputs/Fig1"

# Comb_fig1.py — argv order: A, B, C, D, E, out_path  (optional)
rule fig1_combined:
    input:
        a=f"{BASE}/Fig_outputs/Fig1/Figure1A_Length_Heatmap.png",
        b=f"{BASE}/Fig_outputs/Fig1/Figure1B_Coverage_Heatmap.png",
        c=f"{BASE}/Fig_outputs/Fig1/Figure1C_CountDistribution.png",
        d=f"{BASE}/Fig_outputs/Fig1/Figure1D_AminoAcidComposition.png",
        e=f"{BASE}/Fig_outputs/Fig1/Figure1E_Entropy_Boxplot.png"
    output:
        f"{BASE}/Fig_outputs/Fig1/Figure1_Combined.png"
    shell:
        "python {S}/Fig1/comb_fig1.py {input.a} {input.b} {input.c} {input.d} {input.e} {output}"


# ---------------------------------------------------------------------------
# Fig 2 — consensus regions & motifs
# ---------------------------------------------------------------------------

# multiinter_01.py (your python port of `bedtools multiinter`)
# argv order: bed_bedtools_dir, out_tsv
rule fig2_multiinter:
    input:
        BED_BT_DIR
    output:
        f"{BASE}/dataforFig2/multiinter.tsv"
    shell:
        "python {S}/Fig2/multiinter_01.py {input} {output}"

# 02_consensus_beds.py — argv order: multiinter_tsv, out_dir
rule fig2_consensus_beds:
    input:
        f"{BASE}/dataforFig2/multiinter.tsv"
    output:
        directory(f"{BASE}/dataforFig2/02_consensus_beds")
    shell:
        "python {S}/Fig2/02_consensus_beds.py {input} {output}"

# 03_getfasta.py (your python port of `bedtools getfasta` + extract_sequences.sh)
# argv order: cleaned_fasta, consensus_beds_dir, out_dir
rule fig2_consensus_fastas:
    input:
        beds=f"{BASE}/dataforFig2/02_consensus_beds",
        fasta=CLEANED_FASTA
    output:
        directory(f"{BASE}/dataforFig2/03_consensus_fastas")
    shell:
        "python {S}/Fig2/getfasta.py {input.beds} {input.fasta} {output}"

# 04_substring_motif.py — argv order: consensus_fastas_dir, out_dir
rule fig2_substring_motifs:
    input:
        f"{BASE}/dataforFig2/03_consensus_fastas"
    output:
        directory(f"{BASE}/dataforFig2/04_substring_motifs")
    shell:
        "python {S}/Fig2/04_substring_motif.py {input} {output}"

# 05_subfilter.py — argv order: substring_motifs_dir, out_dir
rule fig2_filtered:
    input:
        f"{BASE}/dataforFig2/04_substring_motifs"
    output:
        directory(f"{BASE}/dataforFig2/05_filtered")
    shell:
        "python {S}/Fig2/05_subfilter.py {input} {output}"

# 06_peptide_counts.py — argv order: filtered_dir, out_dir
rule fig2_peptide_counts:
    input:
        f"{BASE}/dataforFig2/05_filtered"
    output:
        directory(f"{BASE}/dataforFig2/06_peptide_counts")
    shell:
        "python {S}/Fig2/06_peptide_counts.py {input} {output}"

# entropy.py — argv order: consensus_fastas_dir, out_dir
rule fig2_entropy:
    input:
        f"{BASE}/dataforFig2/03_consensus_fastas"
    output:
        f"{BASE}/dataforFig2/entropy.tsv"
    shell:
        "python {S}/Fig2/entropy.py {input} {output}"

# purity.py — argv order: consensus_fastas_dir, out_dir
rule fig2_purity:
    input:
        f"{BASE}/dataforFig2/03_consensus_fastas"
    output:
        f"{BASE}/dataforFig2/purity.tsv"
    shell:
        "python {S}/Fig2/purity.py {input} {output}"

# fig2a.py — argv order: peptide_counts_dir, out_dir
rule fig2a_plot:
    input:
        f"{BASE}/dataforFig2/06_peptide_counts"
    output:
        f"{BASE}/Fig_outputs/Fig2/Fig2A_PeptideMotifs.png"
    shell:
        "python {S}/Fig2/fig2a.py {input} {BASE}/Fig_outputs/Fig2"

# fig2b.py — argv order: entropy_dir, out_dir
rule fig2b_plot:
    input:
        f"{BASE}/dataforFig2/entropy.tsv"
    output:
        f"{BASE}/Fig_outputs/Fig2/Fig2B_Entropy.png"
    shell:
        "python {S}/Fig2/fig2b.py {input} {BASE}/Fig_outputs/Fig2"

# fig2c.py — argv order: purity_dir, out_dir
rule fig2c_plot:
    input:
        f"{BASE}/dataforFig2/purity.tsv"
    output:
        f"{BASE}/Fig_outputs/Fig2/Fig2C_Purity.png"
    shell:
        "python {S}/Fig2/fig2c.py {input} {BASE}/Fig_outputs/Fig2"


# ---------------------------------------------------------------------------
# Fig 3 — purity across tools
# ---------------------------------------------------------------------------

# getfastaforfig3.py (your python port of makefastas.sh)
# argv order: bed_bedtools_dir, cleaned_fasta, out_dir
rule fig3_fastas:
    input:
        beds=BED_BT_DIR,
        fasta=CLEANED_FASTA
    output:
        directory(f"{BASE}/dataforFig3/01_fastas")
    shell:
        "python {S}/Fig3/getfastaforfig3.py {input.beds} {input.fasta} {output}"

# Purity_calculation.py — argv order: fastas_dir, out_dir
rule fig3_purity_metrics:
    input:
        f"{BASE}/dataforFig3/01_fastas"
    output:
        directory(f"{BASE}/dataforFig3/02_metrics")
    shell:
        "python {S}/Fig3/purity_calculation.py {input} {output}"

# fig3_plot.py — argv order: metrics_dir, out_dir
rule fig3_plot:
    input:
        (f"{BASE}/dataforFig3/02_metrics")
    output:
        f"{BASE}/Fig_outputs/Fig3/Fig3_Purity.png"
    shell:
        "python {S}/Fig3/fig3_plot.py {input} {BASE}/Fig_outputs/Fig3"


# ---------------------------------------------------------------------------
# Fig 4 — pairwise tool overlap (Jaccard)
# ---------------------------------------------------------------------------

# jaccard_similarity.py (your python port of jaccardsimilarity.sh)
# argv order: bed_bedtools_dir, out_tsv
rule fig4_jaccard_matrix:
    input:
        BED_BT_DIR
    output:
        f"{BASE}/dataforFig4/jaccard_matrix.tsv"
    shell:
        "python {S}/Fig4/jaccard_similarity.py {input} {output}"

# Fig4_plot.py — argv order: jaccard_tsv, out_dir
rule fig4_plot:
    input:
        f"{BASE}/dataforFig4/jaccard_matrix.tsv"
    output:
        f"{BASE}/Fig_outputs/Fig4/Fig4_Jaccard.png"
    shell:
        "python {S}/Fig4/fig4_plot.py {input} {BASE}/Fig_outputs/Fig4"


# ---------------------------------------------------------------------------
# Fig 5 — mutation metrics
# ---------------------------------------------------------------------------

# mutation_metrices.py — argv order: fastas_dir (dataforFig3/01_fastas), out_dir
rule fig5_metrics:
    input:
        f"{BASE}/dataforFig3/01_fastas"
    output:
        directory(f"{BASE}/dataforFig5/02_metrics")
    shell:
        "python {S}/Fig5/mutation_metrices.py {input} {output}"

# fig5_hm.py produces one PNG per tool — argv order: metrics_dir, out_dir
# Using a flag file since the exact per-tool filenames vary with your BED inputs.
rule fig5_heatmaps:
    input:
        f"{BASE}/dataforFig5/02_metrics"
    output:
        touch(f"{BASE}/Fig_outputs/Fig5/Fig5_heatmaps_done.flag")
    shell:
        "python {S}/Fig5/fig5_hm.py {input} {BASE}/Fig_outputs/Fig5"

# Fig5_combined.py — argv order: fig5_dir, out_path  (optional)
rule fig5_combined:
    input:
        f"{BASE}/Fig_outputs/Fig5/Fig5_heatmaps_done.flag"
    output:
        f"{BASE}/Fig_outputs/Fig5/Fig5_Combined.png"
    shell:
        "python {S}/Fig5/Fig5_combined.py {BASE}/Fig_outputs/Fig5 {output}"


# ---------------------------------------------------------------------------
# Fig 6 — windowed reference classification & TPR/FPR
# ---------------------------------------------------------------------------

# 01_seq_parser.py — argv order: cleaned_fasta, out_bed
rule fig6_windows:
    input:
        CLEANED_FASTA
    output:
        f"{BASE}/dataforFig6/ecoli_windows.bed"
    shell:
        "python {S}/Fig6/01_seq_parser.py {input} {output}"

# 02_complexity_plot.py — argv order: cleaned_fasta, windows_bed, out_tsv
rule fig6_complexity:
    input:
        fasta=CLEANED_FASTA,
        windows=f"{BASE}/dataforFig6/ecoli_windows.bed"
    output:
        f"{BASE}/dataforFig6/ecoli_windows_out.tsv"
    shell:
        "python {S}/Fig6/02_complexity_plot.py {input.fasta} {input.windows} {output}"

# 03_classification.py — argv order: windows_out_tsv, out_tsv
rule fig6_classification:
    input:
        f"{BASE}/dataforFig6/ecoli_windows_out.tsv"
    output:
        f"{BASE}/dataforFig6/ecoli_windows_classified.tsv"
    shell:
        "python {S}/Fig6/03_classification.py {input} {output}"

# 04_annotate.py — argv order: classified_tsv, out_bed
rule fig6_annotate:
    input:
        f"{BASE}/dataforFig6/ecoli_windows_classified.tsv"
    output:
        f"{BASE}/dataforFig6/ecoli_windows_real.bed"
    shell:
        "python {S}/Fig6/04_annotate.py {input} {output}"

# Protein_confusion.py — argv order: cleaned_fasta, reference_bed, tool_bed_dir, out_tsv
rule fig6_confusion:
    input:
        fasta=CLEANED_FASTA,
        reference=f"{BASE}/dataforFig6/ecoli_windows_real.bed",
        tools=BED_BT_DIR
    output:
        f"{BASE}/dataforFig6/protein_confusion.tsv"
    shell:
        "python {S}/Fig6/protein_confusion.py {input.fasta} {input.reference} {input.tools} {output}"

# references_metrices.py — argv order: reference_bed, cleaned_fasta, out_tsv
rule fig6_reference_metrics:
    input:
        reference=f"{BASE}/dataforFig6/ecoli_windows_real.bed",
        fasta=CLEANED_FASTA
    output:
        f"{BASE}/dataforFig6/reference_metrics.tsv"
    shell:
        "python {S}/Fig6/reference_metrices.py {input.reference} {input.fasta} {output}"

# Tpr_fpr.py — argv order: reference_metrics_tsv, confusion_tsv, out_dir
rule fig6_tpr_fpr:
    input:
        metrics=f"{BASE}/dataforFig6/reference_metrics.tsv",
        confusion=f"{BASE}/dataforFig6/protein_confusion.tsv"
    output:
        directory(f"{BASE}/dataforFig6/plot_tables")
    shell:
        "python {S}/Fig6/tpr_fpr.py {input.metrics} {input.confusion} {output}"

# Plot_fig6.py — argv order: plot_tables_dir, out_dir
rule fig6_plots:
    input:
        f"{BASE}/dataforFig6/plot_tables"
    output:
        a=f"{BASE}/Fig_outputs/Fig6/Fig6A_Genelength.png",
        b=f"{BASE}/Fig_outputs/Fig6/Fig6B_LCRCount.png",
        c=f"{BASE}/Fig_outputs/Fig6/Fig6C_Coverage.png",
        d=f"{BASE}/Fig_outputs/Fig6/Fig6D_EntropyRatio.png"
    shell:
        "python {S}/Fig6/plot_fig6.py {input} {BASE}/Fig_outputs/Fig6"

# Fig6_combined.py — argv order: plot_tables_dir, out_path  (optional)
rule fig6_combined:
    input:
        f"{BASE}/dataforFig6/plot_tables"
    output:
        f"{BASE}/Fig_outputs/Fig6/Fig6_Combined.png"
    shell:
        "python {S}/Fig6/fig6_combined.py {input} {output}"


# ---------------------------------------------------------------------------
# Fig 7 — TPR/FPR summary bar plots
# ---------------------------------------------------------------------------

# buildsum.py — argv order: plot_tables_dir, out_dir
rule fig7_summaries:
    input:
        f"{BASE}/dataforFig6/plot_tables"
    output:
        directory(f"{BASE}/dataforFig7")
    shell:
        "python {S}/Fig7/buildsum.py {input} {output}"

# Plot_fig7.py — argv order: dataforFig7_dir, out_dir
rule fig7_plot:
    input:
        f"{BASE}/dataforFig7"
    output:
        f"{BASE}/Fig_outputs/Fig7/Fig7_summary.png"
    shell:
        "python {S}/Fig7/plot_fig7.py {input} {BASE}/Fig_outputs/Fig7"

###########################################
# Final Benchmark Report
###########################################
rule final_report:
    input:
        fig1 = f"{BASE}/Fig_outputs/Fig1/Figure1_Combined.png",
        fig2 = f"{BASE}/Fig_outputs/Fig2/Fig2A_PeptideMotifs.png",
        fig2b = f"{BASE}/Fig_outputs/Fig2/Fig2B_Entropy.png",
        fig2c = f"{BASE}/Fig_outputs/Fig2/Fig2C_Purity.png",
        fig3 = f"{BASE}/Fig_outputs/Fig3/Fig3_Purity.png",
        fig4 = f"{BASE}/Fig_outputs/Fig4/Fig4_Jaccard.png",
        fig5 = f"{BASE}/Fig_outputs/Fig5/Fig5_Combined.png",
        fig6 = f"{BASE}/Fig_outputs/Fig6/Fig6_Combined.png",
        fig7 = f"{BASE}/Fig_outputs/Fig7/Fig7_summary.png"

    output:
        pdf = f"{BASE}/reports/Benchmark_Report.pdf"
    shell:
        r"""
        python "{S}/reports/generate_report.py" "{BASE}" "{output.pdf}"
        """