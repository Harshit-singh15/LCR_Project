"""
LCR Benchmarking Tool — Steps 1 & 2: Upload page + background pipeline runner + status page.

Run this with:
    python app.py

Then open http://127.0.0.1:5000 in your browser.
"""

import json
import re
import subprocess
import threading
import uuid
from datetime import datetime
from pathlib import Path

from flask import Flask, render_template, request, redirect, url_for, send_file
from werkzeug.utils import secure_filename

app = Flask(__name__)

# --- Configuration you may need to change ---
# Absolute path to your existing Snakemake project (Snakefile + Fig1..Fig7 folders).
LCR_PROJECT_DIR = Path("C:/Users/91892/Documents/Project_LCR")

# Folder where every submitted job gets its own subfolder. Lives next to this app.py.
JOBS_DIR = Path(__file__).parent / "jobs"
JOBS_DIR.mkdir(exist_ok=True)

ALLOWED_BED_EXTENSIONS = {".bed"}
ALLOWED_FASTA_EXTENSIONS = {".fasta", ".fa"}

# Human-friendly labels shown on the status page as each Snakemake rule starts.
# Any rule not listed here just falls back to showing its raw name.
RULE_LABELS = {
    "step1": "Fake step 1 (testing)",
    "step2": "Fake step 2 (testing)",
    "clean_fasta": "Cleaning FASTA headers",
    "prepare_bedtools_beds": "Preparing BED files for BEDTools",
    "fig1_extract_sequences": "Fig 1: extracting LCR sequences",
    "fig1_length_counts": "Fig 1: computing length & count stats",
    "fig1_entropy": "Fig 1: computing Shannon entropy",
    "fig1_coverage": "Fig 1: computing coverage stats",
    "fig1_aa_composition": "Fig 1: computing amino acid composition",
    "fig1A_plot": "Fig 1A: plotting length heatmap",
    "fig1B_plot": "Fig 1B: plotting coverage heatmap",
    "fig1C_plot": "Fig 1C: plotting count distribution",
    "fig1D_plot": "Fig 1D: plotting amino acid composition",
    "fig1E_plot": "Fig 1E: plotting entropy boxplot",
    "fig1_combined": "Fig 1: building combined figure",
    "fig2_multiinter": "Fig 2: computing multi-tool intersections",
    "fig2_consensus_beds": "Fig 2: building consensus BED files",
    "fig2_consensus_fastas": "Fig 2: extracting consensus sequences",
    "fig2_substring_motifs": "Fig 2: computing substring motifs",
    "fig2_filtered": "Fig 2: filtering best motifs",
    "fig2_peptide_counts": "Fig 2: counting peptide motifs",
    "fig2_entropy": "Fig 2: computing consensus entropy",
    "fig2_purity": "Fig 2: computing consensus purity",
    "fig2a_plot": "Fig 2A: plotting peptide motifs",
    "fig2b_plot": "Fig 2B: plotting entropy",
    "fig2c_plot": "Fig 2C: plotting purity",
    "fig3_fastas": "Fig 3: extracting per-tool sequences",
    "fig3_purity_metrics": "Fig 3: computing purity metrics",
    "fig3_plot": "Fig 3: plotting purity across tools",
    "fig4_jaccard_matrix": "Fig 4: computing Jaccard similarity",
    "fig4_plot": "Fig 4: plotting pairwise tool overlap",
    "fig5_metrics": "Fig 5: computing mutation metrics",
    "fig5_heatmaps": "Fig 5: plotting per-tool heatmaps",
    "fig5_combined": "Fig 5: building combined figure",
    "fig6_windows": "Fig 6: building sliding windows",
    "fig6_complexity": "Fig 6: computing window complexity",
    "fig6_classification": "Fig 6: classifying windows",
    "fig6_annotate": "Fig 6: annotating HCR/LCR/CBR regions",
    "fig6_confusion": "Fig 6: computing TP/FP/FN/TN per tool",
    "fig6_reference_metrics": "Fig 6: computing reference metrics",
    "fig6_tpr_fpr": "Fig 6: computing TPR/FPR tables",
    "fig6_plots": "Fig 6: plotting TPR/FPR figures",
    "fig6_combined": "Fig 6: building combined figure",
    "fig7_summaries": "Fig 7: building summary tables",
    "fig7_plot": "Fig 7: plotting summary bar charts",
    "generate_report": "Building final PDF report",
    "all": "Finishing up",
}


def has_allowed_extension(filename, allowed_extensions):
    return Path(filename).suffix.lower() in allowed_extensions


def update_status(job_root, **fields):
    """Read-modify-write the job's status.json with whatever fields are passed in."""
    status_path = job_root / "status.json"
    status = {}
    if status_path.exists():
        status = json.loads(status_path.read_text())
    status.update(fields)
    status["updated_at"] = datetime.now().strftime("%H:%M:%S")
    status_path.write_text(json.dumps(status, indent=2))


def write_job_config(job_root, organism, base_dir):
    """Write a per-job config.yaml Snakemake will use, with absolute paths."""
    config_path = job_root / "config.yaml"
    config_text = (
        f"organism: {organism}\n"
        f"base_dir: {base_dir.as_posix()}\n"
        f"scripts_dir: {LCR_PROJECT_DIR.as_posix()}\n"
        f"combined_figures: true\n"
    )
    config_path.write_text(config_text)
    return config_path


def run_pipeline(job_root, config_path):
    """Runs in a background thread: launches Snakemake, watches its output for
    'rule <name>:' lines, and updates status.json as each stage starts."""

    update_status(job_root, state="running", current_rule=None,
                  current_label="Starting pipeline...")

    cmd = [
        "snakemake",
        "--snakefile", str(LCR_PROJECT_DIR / "Snakefile"),
        "--configfile", str(config_path),
        "--cores", "4",
    ]

    log_path = job_root / "log.txt"
    rule_pattern = re.compile(r"^(?:local)?rule (\w+):")

    try:
        with open(log_path, "w") as log_file:
            process = subprocess.Popen(
                cmd,
                cwd=str(LCR_PROJECT_DIR),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            for line in process.stdout:
                log_file.write(line)
                log_file.flush()
                match = rule_pattern.match(line.strip())
                if match:
                    rule_name = match.group(1)
                    label = RULE_LABELS.get(rule_name, f"Running step: {rule_name}")
                    update_status(job_root, current_rule=rule_name, current_label=label)
            process.wait()

        if process.returncode == 0:
            update_status(job_root, state="done", current_label="Pipeline complete")
        else:
            update_status(job_root, state="failed",
                          current_label="Pipeline failed — check the log")
    except Exception as e:
        update_status(job_root, state="failed", current_label=f"Error launching pipeline: {e}")


@app.route("/", methods=["GET"])
def upload_form():
    return render_template("upload.html")

# --- add these three ---
@app.route("/about")
def index():
    return render_template("index.html")

@app.route("/docs")
def docs_page():
    return render_template("docs.html")

@app.route("/results")
def results_page():
    return render_template("results.html")
# ------------------------
@app.route("/submit", methods=["POST"])
def submit_job():
    organism = request.form.get("organism", "").strip()
    if not organism:
        return "Error: organism name is required. Go back and try again.", 400

    fasta_file = request.files.get("fasta_file")
    bed_files = request.files.getlist("bed_files")

    if fasta_file is None or fasta_file.filename == "":
        return "Error: a FASTA file is required. Go back and try again.", 400
    if not has_allowed_extension(fasta_file.filename, ALLOWED_FASTA_EXTENSIONS):
        return f"Error: '{fasta_file.filename}' is not a .fasta/.fa file.", 400
    if not bed_files or bed_files[0].filename == "":
        return "Error: at least one BED file is required.", 400
    for bed_file in bed_files:
        if not has_allowed_extension(bed_file.filename, ALLOWED_BED_EXTENSIONS):
            return f"Error: '{bed_file.filename}' is not a .bed file.", 400

    job_id = uuid.uuid4().hex[:10]
    job_root = JOBS_DIR / job_id
    organism_dir = job_root / organism
    bed_dir = organism_dir / f"bed_{organism}"
    bed_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = organism_dir / f"{organism}.fasta"
    fasta_file.save(fasta_path)

    for bed_file in bed_files:
        safe_name = secure_filename(bed_file.filename)
        bed_file.save(bed_dir / safe_name)

    config_path = write_job_config(job_root, organism, organism_dir)
    update_status(job_root, state="queued", current_label="Waiting to start...")

    thread = threading.Thread(target=run_pipeline, args=(job_root, config_path), daemon=True)
    thread.start()

    return redirect(url_for("job_status", job_id=job_id))


@app.route("/status/<job_id>")
def job_status(job_id):
    job_root = JOBS_DIR / job_id
    status_path = job_root / "status.json"
    if not status_path.exists():
        return "Job not found.", 404
    status = json.loads(status_path.read_text())
    return render_template("status.html", job_id=job_id, status=status)


@app.route("/retry/<job_id>", methods=["POST"])
def retry_job(job_id):
    """Re-run the same job using its already-saved config.yaml and uploaded files —
    no re-upload needed."""
    job_root = JOBS_DIR / job_id
    config_path = job_root / "config.yaml"
    if not config_path.exists():
        return "Job not found.", 404

    update_status(job_root, state="queued", current_label="Waiting to restart...")
    thread = threading.Thread(target=run_pipeline, args=(job_root, config_path), daemon=True)
    thread.start()
    return redirect(url_for("job_status", job_id=job_id))


@app.route("/report/<job_id>")
def job_report(job_id):
    job_root = JOBS_DIR / job_id
    organism_dirs = [p for p in job_root.iterdir() if p.is_dir()]
    if not organism_dirs:
        return "Report not found.", 404
    organism_dir = organism_dirs[0]
    report_path = organism_dir / "reports" / "Benchmark_Report.pdf"
    if not report_path.exists():
        return "Report not found — the pipeline may not have finished generating it.", 404
    return send_file(report_path, as_attachment=True, download_name="Benchmark_Report.pdf")


@app.route("/log/<job_id>")
def job_log(job_id):
    log_path = JOBS_DIR / job_id / "log.txt"
    if not log_path.exists():
        return "No log yet.", 404
    return f"<pre>{log_path.read_text()}</pre>"


if __name__ == "__main__":
    app.run(debug=True)