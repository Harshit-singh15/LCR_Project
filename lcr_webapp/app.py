"""
LCR Benchmarking Tool — Steps 1 & 2: Upload page + background pipeline runner + status page.

Run this with:
    python app.py

Then open http://127.0.0.1:5000 in your browser.
"""
 
import io
import zipfile
import json
import re
import subprocess
import threading
import uuid
from datetime import datetime
from pathlib import Path
import os
import shutil
import tempfile
import sys
 

from flask import Flask, render_template, request, redirect, url_for, send_file
from werkzeug.utils import secure_filename

app = Flask(__name__)

# --- Configuration you may need to change ---
# ✅ FIXED: Dynamically finds the Project_LCR folder relative to this file
# __file__ is 'app.py' -> .parent is 'lcr_webapp' -> .parent.parent is 'Project_LCR'
LCR_PROJECT_DIR = Path(__file__).resolve().parent.parent

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


# --- Config: where your converter scripts live ---
CONVERTER_SCRIPTS_DIR = LCR_PROJECT_DIR / "converter_to_bed"
 
 
def run_converter(script_name, file_obj, cmd_args, output_filename):
    """
    script_name:     e.g. "convert_seg_alcor.py"
    file_obj:        The uploaded Flask file object (e.g., f)
    cmd_args:        List of extra CLI args to pass AFTER input/output paths if needed
    output_filename: The name of the file expected to be created
    """
    # 1. Create the authoritative temp directory inside the helper
    work_dir = Path(tempfile.mkdtemp(prefix="lcrbench_conv_"))
    
    # 2. Save the file right inside this newly created directory
    input_file_path = work_dir / secure_filename(file_obj.filename)
    file_obj.save(input_file_path)
    
    # 3. Establish the destination output path
    output_path = work_dir / output_filename
 
    # 4. Resolve the absolute path to your converter script
    script_path = CONVERTER_SCRIPTS_DIR / script_name
 
    # 5. Build your structural sys.argv mapping array:
    # sys.argv[1] = input file path, sys.argv[2] = output file path
    cmd = ["python", str(script_path), str(input_file_path), str(output_path)] + cmd_args
 
    # Execute the command inside the working directory scope
    result = subprocess.run(
        cmd, cwd=str(work_dir),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    
    if result.returncode != 0:
        shutil.rmtree(work_dir, ignore_errors=True)
        raise RuntimeError(result.stderr or result.stdout or "Converter script failed.")
 
    if not output_path.exists():
        shutil.rmtree(work_dir, ignore_errors=True)
        raise RuntimeError(f"Expected output '{output_filename}' was not created.")
 
    return output_path, work_dir

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

def find_organism_dir(job_id):
    """Same lookup job_report() already uses — the one subfolder under jobs/<job_id>/."""
    job_root = JOBS_DIR / job_id
    if not job_root.exists():
        return None
    organism_dirs = [p for p in job_root.iterdir() if p.is_dir()]
    return organism_dirs[0] if organism_dirs else None
 

def run_pipeline(job_root, config_path):
    """Runs in a background thread: launches Snakemake, watches its output for
    'rule <name>:' lines, and updates status.json as each stage starts."""

    update_status(job_root, state="running", current_rule=None,
                  current_label="Starting pipeline...")

    # --- temporary debug ---
    debug_msg = (
        f"LCR_PROJECT_DIR={LCR_PROJECT_DIR} exists={LCR_PROJECT_DIR.exists()} | "
        f"Snakefile exists={ (LCR_PROJECT_DIR / 'Snakefile').exists() }"
    )
    print(debug_msg, flush=True)
    update_status(job_root, current_label=debug_msg)
    # --- end debug ---

    # ✅ FIX 1: Dynamically set cores. Use 1 core on cloud hosting, 4 locally.
    # Render environments usually set standard cloud variables, or we can just default to 1 for safety.
    cores = "1" if os.environ.get("RENDER") else "4"

    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(LCR_PROJECT_DIR / "Snakefile"),
        "--configfile", str(config_path),
        "--cores", cores,
    ]

    log_path = job_root / "log.txt"
    # ✅ FIX 2: Refined regex pattern to match Snakemake's true terminal logging format
    rule_pattern = re.compile(r"^(?:local)?rule\s+(\w+):")

    try:
        with open(log_path, "w") as log_file:
            process = subprocess.Popen(
                cmd,
                cwd=str(LCR_PROJECT_DIR),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1, # Line buffered so it reads stdout in real-time
            )
            
            # Read stdout line by line as it prints
            for line in process.stdout:
                log_file.write(line)
                log_file.flush()
                
                match = rule_pattern.search(line.strip())
                if match:
                    rule_name = match.group(1)
                    label = RULE_LABELS.get(rule_name, f"Running step: {rule_name}")
                    update_status(job_root, current_rule=rule_name, current_label=label)
            
            process.wait()

        if process.returncode == 0:
            update_status(job_root, state="done", current_label="Pipeline complete")
        else:
            update_status(job_root, state="failed",
                          current_label="Pipeline failed — check log.txt inside your job folder")
    except FileNotFoundError as e:
        update_status(job_root, state="failed", current_label=f"Error: {e}")
        
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

'''@app.route("/results")
def results_page():
    return render_template("results.html")'''


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


@app.route("/convert")
def convert_page():
    return render_template("convert.html")
 
@app.route("/convert/seg_alcor", methods=["POST"])
def convert_seg_alcor():
    f = request.files.get("masked_fasta")
    if not f or f.filename == "":
        return "Error: upload a soft-masked FASTA file.", 400
 
    output_name = "seg_alcor_output.bed"
 
    try:
        # Pass the script name, the raw file object, empty extra args list, and output name
        output_path, _ = run_converter("convert_seg_alcor.py", f, [], output_name)
    except RuntimeError as e:
        return f"Conversion failed: {e}", 500
 
    return send_file(output_path, as_attachment=True, download_name="seg_alcor.bed")

@app.route("/convert/xstream", methods=["POST"])
def convert_xstream():
    f1 = request.files.get("xstream_html")
    if not f1 or f1.filename == "":
        return "Error: No file uploaded.", 400
        
    if not f1.filename.lower().endswith('.html'):
        return "Error: Please upload a valid .html file.", 400
 
    output_name = "xstream_output.bed"
 
    try:
        output_path, _ = run_converter("convert_xstream.py", f1, [], output_name)
    except RuntimeError as e:
        return f"Conversion failed: {e}", 500
 
    return send_file(output_path, as_attachment=True, download_name="xstream.bed")

@app.route("/convert/flps", methods=["POST"])
def convert_flps():
    f = request.files.get("flps_out")
    if not f or f.filename == "":
        return "Error: upload an fLPS .out file.", 400
 
    output_name = "flps_output.bed"
 
    try:
        output_path, _ = run_converter("convert_flps.py", f, [], output_name)
    except RuntimeError as e:
        return f"Conversion failed: {e}", 500
 
    return send_file(output_path, as_attachment=True, download_name="flps.bed")

@app.route("/convert/treks", methods=["POST"])
def convert_treks():
    f = request.files.get("treks_tsv")
    if not f or f.filename == "":
        return "Error: upload a T-REKS .tsv file.", 400
 
    output_name = "treks_output.bed"
 
    try:
        output_path, _ = run_converter("convert_treks.py", f, [], output_name)
    except RuntimeError as e:
        return f"Conversion failed: {e}", 500
 
    return send_file(output_path, as_attachment=True, download_name="treks.bed")


@app.route("/results")
@app.route("/results/<job_id>")
def results_page(job_id=None):
    return render_template("results.html", job_id=job_id)
 
 
@app.route("/download/figure/<job_id>/<int:fig_num>")
def download_figure(job_id, fig_num):
    organism_dir = find_organism_dir(job_id)
    if organism_dir is None:
        return "Job not found.", 404
 
    fig_dir = organism_dir / "Fig_outputs" / f"Fig{fig_num}"
    if not fig_dir.exists():
        return f"No output found for Figure {fig_num}.", 404
 
    pngs = sorted(fig_dir.glob("*.png"))
    if not pngs:
        return f"No PNGs found for Figure {fig_num}.", 404
 
    if len(pngs) == 1:
        return send_file(pngs[0], as_attachment=True, download_name=pngs[0].name)

    # Preview requests (the <img> tag sends ?inline=1) get just the combined PNG,
    # not a zip — pick the one with "combined" in its name, or fall back to the last one.
    if request.args.get("inline"):
        if fig_num in (1, 5):
            combined = next((p for p in pngs if "combined" in p.name.lower()), pngs[0])
        elif fig_num == 6:
        # Takes the 2nd png in the folder, falling back to the 1st if fewer than 2 exist
            combined = pngs[1] if len(pngs) > 1 else pngs[0]
        else:
            combined = pngs[0]
        return send_file(combined, as_attachment=False, download_name=combined.name)
        
    # Real downloads still get every panel PNG for this figure, zipped.
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in pngs:
            zf.write(p, arcname=p.name)
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name=f"Fig{fig_num}.zip", mimetype="application/zip")
 
@app.route("/download/figures/<job_id>")
def download_all_figures(job_id):
    organism_dir = find_organism_dir(job_id)
    if organism_dir is None:
        return "Job not found.", 404
 
    fig_outputs = organism_dir / "Fig_outputs"
    if not fig_outputs.exists():
        return "No figures found for this job.", 404
 
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in range(1, 8):
            fig_dir = fig_outputs / f"Fig{n}"
            if fig_dir.exists():
                for p in fig_dir.glob("*.png"):
                    zf.write(p, arcname=f"Fig{n}/{p.name}")
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name="LCRBench_figures.zip", mimetype="application/zip")
 
 
@app.route("/download/intermediate/<job_id>")
def download_intermediate(job_id):
    organism_dir = find_organism_dir(job_id)
    if organism_dir is None:
        return "Job not found.", 404
 
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        found_any = False
        for n in range(1, 8):
            data_dir = organism_dir / f"dataforFig{n}"
            if data_dir.exists():
                for p in data_dir.rglob("*"):
                    if p.is_file():
                        zf.write(p, arcname=f"dataforFig{n}/{p.relative_to(data_dir)}")
                        found_any = True
        if not found_any:
            buf.close()
            return "No intermediate files found for this job.", 404
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name="LCRBench_intermediate_files.zip", mimetype="application/zip")
 
 
@app.route("/download/complete/<job_id>")
def download_complete(job_id):
    """Report + all figures + all intermediate data + log + config, in one zip."""
    organism_dir = find_organism_dir(job_id)
    if organism_dir is None:
        return "Job not found.", 404
 
    job_root = JOBS_DIR / job_id
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        report_path = organism_dir / "reports" / "Benchmark_Report.pdf"
        if report_path.exists():
            zf.write(report_path, arcname="Benchmark_Report.pdf")
 
        fig_outputs = organism_dir / "Fig_outputs"
        for n in range(1, 8):
            fig_dir = fig_outputs / f"Fig{n}"
            if fig_dir.exists():
                for p in fig_dir.glob("*.png"):
                    zf.write(p, arcname=f"figures/Fig{n}/{p.name}")
 
            data_dir = organism_dir / f"dataforFig{n}"
            if data_dir.exists():
                for p in data_dir.rglob("*"):
                    if p.is_file():
                        zf.write(p, arcname=f"intermediate/dataforFig{n}/{p.relative_to(data_dir)}")
 
        log_path = job_root / "log.txt"
        if log_path.exists():
            zf.write(log_path, arcname="log.txt")
 
        config_path = job_root / "config.yaml"
        if config_path.exists():
            zf.write(config_path, arcname="config.yaml")
 
        status_path = job_root / "status.json"
        if status_path.exists():
            zf.write(status_path, arcname="status.json")
 
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name="LCRBench_complete_results.zip", mimetype="application/zip")

if __name__ == "__main__":
    app.run(debug=True)