# LCRBenchmark

**LCRBenchmark** is a benchmarking toolkit for comparing Low Complexity Region (LCR) detection tools (SEG, ALCOR, XSTREAM, T-REKS, fLPS, fLPS2, LCR-FINDER, and others) across a proteome. It produces a full suite of comparison figures (length, coverage, entropy, motif purity, Jaccard overlap, TPR/FPR against a reference set, etc.) and a final PDF benchmark report.

The tool can be used in **two ways**:

1. **Command line** — via a Snakemake pipeline (full control, reproducible, scriptable)
2. **Web app (GUI)** — a local Flask-based interface for uploading files and getting outputs without touching the command line

> The work you need is on the **`website`** branch — this README lives there. The `main` branch contains the original command-line-only version and additional supplementary figures of this project. **`website` is not the default branch**, so make sure to check it out explicitly using the commands below.

---

## Table of Contents

- [Cloning the repo](#cloning-the-repo)
- [Option 1: Command-line (Snakemake) pipeline](#option-1-command-line-snakemake-pipeline)
- [Option 2: Web app (GUI)](#option-2-web-app-gui)
- [Try it with sample data (8 organisms)](#try-it-with-sample-data-8-organisms)
- [Live demo](#live-demo)
- [Project report](#project-report)

---

## Cloning the repo

This repository uses **Git LFS** (Large File Storage) for large data/output files, so a full clone can be large. Pick whichever combination suits your storage/bandwidth needs. Since `website` is **not** the default branch, always specify it with `-b website`.

### 1. Full clone (all branches + all LFS files)

```bash
git lfs install
git clone https://github.com/Harshit-singh15/LCR_Project
```

### 2. Full clone, without LFS files (pointers only, saves storage)

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/Harshit-singh15/LCR_Project
```

### 3. Clone only the `website` branch (all LFS files)

```bash
git lfs install
git clone -b website --single-branch https://github.com/Harshit-singh15/LCR_Project
```

### 4. Clone only the `website` branch, without LFS files (smallest possible clone)

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone -b website --single-branch https://github.com/Harshit-singh15/LCR_Project
```

> If you used one of the "without LFS" options and later need the actual large files, run this from inside the repo:
> ```bash
> git lfs pull                                     # pulls all LFS files tracked on the current branch
> git lfs pull --include="path/to/specific/folder"  # pulls only what you need
> ```
> This is also how you download the pre-computed BED files for the [8 sample organisms](#try-it-with-sample-data-8-organisms) described below.

---

## Option 1: Command-line (Snakemake) pipeline

The CLI pipeline is defined in the `Snakefile` in the project root, and produces Figures 1–7 plus a final `Benchmark_Report.pdf`. (Fig 8 is excluded — it depends on external DisProt/SIFTS data and is a separate workflow.)

It's recommended to run this inside a **virtual environment**:

```bash
python -m venv venv
source venv/bin/activate        # Windows CP: venv\Scripts\activate
pip install -r requirements.txt # Windows PS : venv\Scripts\Activate.ps1
```

### Step 1 — Add your data folder

In the **root directory**, add a folder named after your organism (this name is used throughout the pipeline as a prefix), for example:

```
ecoli/
├── ecoli.fasta                # proteome FASTA file for the organism
└── bed_ecoli/                 # 1-based BED outputs from each LCR detection tool
    ├── seg.bed
    ├── alcor.bed
    ├── xstream.bed
    ├── treks.bed
    ├── flps.bed
    ├── flps2.bed
    ├── lcrfinder.bed
    └── ...                    # any other tool you want to benchmark
```

**BED filenames don't need to follow any fixed convention.** Name each file after the tool (e.g. `seg.bed`, `alcor_mode1.bed`, `myowntool.bed`) — the filename is used directly as the **legend label** in the figures, so pick something short and readable. You're not limited to the tools listed above; you can drop in BED files from **any** LCR detection tool you want compared, as long as it's in the same 1-based BED format inside `bed_<organism>/`.

### Step 2 — Converting raw tool output to 1-based BED format

If your tool's raw output isn't already in the required 1-based BED format, use the converters in `converter_to_bed/`:

| Script | Converts output from |
|---|---|
| `convert_seg_alcor.py` | SEG / ALCOR |
| `convert_flps.py` | fLPS / fLPS2 |
| `convert_treks.py` | T-REKS |
| `convert_xstream.py` | XSTREAM |
| `convert_dotplot.py` | Not a detection tool — data referenced from a paper; see the [project report](#project-report) for details |

**LCR-FINDER** uses a separate MATLAB-based conversion path (see below).

### Step 3 — Converting LCR-FINDER output (MATLAB scripts)

The `matlab_scripts/` folder contains two functions used specifically for LCR-FINDER:

**1. `processFasta`** — prepares your FASTA file for LCR-FINDER:

```matlab
processFasta(fastaPath, targetFolder, outputName)
```

- `fastaPath` — path to your proteome FASTA file
- `targetFolder` — output folder for the processed files
- `outputName` — should match your organism name
The scripts contains the variable `toolsFolder` — Add the path where LcrFinder is located in your system.

Run LCR-FINDER on the output of this step.

**2. `createLcrBedFile`** — converts LCR-FINDER's output into a 1-based BED file (`ProteinID  Start  End`):

```matlab
createLcrBedFile( ...
    'mouse_proteome.fasta', ...
    'LCRFinder_Output', ...
    'lcrfinder.bed');
```

Place the resulting `lcrfinder.bed` inside your `bed_<organism>/` folder alongside the other tool outputs.

### Step 4 — Reference material (optional, for orientation)

- **Command sheet (PNG)** — a quick-reference image of the commands used to run each detection tool:

  `![Command sheet](cheatsheet.png)`


- **Workflow DAG (PDF)** — the Snakemake DAG graph showing how all rules/figures connect:

  `[Workflow DAG](dag.pdf)`

### Step 5 — Edit `config.yaml`

Update `config.yaml` in the project root to match your organism and folder setup:

```yaml
# Organism name — matches the folder name and file prefixes throughout
# e.g. ecoli.fasta, bed_ecoli/, dataforFig1/, etc.
organism: yeast

# Root folder for this organism's data (relative to the Snakefile, or absolute)
base_dir: yeast

# Folder containing your converted Python scripts, organized same as before:
scripts_dir: .

# Whether to also build the combined (multi-panel) figures for Fig1, Fig5, Fig6
combined_figures: true
```

### Step 6 — Run the pipeline

**Dry run** (see the execution plan without running anything):

```bash
snakemake -n
```

**Visualize the DAG** (requires graphviz):

```bash
snakemake --dag | dot -Tpng > dag.png
```

**Run using all available cores:**

```bash
snakemake --cores all
```

**Run using a specific number of cores** (e.g. 4):

```bash
snakemake --cores 4
```

**Run only a specific target** (e.g. just Fig 7):

```bash
snakemake --cores 4 fig7
```

Outputs are written under `<base_dir>/Fig_outputs/Fig1` through `Fig7`, and the final combined report at `<base_dir>/reports/Benchmark_Report.pdf`.

---

## Option 2: Web app (GUI)

If you'd rather not use the command line, use the local web interface.

It's recommended to run this inside a **virtual environment** as well:

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
```

### Step 1 — Navigate to the web app folder

From the project root:

```bash
cd LCR_Project/lcr_webapp
```

### Step 2 — Install dependencies (if not already installed)

```bash
pip install -r requirements.txt
```

### Step 3 — Run the app

```bash
python app.py
```

### Step 4 — Open it in your browser

Go to the local host address shown in the terminal (typically `http://127.0.0.1:5000`).

### Step 5 — Upload and run

Navigate to the **Upload** section in the app, upload your FASTA and BED files, and the app will generate the benchmark outputs for you — no command line needed.

---

## Try it with sample data (8 organisms)

If you don't have your own data yet and just want to try the tool, pre-computed data is provided for **8 organisms**:

- Zebrafish
- *C. elegans*
- *E. coli*
- Yeast
- Mouse
- Human
- *Arabidopsis*
- Fruit fly (*Drosophila*)

For each organism, the full set of intermediate pipeline files is available in the repo (fetched via LFS). For quick practice, you only need the **FASTA and BED files** — pull just those instead of everything:

```bash
git lfs pull --include="<organism_folder>/<organism>.fasta,<organism_folder>/bed_<organism>/*"
```

In the **web app**, the sample-data section lets you directly download the FASTA and BED files for any of these 8 organisms and try the tool immediately without preparing your own dataset.

---

## Live demo

A hosted version of the web app is available here:

`https://project-lcr.onrender.com/`

---

## Project report

The full project/internship report is included in this repository. Add a relative link to it once it's in place, for example:

```markdown
[Project Report (PDF)](LCR_Benchmarking_Internship_Report.pdf)
```

GitHub will render this as a clickable link that opens/downloads the file directly from the repo — no external hosting needed.

---

## Notes

- If cloning with `GIT_LFS_SKIP_SMUDGE=1`, remember to run `git lfs pull` before running the pipeline if any required data/scripts are tracked via LFS.
- Make sure `organism`, `base_dir`, and the `bed_<organism>` folder naming stay consistent with each other — the pipeline relies on this naming convention throughout.
- BED filenames are free-form and used as figure legends — name them clearly after the tool that produced them.
