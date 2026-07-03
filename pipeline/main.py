"""
main.py

Main pipeline for running all LCR detection tools.

Only edit:

1. INPUT_FASTA
2. OUTPUT_DIR
"""

from pathlib import Path

from tools.tools import (
    run_seg,
    run_flps,
    run_flps2,
    run_alcor,
    run_xstream,
    run_treks,
)

# ==========================================================
# Input / Output
# ==========================================================

BASE_DIR = Path.cwd().resolve()

INPUT_FASTA = (
    BASE_DIR
    / "ecoli"
    / "ecoli.fasta"
)

OUTPUT_DIR = (
    BASE_DIR
    / "pipeline"
    / "toolsoutput"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

# ==========================================================
# Run Pipeline
# ==========================================================

results = []

print("\nStarting LCR Benchmarking Pipeline...\n")

# ----------------------------------------------------------
# SEG
# ----------------------------------------------------------

results.append(
    run_seg(
        INPUT_FASTA,
        OUTPUT_DIR,
        window=12,
        locut=2.2,
        hicut=2.5,
    )
)

results.append(
    run_seg(
        INPUT_FASTA,
        OUTPUT_DIR,
        window=15,
        locut=1.9,
        hicut=2.5,
    )
)

results.append(
    run_seg(
        INPUT_FASTA,
        OUTPUT_DIR,
        window=15,
        locut=1.5,
        hicut=1.8,
    )
)

# ----------------------------------------------------------
# fLPS
# ----------------------------------------------------------

results.append(
    run_flps(
        INPUT_FASTA,
        OUTPUT_DIR,
    )
)

results.append(
    run_flps(
        INPUT_FASTA,
        OUTPUT_DIR,
        m=5,
        M=25,
        t=0.00001,
    )
)

# ----------------------------------------------------------
# fLPS2
# ----------------------------------------------------------

results.append(
    run_flps2(
        INPUT_FASTA,
        OUTPUT_DIR,
    )
)

results.append(
    run_flps2(
        INPUT_FASTA,
        OUTPUT_DIR,
        m=5,
        M=25,
        t=0.00001,
    )
)

# ----------------------------------------------------------
# AlcoR
# ----------------------------------------------------------

results.append(
    run_alcor(
        INPUT_FASTA,
        OUTPUT_DIR,
        mode=1,
    )
)

results.append(
    run_alcor(
        INPUT_FASTA,
        OUTPUT_DIR,
        mode=2,
    )
)

# ----------------------------------------------------------
# XSTREAM
# ----------------------------------------------------------

results.append(
    run_xstream(
        INPUT_FASTA,
        OUTPUT_DIR,
        mode=1,
    )
)

# ----------------------------------------------------------
# T-REKS
# ----------------------------------------------------------

results.append(
    run_treks(
        INPUT_FASTA,
        OUTPUT_DIR,
    )
)

# ==========================================================
# Summary
# ==========================================================
completed = sum(1 for result in results if result.get("status") == "SUCCESS")
failed = sum(1 for result in results if result.get("status") in {"FAILED", "ERROR"})
timed_out = sum(1 for result in results if result.get("status") == "TIMEOUT")

for result in results:
    tool = result.get("tool", "UNKNOWN")
    status = result.get("status", "UNKNOWN")
    runtime = result.get("runtime", 0.0)

    print(f"{tool:<24} {status:<10} {runtime:>8.2f} s")
    if status != "SUCCESS":
        print(f"    Reason: {result.get('reason', 'unknown')}")

print("\n" + "=" * 60)
print("PIPELINE SUMMARY")
print("=" * 60)

print(f"Completed : {completed}")
print(f"Failed    : {failed}")
print(f"Timed Out : {timed_out}")

print("\nPipeline Finished.")