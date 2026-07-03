"""
tools.py

Wrapper functions for running external LCR detection tools from Python.

Supported tools
---------------
1. SEG
2. fLPS
3. fLPS2
4. AlcoR
5. XSTREAM
6. T-REKS

Features
--------
✓ Automatic Windows → WSL path conversion
✓ Runtime measurement
✓ Timeout support
✓ Error handling
✓ Pipeline never crashes because one tool fails
✓ Automatic logging
"""

from pathlib import Path
from shlex import shlex
import subprocess
import time
import traceback
from tools.config import (
    SEG,
    FLPS,
    FLPS2,
    ALCOR,
    XSTREAM,
    TREKS,
    CLUSTALW,
    DEFAULT_TIMEOUT,
)


# ==========================================================
# Helper Functions
# ==========================================================

def convert_to_wsl_path(windows_path: Path) -> str:
    """Converts a Windows absolute Path object into a WSL /mnt/c/... string."""
    # 1. Force the path to use forward slashes (e.g., "C:/Users/...")
    posix_path = windows_path.resolve().as_posix()
    
    # 2. Replace the drive letter "C:" with "/mnt/c"
    # This handles both uppercase 'C:' and lowercase 'c:'
    if posix_path.startswith("C:"):
        wsl_path = posix_path.replace("C:", "/mnt/c", 1)
    elif posix_path.startswith("c:"):
        wsl_path = posix_path.replace("c:", "/mnt/c", 1)
    else:
        wsl_path = posix_path
        
    return wsl_path

def get_prefix(input_fasta: Path) -> str:
    """
    Returns filename without extension.

    Example
    -------
    ecoli.fasta

    becomes

    ecoli
    """

    return Path(input_fasta).stem


# ----------------------------------------------------------


def build_output(
    output_dir: Path,
    input_fasta: Path,
    suffix: str,
) -> Path:
    """
    Automatically build output filename.

    Example
    -------

    input

        ecoli.fasta

    suffix

        seg.fa

    output

        output_dir/ecoli_seg.fa
    """

    prefix = get_prefix(input_fasta)

    return output_dir / f"{prefix}_{suffix}"


# ----------------------------------------------------------
def run_command(
    tool_name: str,
    command,
    timeout: int = DEFAULT_TIMEOUT,
):
    print("\n" + "=" * 70)
    print(f"Running : {tool_name}")
    print("=" * 70)

    print("\nCommand\n")
    if isinstance(command, (list, tuple)):
        print(" ".join(map(str, command)))
    else:
        print(command)

    start = time.time()

    try:
        subprocess.run(
            command,
            check=True,
            timeout=timeout,
        )

        elapsed = time.time() - start

        print("\nStatus : SUCCESS")
        print(f"Runtime : {elapsed:.2f} sec")

        return {
            "tool": tool_name,
            "status": "SUCCESS",
            "runtime": elapsed,
        }

    except subprocess.TimeoutExpired as e:

        elapsed = time.time() - start

        print("\nStatus : TIMEOUT")
        print(f"Exceeded {timeout} seconds")
        print(f"Runtime : {elapsed:.2f} sec")

        return {
            "tool": tool_name,
            "status": "TIMEOUT",
            "runtime": elapsed,
            "reason": str(e),
        }

    except subprocess.CalledProcessError as e:

        elapsed = time.time() - start

        print("\nStatus : FAILED")
        print(f"Runtime : {elapsed:.2f} sec")
        print("\nReason\n")
        print(e)

        return {
            "tool": tool_name,
            "status": "FAILED",
            "runtime": elapsed,
            "reason": str(e),
        }

    except Exception as e:

        elapsed = time.time() - start

        print("\nUnexpected Error\n")
        traceback.print_exc()
        print(f"\nRuntime : {elapsed:.2f} sec")

        return {
            "tool": tool_name,
            "status": "ERROR",
            "runtime": elapsed,
            "reason": str(e),
        }
# ----------------------------------------------------------


def run_bash_command(
    tool_name: str,
    bash_command: str,
    timeout: int = DEFAULT_TIMEOUT,
):
    """
    Execute a bash command.

    Used for tools which need

        >
        pipes
        shell operators

    Example

        fLPS

        XSTREAM
    """

    cmd = [
        "wsl",
        "bash",
        "-c",
        bash_command,
    ]

    return run_command(
        tool_name,
        cmd,
        timeout,
    )

# ==========================================================
# SEG
# ==========================================================

def run_seg(
    input_fasta: Path,
    output_dir: Path,
    window: int,
    locut: float,
    hicut: float,
    timeout: int = DEFAULT_TIMEOUT,
):
    """
    Run NCBI SEG (segmasker).

    Parameters
    ----------
    input_fasta : Path
        Input FASTA file.

    output_dir : Path
        Directory where output will be written.

    window : int
        SEG window size.

    locut : float
        Low complexity cutoff.

    hicut : float
        High complexity cutoff.
    """

    output_file = build_output(
        output_dir,
        input_fasta,
        f"seg_w{window}_l{locut}_h{hicut}.fa"
    )

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output = convert_to_wsl_path(output_file)

    cmd = [
        "wsl",
        SEG,
        "-in", wsl_input,
        "-window", str(window),
        "-locut", str(locut),
        "-hicut", str(hicut),
        "-out", wsl_output,
        "-outfmt", "fasta",
    ]

    return run_command(
        tool_name=f"SEG (window={window}, locut={locut}, hicut={hicut})",
        command=cmd,
        timeout=timeout,
    )


# ==========================================================
# fLPS
# ==========================================================
def run_flps(
    input_fasta: Path,
    output_dir: Path,
    m=None,
    M=None,
    t=None,
    timeout: int = DEFAULT_TIMEOUT,
):
    if m is None:
        suffix = "flps_default.out"
    else:
        suffix = "flps_strict.out"

    output_file = build_output(output_dir, input_fasta, suffix)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output = convert_to_wsl_path(output_file)

    command = f"{FLPS}"
    if m is not None:
        command += f" -m {m}"
    if M is not None:
        command += f" -M {M}"
    if t is not None:
        command += f" -t {t}"
    command += f" {shlex.quote(wsl_input)} > {shlex.quote(wsl_output)}"

    return run_bash_command(
        tool_name=suffix,
        bash_command=command,
        timeout=timeout,
    )
# ==========================================================
# fLPS2
# ==========================================================

def run_flps2(
    input_fasta: Path,
    output_dir: Path,
    m=None,
    M=None,
    t=None,
    timeout: int = DEFAULT_TIMEOUT,
):
    if m is None:
        suffix = "flps2_default.out"
    else:
        suffix = "flps2_strict.out"

    output_file = build_output(output_dir, input_fasta, suffix)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output = convert_to_wsl_path(output_file)

    command = f"{FLPS2}"
    if m is not None:
        command += f" -m {m}"
    if M is not None:
        command += f" -M {M}"
    if t is not None:
        command += f" -t {t}"
    command += f" {shlex.quote(wsl_input)} > {shlex.quote(wsl_output)}"

    return run_bash_command(
        tool_name=suffix,
        bash_command=command,
        timeout=timeout,
    )
# ==========================================================
# AlcoR
# ==========================================================

def run_alcor(
    input_fasta: Path,
    output_dir: Path,
    mode: int,
    timeout: int = DEFAULT_TIMEOUT,
):
    """
    Run AlcoR mapper.

    Parameters
    ----------
    input_fasta : Path
        Input FASTA file.

    output_dir : Path
        Output directory.

    mode : int
        1 or 2.

    timeout : int
        Timeout in seconds.
    """

    if mode == 1:

        mapper = "5:20:0:0:10:0.9/3:10:0.9"
        suffix = "alcor_mode1_masked.fa"

    elif mode == 2:

        mapper = "5:10:0:0:10:0.9/1:1:0.9"
        suffix = "alcor_mode2_masked.fa"

    else:

        raise ValueError("AlcoR mode must be 1 or 2.")

    output_file = build_output(
        output_dir,
        input_fasta,
        suffix,
    )

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output = convert_to_wsl_path(output_file)

    cmd = [
        "wsl",
        ALCOR,
        "mapper",
        "-v",
        "-m",
        mapper,
        "-w",
        "5",
        "-k",
        "-o",
        wsl_output,
        wsl_input,
    ]

    return run_command(
        tool_name=f"AlcoR Mode {mode}",
        command=cmd,
        timeout=timeout,
    )


# ==========================================================
# XSTREAM
# ==========================================================
def run_xstream(
    input_fasta: Path,
    output_dir: Path,
    mode: int = 1,
    timeout: int = DEFAULT_TIMEOUT,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output_dir = convert_to_wsl_path(output_dir)

    command = (
        f"cd {shlex.quote(wsl_output_dir)} && "
        f"java -jar {XSTREAM} "
        f"{shlex.quote(wsl_input)} "
        f"-m{mode} -a_m{mode}"
    )

    return run_bash_command(
        tool_name=f"XSTREAM m{mode}",
        bash_command=command,
        timeout=timeout,
    )
# ==========================================================
# T-REKS
# ==========================================================

def run_treks(
    input_fasta: Path,
    output_dir: Path,
    timeout: int = DEFAULT_TIMEOUT,
):
    """
    Run T-REKS.

    Alignment file is NOT generated.
    """

    output_file = build_output(
        output_dir,
        input_fasta,
        "treks.tsv",
    )

    wsl_input = convert_to_wsl_path(input_fasta)
    wsl_output = convert_to_wsl_path(output_file)

    cmd = [
        "wsl",
        "java",
        "-Xmx4G",
        "-jar",
        TREKS,
        "-f",
        wsl_input,
        "-t",
        wsl_output,
        "-c",
        CLUSTALW,
    ]

    return run_command(
        tool_name="T-REKS",
        command=cmd,
        timeout=timeout,
    )