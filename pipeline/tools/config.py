"""
config.py

Central configuration file for the LCR benchmarking pipeline.

Edit this file only if the location of the external tools changes.
"""

# ==========================================================
# Tool Locations (WSL)
# ==========================================================

SEG = "segmasker"

FLPS = (
    "~/Lcr_Project/tools/flps/flps/fLPS/bin/linux/fLPS"
)

FLPS2 = (
    "~/Lcr_Project/tools/flps2/fLPS2programs/src/fLPS2"
)

ALCOR = (
    "~/Lcr_Project/tools/alcor/alcor/bin/AlcoR"
)

XSTREAM = (
    "~/Lcr_Project/tools/XSTREAM/xstream.jar"
)

TREKS = (
    "~/Lcr_Project/tools/treks-hpc/target/T-ReksHPC_0.1-SNAPSHOT.jar"
)

CLUSTALW = "/usr/bin/clustalw"

# ==========================================================
# Pipeline Settings
# ==========================================================

DEFAULT_TIMEOUT = 3600