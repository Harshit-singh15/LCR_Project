import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==========================================================
# INPUT / OUTPUT
# ==========================================================

input_file = r"fruitfly\dataforFig8\fig8A_disprot_complexity.tsv"

output_dir = Path(r"fruitfly\Fig_outputs\Fig8")
output_dir.mkdir(
    parents=True,
    exist_ok=True
)

# ==========================================================

df = pd.read_csv(
    input_file,
    sep="\t"
)


# ==========================================================
# GROUP EXPERIMENTAL METHODS
# ==========================================================

def simplify_method(ec):

    if pd.isna(ec):
        return "Other"

    ec = str(ec).lower()

    # ------------------------------------------------------
    # NMR spectroscopy
    # ------------------------------------------------------

    if any(x in ec for x in [

        "nuclear magnetic resonance",

        "heteronuclear single quantum coherence",

        "proton-based nuclear magnetic resonance"

    ]):

        return "NMR spectroscopy"

    # ------------------------------------------------------
    # CD / IR
    # ------------------------------------------------------

    if any(x in ec for x in [

        "circular dichroism",

        "infrared spectroscopy"

    ]):

        return "CD / IR spectroscopy"

    # ------------------------------------------------------
    # SAXS / SANS
    # ------------------------------------------------------

    if any(x in ec for x in [

        "small-angle x-ray scattering",

        "small-angle neutron scattering"

    ]):

        return "Small-angle scattering"

    # ------------------------------------------------------

    if "x-ray crystallography" in ec:

        return "X-ray crystallography"

    # ------------------------------------------------------

    if "cryogenic electron microscopy" in ec:

        return "Cryo-EM"

    # ------------------------------------------------------

    if "mass spectrometry" in ec:

        return "Mass spectrometry"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "microscopy",

        "confocal"

    ]):

        return "Microscopy"

    # ------------------------------------------------------

    if "fluorescence" in ec:

        return "Fluorescence methods"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "chromatography",

        "gel-filtration",

        "electrophoresis",

        "electrophoretic mobility",

        "gel electrophoresis",

        "sodium dodecyl sulfate"

    ]):

        return "Chromatography / Electrophoresis"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "protein binding",

        "co-immunoprecipitation",

        "immunoprecipitation",

        "pull-down assay",

        "yeast 2-hybrid"

    ]):

        return "Interaction assays"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "mutant phenotype",

        "loss-of-function",

        "physiological response"

    ]):

        return "Genetic / Transcriptional assays"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "protein kinase assay",

        "enzyme-linked immunoabsorbent",

        "immunodetection",

        "immunohistochemistry",

        "western immunoblotting",

        "blocking monoclonal antibody"

    ]):

        return "Cell-based / Phenotypic assays"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "calorimetry",

        "thermal shift",

        "temperature-induced",

        "light scattering"

    ]):

        return "Thermodynamics"

    # ------------------------------------------------------

    if any(x in ec for x in [

        "cleavage assay",

        "deglycosylation assay"

    ]):

        return "Biochemical activity assays"

    # ------------------------------------------------------

    if "electron paramagnetic resonance" in ec:

        return "Electron paramagnetic resonance"

    return "Other"


df["Method_Group"] = df["Experimental_Method"].apply(
    simplify_method
)

# ==========================================================
# PLOT
# ==========================================================

plt.figure(figsize=(10, 8))

methods = sorted(
    df["Method_Group"].unique()
)

for method in methods:

    temp = df[
        df["Method_Group"] == method
    ]

    plt.scatter(

        temp["Mutation_Percent"],

        temp["Most_Frequent_AA_Percent"],

        s=30,

        alpha=0.75,

        label=method

    )

plt.xlabel(
    "Mutation Percentage (%)",
    fontsize=12
)

plt.ylabel(
    "Most Frequent Amino Acid Percentage (%)",
    fontsize=12
)

plt.title(
    "Fig. 8A  DisProt Intrinsically Disordered Regions",
    fontsize=14
)

plt.xlim(0, 100)
plt.ylim(0, 100)

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.grid(
    alpha=0.3
)

plt.legend(

    fontsize=8,

    frameon=False,

    markerscale=1.3,

    bbox_to_anchor=(1.02, 1),

    loc="upper left"

)

plt.tight_layout()

plt.savefig(
    output_dir / "Fig8A_disprot_complexity.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()