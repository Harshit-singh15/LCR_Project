import pandas as pd
import matplotlib.pyplot as plt

input_file = r"Fig8\fig8A_disprot_complexity.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

def simplify_method(ec):

    ec = str(ec).lower()

    # =================================================
    # NMR spectroscopy
    # =================================================
    if any(x in ec for x in [
        "nuclear magnetic resonance",
        "heteronuclear single quantum coherence",
        "proton-based nuclear magnetic resonance"
    ]):
        return "NMR spectroscopy"

    # =================================================
    # CD / IR spectroscopy
    # =================================================
    if any(x in ec for x in [
        "circular dichroism",
        "infrared spectroscopy"
    ]):
        return "CD/IR spectroscopy"

    # =================================================
    # Small-angle scattering
    # =================================================
    if any(x in ec for x in [
        "small-angle x-ray scattering",
        "small-angle neutron scattering"
    ]):
        return "Small-angle scattering"

    # =================================================
    # X-ray crystallography
    # =================================================
    if "x-ray crystallography" in ec:
        return "X-ray crystallography"

    # =================================================
    # Cryo-EM
    # =================================================
    if "cryogenic electron microscopy" in ec:
        return "Cryo-EM"

    # =================================================
    # Mass spectrometry
    # =================================================
    if "mass spectrometry" in ec:
        return "Mass spectrometry"

    # =================================================
    # Microscopy
    # =================================================
    if any(x in ec for x in [
        "microscopy",
        "confocal"
    ]):
        return "Microscopy"

    # =================================================
    # Fluorescence methods
    # =================================================
    if "fluorescence" in ec:
        return "Fluorescence methods"

    # =================================================
    # Chromatography / electrophoresis
    # =================================================
    if any(x in ec for x in [
        "chromatography",
        "gel-filtration",
        "electrophoresis",
        "electrophoretic mobility",
        "gel electrophoresis",
        "sodium dodecyl sulfate"
    ]):
        return "Chromatography/electrophoresis"

    # =================================================
    # Interaction assays
    # =================================================
    if any(x in ec for x in [
        "protein binding",
        "co-immunoprecipitation",
        "immunoprecipitation",
        "pull-down assay",
        "yeast 2-hybrid"
    ]):
        return "Interaction assays"

    # =================================================
    # Genetic / transcriptional assays
    # =================================================
    if any(x in ec for x in [
        "mutant phenotype",
        "loss-of-function",
        "physiological response"
    ]):
        return "Genetic/transcriptional assays"

    # =================================================
    # Cell-based / phenotypic assays
    # =================================================
    if any(x in ec for x in [
        "protein kinase assay",
        "enzyme-linked immunoabsorbent",
        "immunodetection",
        "immunohistochemistry",
        "western immunoblotting",
        "blocking monoclonal antibody"
    ]):
        return "Cell-based/phenotypic assays"

    # =================================================
    # Thermodynamics
    # =================================================
    if any(x in ec for x in [
        "calorimetry",
        "thermal shift",
        "temperature-induced",
        "light scattering"
    ]):
        return "Thermodynamics"

    # =================================================
    # Biochemical activity assays
    # =================================================
    if any(x in ec for x in [
        "cleavage assay",
        "deglycosylation assay"
    ]):
        return "Biochemical activity assays"

    # =================================================
    # Electron microscopy
    # =================================================
    if "electron paramagnetic resonance" in ec:
        return "Electron microscopy"

    return "Other"

df["Method_Group"] = df["ec_name"].apply(
    simplify_method
)

for method in sorted(df["Method_Group"].unique()):

    temp = df[
        df["Method_Group"] == method
    ]

    plt.scatter(
        temp["Mutation_Percent"],
        temp["Most_Frequent_AA_Percent"],
        label=method,
        s=35,
        alpha=0.8
    )
plt.figure(
    figsize=(10,8)
)

methods = df["ec_name"].dropna().unique()

for method in methods:

    temp = df[
        df["ec_name"] == method
    ]

    plt.scatter(
        temp["Mutation_Percent"],
        temp["Most_Frequent_AA_Percent"],
        s=15,
        alpha=0.6,
        label=method
    )

plt.xlabel(
    "Mutation Percentage"
)

plt.ylabel(
    "Most Frequent Amino Acid Percentage"
)

plt.title(
    "Fig. 8A - DisProt Intrinsically Disordered Regions"
)

plt.xlim(0,100)
plt.ylim(0,100)

plt.grid(alpha=0.3)

plt.legend(
    fontsize=6,
    bbox_to_anchor=(1.05,1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    "Figure8A.png",
    dpi=300
)

plt.show()