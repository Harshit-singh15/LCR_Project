import pandas as pd
from plot_fig8A import simplify_method
# =====================================
# INPUT
# =====================================

input_file = r"celegans\dataforFig8\fig8A_disprot_complexity.tsv"

# =====================================

df = pd.read_csv(
    input_file,
    sep="\t"
)

df["Method_Group"] = df["Experimental_Method"].apply(
    simplify_method
)

print("\nGrouped Experimental Methods")
print("=" * 60)

counts = (
    df["Method_Group"]
    .value_counts()
)

percent = (
    counts / len(df) * 100
).round(2)

summary = pd.DataFrame({
    "Count": counts,
    "Percent": percent
})

print(summary)

print("\n")

print("=" * 60)
print("Original Experimental Methods")
print("=" * 60)

print(
    df["Experimental_Method"]
    .value_counts()
)

print(len(df))
print(df["Protein_ID"].nunique())
print(
df[
[
"Protein_ID",
"Region_ID",
"Mutation_Percent",
"Most_Frequent_AA_Percent"
]
].duplicated().sum()
)