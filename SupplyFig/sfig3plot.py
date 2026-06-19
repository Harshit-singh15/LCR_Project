import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

all_df=[]

for k in range(1,14):

    df=pd.read_csv(
        f"SupplyFig\\complexity\\complexity_{k}.tsv",
        sep="\t"
    )

    df["k"]=k

    all_df.append(df)

df=pd.concat(all_df)

df["xbin"]=(df["Mutation_Percent"]//2).astype(int)
df["ybin"]=(df["Most_Frequent_AA_Percent"]//2).astype(int)

heat=(
    df
    .groupby(["xbin","ybin"])
    ["k"]
    .max()
    .reset_index()
)

matrix=np.full((51,51),np.nan)

for _,row in heat.iterrows():

    matrix[
        int(row["ybin"]),
        int(row["xbin"])
    ]=row["k"]

plt.figure(figsize=(8,8))

plt.imshow(
    matrix,
    origin="lower",
    aspect="auto",
    extent=[0,100,0,100]
)

plt.colorbar(label="Maximum Consensus k")

plt.xlabel("Mutation Percent")
plt.ylabel("Most Frequent AA Percent")

plt.title("Suppl Figure 3 - Consensus Map")

plt.savefig(
    "SupplyFig\\plots\\Suppl_Fig3_consensus_map.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()