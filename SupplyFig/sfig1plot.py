import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pathlib import Path

plt.figure(figsize=(10,8))

colors = plt.cm.viridis(np.linspace(0,1,13))

for k in range(1,14):

    file = f"SupplyFig\\complexity\\complexity_{k}.tsv"

    df = pd.read_csv(file, sep="\t")

    x = df["Mutation_Percent"].values
    y = df["Most_Frequent_AA_Percent"].values

    xy = np.vstack([x,y])

    kde = gaussian_kde(xy)

    xx,yy = np.mgrid[0:100:100j,0:100:100j]

    positions = np.vstack([xx.ravel(),yy.ravel()])

    zz = kde(positions).reshape(xx.shape)

    plt.contour(
        xx,
        yy,
        zz,
        levels=5,
        colors=[colors[k-1]]
    )

plt.xlim(0,100)
plt.ylim(0,100)

plt.xlabel("Mutation Percent")
plt.ylabel("Most Frequent AA Percent")

plt.title("Suppl Figure 1")

plt.savefig(
    "SupplyFig\\plots\\Suppl_Fig1_density_contours.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()