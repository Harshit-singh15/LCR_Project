import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. DATA LOADING AND PROCESSING
# ==========================================
summary = []

for k in range(1, 14):
    df = pd.read_csv(
        f"SupplyFig\\complexity\\complexity_{k}.tsv",
        sep="\t"
    )

    x = df["Mutation_Percent"]
    y = df["Most_Frequent_AA_Percent"]

    summary.append({
        "k": k,
        "x_med": x.median(),
        "y_med": y.median(),
        "x_q1": x.quantile(0.25),
        "x_q3": x.quantile(0.75),
        "y_q1": y.quantile(0.25),
        "y_q3": y.quantile(0.75)
    })

# Convert the list of dictionaries into a clean DataFrame
summary = pd.DataFrame(summary)

# ==========================================
# 2. PLOTTING SECTION (WITH LEGEND)
# ==========================================
# Generate a distinct color palette for the 13 tiers using Set2
colors = plt.cm.Set2(np.linspace(0, 1, len(summary)))

# Slightly widened figure width (11 instead of 10) to make room for the legend on the right
plt.figure(figsize=(11, 8))

# Draw the background grey trajectory line connecting the medians
plt.plot(
    summary["x_med"],
    summary["y_med"],
    color="lightgrey",
    linewidth=2,
    zorder=1
)

# Loop through each row to plot individual IQR bars, points, and labels
for i, row in summary.iterrows():
    c = colors[i]
    k_val = int(row["k"])

    # Horizontal IQR bar
    plt.plot(
        [row["x_q1"], row["x_q3"]],
        [row["y_med"], row["y_med"]],
        color=c,
        linewidth=2
    )

    # Vertical IQR bar
    plt.plot(
        [row["x_med"], row["x_med"]],
        [row["y_q1"], row["y_q3"]],
        color=c,
        linewidth=2
    )

    # Centroid point (Includes 'label' so it shows up cleanly in the legend)
    plt.scatter(
        row["x_med"],
        row["y_med"],
        s=180,
        color=c,
        edgecolor="black",
        label=f"k = {k_val}",
        zorder=3
    )

    # Offset text label placed just above the point
    plt.text(
        row["x_med"],
        row["y_med"] + 1.5,
        str(k_val),
        color=c,
        fontsize=11,
        ha="center",
        weight="bold"
    )

# Axis boundaries
plt.xlim(0, 100)
plt.ylim(0, 100)

# Labels and Styling
plt.xlabel("Mutation (%)", fontsize=16)
plt.ylabel("Most frequent amino acid (%)", fontsize=16)

plt.title(
    "Centroid flow across consensus tiers (median and IQR)",
    fontsize=18,
    pad=15
)

# Position the legend neatly outside the plot grid on the upper right side
plt.legend(
    title="Consensus Tiers (k)", 
    bbox_to_anchor=(1.02, 1), 
    loc='upper left', 
    borderaxespad=0,
    fontsize=11,
    title_fontsize=12
)

plt.tight_layout()

# Save image (bbox_inches="tight" ensures the legend doesn't get clipped)
plt.savefig(
    "SupplyFig\\plots\\Suppl_Fig2_centroid_flow.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()