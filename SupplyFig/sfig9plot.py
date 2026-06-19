import pandas as pd
import matplotlib.pyplot as plt

# =====================================
# INPUT / OUTPUT
# =====================================

INPUT_FILE = r"SupplyFig\SupplFig9\retention_table.tsv"

OUTPUT_FILE = r"SupplyFig\SupplFig9\Suppl_Fig9.png"

# =====================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print("Rows:", len(df))
print("\nColumns:")
print(df.columns.tolist())

print("\nTools found:")
print(df["Tool"].unique())

# =====================================
# COLORS
# =====================================

palette = [
    "#5a5050",   # dark grey
    "#d8d8d8",   # light grey
    "#ff0022",   # red
    "#e600ff",   # magenta
    "#20e020",   # green
    "#0066ff",   # blue
    "#ff9900",   # orange
    "#00cccc",   # cyan
]

tools = sorted(df["Tool"].unique())

colors = {
    tool: palette[i % len(palette)]
    for i, tool in enumerate(tools)
}

# =====================================
# PLOT
# =====================================

plt.figure(figsize=(10, 8))

for tool in tools:

    sub = (
        df[df["Tool"] == tool]
        .sort_values("Purity")
    )

    print(
        f"{tool}: {len(sub)} points"
    )

    plt.plot(
        sub["Purity"],
        sub["Proportion"],
        marker="o",
        linewidth=3,
        markersize=8,
        color=colors[tool],
        label=tool
    )

# =====================================
# FORMATTING
# =====================================

plt.title(
    "Effect of Purity on Proportion of Retained Entities",
    fontsize=22,
    pad=20
)

plt.xlabel(
    "Purity Level (%)",
    fontsize=16
)

plt.ylabel(
    "Proportion Retained",
    fontsize=16
)

plt.xlim(0, 90)
plt.ylim(0, 1.05)

plt.xticks(
    range(0, 100, 10),
    fontsize=12
)

plt.yticks(
    fontsize=12
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend(
    title="Dataset",
    fontsize=12,
    title_fontsize=13
)

plt.tight_layout()

# =====================================
# SAVE
# =====================================

plt.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

print(
    "\nSaved:",
    OUTPUT_FILE
)

plt.show()