import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
import re
import os

# ===================================
# EDIT
# ===================================

COMPLEXITY_DIR = "SupplyFig\\complexity"

BACKGROUND_FILE = "Fig8\\fig8A_disprot_complexity.tsv"

OUTPUT_DIR = "SupplyFig\\SupplFig5"

# ===================================

BIN_SIZE = 2
ALPHA = 1

KMINS = [7,8,9,10,11,12]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ===================================
# LOAD ALL COMPLEXITY FILES
# ===================================

dfs = []

for f in Path(COMPLEXITY_DIR).glob("complexity_*.tsv"):

    k = int(
        re.search(r"(\d+)", f.stem).group(1)
    )

    tmp = pd.read_csv(
        f,
        sep="\t"
    )

    tmp["k"] = k

    dfs.append(tmp)

lcr = pd.concat(
    dfs,
    ignore_index=True
)

# ===================================
# BACKGROUND
# ===================================

bg = pd.read_csv(
    BACKGROUND_FILE,
    sep="\t"
)

bg = bg[
    [
        "Mutation_Percent",
        "Most_Frequent_AA_Percent"
    ]
].copy()

# ===================================
# BIN HELPERS
# ===================================

def add_bins(df):

    df = df.copy()

    df["xb"] = (
        np.floor(
            df["Mutation_Percent"]/BIN_SIZE
        ) * BIN_SIZE
    )

    df["yb"] = (
        np.floor(
            df["Most_Frequent_AA_Percent"]/BIN_SIZE
        ) * BIN_SIZE
    )

    return df


def build_posterior_grid(
    pos_train,
    neg_train
):

    pos_train = add_bins(pos_train)
    neg_train = add_bins(neg_train)

    pos_counts = (
        pos_train
        .groupby(["xb","yb"])
        .size()
        .reset_index(name="n_pos")
    )

    neg_counts = (
        neg_train
        .groupby(["xb","yb"])
        .size()
        .reset_index(name="n_neg")
    )

    grid = pd.merge(
        pos_counts,
        neg_counts,
        how="outer",
        on=["xb","yb"]
    ).fillna(0)

    grid["posterior"] = (
        grid["n_pos"] + ALPHA
    ) / (
        grid["n_pos"] +
        grid["n_neg"] +
        2*ALPHA
    )

    return grid


def score_points(
    df,
    score_lookup
):

    df = add_bins(df)

    return np.array(
        [
            score_lookup.get(
                (x,y),
                0.5
            )
            for x,y in zip(
                df["xb"],
                df["yb"]
            )
        ]
    )


# ===================================
# COMBINED FIGURE
# ===================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(16,10)
)

axes = axes.flatten()

letters = list("ABCDEF")

last_cf = None

# ===================================
# MAIN LOOP
# ===================================

for idx, kmin in enumerate(KMINS):

    print(f"k >= {kmin}")

    pos = lcr[
        lcr["k"] >= kmin
    ][
        [
            "Mutation_Percent",
            "Most_Frequent_AA_Percent"
        ]
    ].copy()

    pos_train, pos_test = train_test_split(
        pos,
        test_size=0.25,
        random_state=1
    )

    neg_train, neg_test = train_test_split(
        bg,
        test_size=0.25,
        random_state=1
    )

    # ------------------------
    # Posterior grid
    # ------------------------

    grid = build_posterior_grid(
        pos_train,
        neg_train
    )

    lookup = {
        (r.xb, r.yb): r.posterior
        for r in grid.itertuples()
    }

    # ------------------------
    # ROC
    # ------------------------

    pos_scores = score_points(
        pos_test,
        lookup
    )

    neg_scores = score_points(
        neg_test,
        lookup
    )

    y_true = np.concatenate([
        np.ones(len(pos_scores)),
        np.zeros(len(neg_scores))
    ])

    y_score = np.concatenate([
        pos_scores,
        neg_scores
    ])

    fpr, tpr, thr = roc_curve(
        y_true,
        y_score
    )

    auc_val = auc(
        fpr,
        tpr
    )

    J = tpr - fpr

    best_idx = np.argmax(J)

    p_thresh = thr[best_idx]

    # ------------------------
    # Full 51x51 surface
    # ------------------------

    xvals = np.arange(
        0,
        102,
        BIN_SIZE
    )

    yvals = np.arange(
        0,
        102,
        BIN_SIZE
    )

    surface = np.zeros(
        (
            len(yvals),
            len(xvals)
        )
    )

    surface[:] = 0.5

    for row in grid.itertuples():

        xi = int(
            row.xb / BIN_SIZE
        )

        yi = int(
            row.yb / BIN_SIZE
        )

        surface[yi,xi] = row.posterior

    X,Y = np.meshgrid(
        xvals,
        yvals
    )

    # ==================================
    # INDIVIDUAL PANEL
    # ==================================

    fig_single, ax_single = plt.subplots(
        figsize=(6,5)
    )

    cf = ax_single.contourf(
        X,
        Y,
        surface,
        levels=np.linspace(0,1,100),
        cmap="viridis"
    )

    ax_single.contour(
        X,
        Y,
        surface,
        levels=[p_thresh],
        colors="black",
        linewidths=2
    )

    ax_single.set_xlim(0,100)
    ax_single.set_ylim(0,100)

    ax_single.set_xlabel(
        "Mutation (%)"
    )

    ax_single.set_ylabel(
        "Most frequent amino acid (%)"
    )

    ax_single.set_title(
        f"k ≥ {kmin}\n"
        f"AUC={auc_val:.3f}  "
        f"p={p_thresh:.3f}"
    )

    plt.colorbar(
        cf,
        ax=ax_single,
        label="P(LCR|bin)"
    )

    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/k{kmin}.png",
        dpi=300
    )

    plt.close()

    # ==================================
    # COMBINED PANEL
    # ==================================

    ax = axes[idx]

    last_cf = ax.contourf(
        X,
        Y,
        surface,
        levels=np.linspace(0,1,100),
        cmap="viridis"
    )

    ax.contour(
        X,
        Y,
        surface,
        levels=[p_thresh],
        colors="black",
        linewidths=2
    )

    ax.set_xlim(0,100)
    ax.set_ylim(0,100)

    ax.set_title(
        f"{letters[idx]}) k≥{kmin}\n"
        f"AUC={auc_val:.3f}"
    )

# ===================================
# COLORBAR
# ===================================

cbar = fig.colorbar(
    last_cf,
    ax=axes,
    shrink=0.8
)

cbar.set_label(
    "P(LCR|bin)"
)

plt.tight_layout()

plt.savefig(
    f"{OUTPUT_DIR}/SupplFig5_combined.png",
    dpi=300
)

plt.show()