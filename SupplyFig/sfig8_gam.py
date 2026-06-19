import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from pygam import LogisticGAM, te
import re
import os

# =====================================
# EDIT THESE
# =====================================

COMPLEXITY_DIR = "SupplyFig\\complexity"

BACKGROUND_FILE = "Fig8\\fig8B_missing_residue_complexity.tsv"

OUTPUT_DIR = "SupplyFig\\SupplFig8"

# =====================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

KMINS = [7,8,9,10,11,12]

# =====================================
# LOAD ALL LCR FILES
# =====================================

dfs = []

for f in Path(COMPLEXITY_DIR).glob("complexity_*.tsv"):

    k = int(
        re.search(
            r"(\d+)",
            f.stem
        ).group(1)
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

# =====================================
# LOAD BACKGROUND
# =====================================

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

# =====================================
# COMBINED FIGURE
# =====================================

fig, axes = plt.subplots(
    2,
    3,
    figsize=(16,10)
)

axes = axes.flatten()

letters = list("ABCDEF")

last_cf = None

# =====================================
# LOOP OVER KMIN
# =====================================

for idx, kmin in enumerate(KMINS):

    print(f"\nProcessing k >= {kmin}")

    pos = lcr[
        lcr["k"] >= kmin
    ][
        [
            "Mutation_Percent",
            "Most_Frequent_AA_Percent"
        ]
    ].copy()

    pos["label"] = 1

    neg = bg.copy()

    neg["label"] = 0

    data = pd.concat(
        [pos, neg],
        ignore_index=True
    )

    train_df, test_df = train_test_split(
        data,
        test_size=0.25,
        random_state=1,
        stratify=data["label"]
    )

    X_train = train_df[
        [
            "Mutation_Percent",
            "Most_Frequent_AA_Percent"
        ]
    ].values

    y_train = train_df["label"].values

    X_test = test_df[
        [
            "Mutation_Percent",
            "Most_Frequent_AA_Percent"
        ]
    ].values

    y_test = test_df["label"].values

    # =====================================
    # FIT GAM
    # =====================================

    gam = LogisticGAM(
        te(
            0,
            1,
            n_splines=20
        )
    )

    gam.fit(
        X_train,
        y_train
    )

    # =====================================
    # ROC
    # =====================================

    test_probs = gam.predict_proba(
        X_test
    )

    fpr, tpr, thr = roc_curve(
        y_test,
        test_probs
    )

    auc_val = auc(
        fpr,
        tpr
    )

    J = tpr - fpr

    best_idx = np.argmax(J)

    p_thresh = thr[best_idx]

    # =====================================
    # PREDICTION GRID
    # =====================================

    xgrid = np.linspace(
        0,
        100,
        201
    )

    ygrid = np.linspace(
        0,
        100,
        201
    )

    Xg, Yg = np.meshgrid(
        xgrid,
        ygrid
    )

    grid_points = np.column_stack(
        [
            Xg.ravel(),
            Yg.ravel()
        ]
    )

    surface = gam.predict_proba(
        grid_points
    )

    surface = surface.reshape(
        Xg.shape
    )

    # =====================================
    # INDIVIDUAL PANEL
    # =====================================

    fig_single, ax_single = plt.subplots(
        figsize=(6,5)
    )

    cf = ax_single.contourf(
        Xg,
        Yg,
        surface,
        levels=100,
        cmap="viridis"
    )

    ax_single.contour(
        Xg,
        Yg,
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
        f"AUC={auc_val:.3f}\n"
        f"p={p_thresh:.3f}"
    )

    plt.colorbar(
        cf,
        ax=ax_single,
        label="P(LCR)"
    )

    plt.tight_layout()

    plt.savefig(
        f"{OUTPUT_DIR}/k{kmin}.png",
        dpi=300
    )

    plt.close()

    # =====================================
    # COMBINED PANEL
    # =====================================

    ax = axes[idx]

    last_cf = ax.contourf(
        Xg,
        Yg,
        surface,
        levels=100,
        cmap="viridis"
    )

    ax.contour(
        Xg,
        Yg,
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

# =====================================
# COLORBAR
# =====================================

cbar = fig.colorbar(
    last_cf,
    ax=axes,
    shrink=0.8
)

cbar.set_label(
    "P(LCR)"
)

plt.tight_layout()

plt.savefig(
    f"{OUTPUT_DIR}/SupplFig8_combined.png",
    dpi=300
)

plt.show()