from pathlib import Path
from PIL import Image
import math
import matplotlib.pyplot as plt
import sys
# =====================================================
# Directories
# =====================================================

INPUT_DIR = Path(sys.argv[1])  # Directory containing individual heatmaps   
OUTPUT_DIR = Path(sys.argv[2])  # Directory to save the combined figure

OUTPUT_DIR.parent.mkdir(
    parents=True,
    exist_ok=True
)

# =====================================================
# Read all PNG heatmaps
# =====================================================

images = sorted(
    [
        f for f in INPUT_DIR.glob("*.png")
        if "combined" not in f.stem.lower()
    ]
)

if len(images) == 0:
    raise RuntimeError("No heatmaps found.")

# =====================================================
# Determine layout automatically
# =====================================================

n = len(images)

cols = min(4, math.ceil(math.sqrt(n)))

rows = math.ceil(n / cols)

# =====================================================
# Create figure
# =====================================================

fig, axes = plt.subplots(
    rows,
    cols,
    figsize=(5*cols,5*rows)
)

# Flatten axes safely

if rows == 1 and cols == 1:
    axes = [axes]
elif rows == 1 or cols == 1:
    axes = axes.flatten()
else:
    axes = axes.ravel()

# =====================================================
# Add each image
# =====================================================

for ax, img_file in zip(axes, images):

    img = Image.open(img_file)

    ax.imshow(img)

    ax.set_title(
        img_file.stem.replace("_Fig5",""),
        fontsize=12
    )

    ax.axis("off")

# Hide unused panels

for ax in axes[len(images):]:
    ax.axis("off")

plt.suptitle(
    "Relationship between amino acid composition and mutation percentage",
    fontsize=18
)

plt.tight_layout()

# =====================================================
# Save
# =====================================================

plt.savefig(
    OUTPUT_DIR,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("Combined figure saved.")