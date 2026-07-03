#!/usr/bin/env python3
"""
=========================================================
Generate Benchmark_Report.pdf
=========================================================

Reads:
    BASE/
        Report/report.yaml
        Fig_outputs/

Produces:
    BASE/
        Benchmark_Report.pdf

Usage:
    python Report/generate_report.py <BASE>

"""

import sys
from pathlib import Path
import yaml

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    PageBreak,
)
from reportlab.pdfbase import pdfmetrics


# =========================================================
# Base directory
# =========================================================

BASE = Path(sys.argv[1]).resolve()

REPORT = BASE / "reports"
FIG = BASE / "Fig_outputs"

YAML_FILE = REPORT / "captions.yaml"
OUTPUT = BASE / "Benchmark_Report.pdf"


# =========================================================
# Load YAML
# =========================================================

with open(YAML_FILE, "r", encoding="utf-8") as f:
    CAP = yaml.safe_load(f)


# =========================================================
# Styles
# =========================================================

styles = getSampleStyleSheet()

title_style = styles["Heading1"]
title_style.alignment = TA_CENTER
title_style.spaceAfter = 12
title_style.fontSize = 22

fig_style = styles["Heading2"]
fig_style.alignment = TA_CENTER
fig_style.spaceBefore = 6
fig_style.spaceAfter = 10
fig_style.fontSize = 16

caption_style = styles["Italic"]
caption_style.fontSize = 10
caption_style.leading = 15
caption_style.spaceBefore = 8
caption_style.spaceAfter = 6

sub_style = styles["BodyText"]
sub_style.fontSize = 10
sub_style.leading = 15
sub_style.leftIndent = 18
sub_style.spaceAfter = 4


# =========================================================
# Footer
# =========================================================

def add_page_number(canvas, doc):
    page = canvas.getPageNumber()

    canvas.setFont("Helvetica", 9)

    canvas.drawCentredString(
        A4[0] / 2,
        0.4 * inch,
        f"Page {page}"
    )


# =========================================================
# Figure helper
# =========================================================

story = []


def add_figure(key, image_path):

    if not image_path.exists():
        print(f"Skipping missing image: {image_path}")
        return

    meta = CAP.get(key, {})

    title = meta.get(
        "title",
        key.replace("figure", "Figure ").upper()
    )

    caption = meta.get("caption", "")

    story.append(Paragraph(title, fig_style))
    story.append(Spacer(1, 0.12 * inch))

    img = Image(str(image_path))

    max_width = 6.8 * inch
    max_height = 8.0 * inch

    scale = min(
        max_width / img.drawWidth,
        max_height / img.drawHeight
    )

    img.drawWidth *= scale
    img.drawHeight *= scale

    story.append(img)

    story.append(Spacer(1, 0.12 * inch))

    if caption:

        story.append(
            Paragraph(
                f"<i>{caption}</i>",
                caption_style
            )
        )

    # -----------------------------------------------------
    # Subfigure captions
    # -----------------------------------------------------

# =========================================================
# Cover page
# =========================================================

story.append(
    Paragraph(
        "Benchmarking of Low-Complexity Region Detection Tools",
        title_style,
    )
)

story.append(Spacer(1, 1.0 * inch))

story.append(
    Paragraph(
        "<b>Automatically Generated Benchmark Report</b>",
        styles["Heading2"],
    )
)

story.append(Spacer(1, 0.25 * inch))

story.append(
    Paragraph(
        f"<b>Dataset:</b> {BASE.name}",
        styles["Normal"],
    )
)

story.append(
    Paragraph(
        f"<b>Project Folder:</b> {BASE}",
        styles["Normal"],
    )
)

story.append(PageBreak())


# =========================================================
# Individual Figures
# =========================================================

figures = [

    ("figure1", FIG / "Fig1" / "Figure1A_Length_Heatmap.png"),
    ("figure1", FIG / "Fig1" / "Figure1B_Coverage_Heatmap.png"),
    ("figure1", FIG / "Fig1" / "Figure1C_CountDistribution.png"),
    ("figure1", FIG / "Fig1" / "Figure1D_AminoAcidComposition.png"),
    ("figure1", FIG / "Fig1" / "Figure1E_Entropy_Boxplot.png"),

    ("figure2a", FIG / "Fig2" / "Fig2A_PeptideMotifs.png"),
    ("figure2b", FIG / "Fig2" / "Fig2B_Entropy.png"),
    ("figure2c", FIG / "Fig2" / "Fig2C_Purity.png"),

    ("figure3", FIG / "Fig3" / "Fig3_Purity.png"),

    ("figure4", FIG / "Fig4" / "Fig4_Jaccard.png"),

    ("figure5", FIG / "Fig5" / "Fig5_Combined.png"),

    ("figure6a", FIG / "Fig6" / "Fig6A_GeneLength.png"),
    ("figure6b", FIG / "Fig6" / "Fig6B_LCRCount.png"),
    ("figure6c", FIG / "Fig6" / "Fig6C_Coverage.png"),
    ("figure6d", FIG / "Fig6" / "Fig6D_EntropyRatio.png"),

    ("figure7", FIG / "Fig7" / "Fig7_summary.png"),
]

for key, img in figures:
    add_figure(key, img)


# =========================================================
# Build PDF
# =========================================================

doc = SimpleDocTemplate(
    str(OUTPUT),
    pagesize=A4,
    rightMargin=0.55 * inch,
    leftMargin=0.55 * inch,
    topMargin=0.6 * inch,
    bottomMargin=0.6 * inch,
)

doc.build(
    story,
    onFirstPage=add_page_number,
    onLaterPages=add_page_number,
)

print("\n====================================")
print("Benchmark report generated")
print(OUTPUT)
print("====================================")