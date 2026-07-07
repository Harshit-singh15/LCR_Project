#!/usr/bin/env python3
"""
=========================================================
Generate Benchmark_Report.pdf
=========================================================

Reads:
    BASE/
        Fig_outputs/
    REPORT/
        captions.yaml

Produces:
    OUTPUT (Benchmark_Report.pdf)

Usage:
    python Report/generate_report.py <BASE> <OUTPUT>

"""

import sys
import datetime
from pathlib import Path

import yaml

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    PageBreak,
    Table,
    TableStyle,
    HRFlowable,
    KeepTogether,
)
from reportlab.pdfgen import canvas as canvas_module

# =========================================================
# Paths
# =========================================================

BASE = Path(sys.argv[1]).resolve()
OUTPUT = Path(sys.argv[2]).resolve()

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "reports"
FIG = BASE / "Fig_outputs"
YAML_FILE = REPORT / "captions.yaml"

with open(YAML_FILE, "r", encoding="utf-8") as f:
    CAP = yaml.safe_load(f)

# =========================================================
# Palette
# =========================================================

NAVY = HexColor("#1F3A54")     # headings / section titles
TEAL = HexColor("#3E7C8C")     # accents / rules / sub-labels
INK = HexColor("#2B2B2B")      # body text
SLATE = HexColor("#5C6773")    # captions, secondary text
PALE = HexColor("#EEF2F5")     # cover info box fill
LINE = HexColor("#C9D3D8")     # thin rules

REPORT_TITLE = "Benchmarking of Low-Complexity Region Detection Tools"

# =========================================================
# Styles
# =========================================================

base = getSampleStyleSheet()

cover_title_style = ParagraphStyle(
    "CoverTitle", parent=base["Title"],
    fontName="Helvetica-Bold", fontSize=25, leading=31,
    textColor=NAVY, alignment=TA_CENTER, spaceAfter=6,
)

cover_subtitle_style = ParagraphStyle(
    "CoverSubtitle", parent=base["Normal"],
    fontName="Helvetica", fontSize=13, leading=18,
    textColor=TEAL, alignment=TA_CENTER, spaceAfter=0,
)

cover_meta_label_style = ParagraphStyle(
    "CoverMetaLabel", parent=base["Normal"],
    fontName="Helvetica-Bold", fontSize=9.5, leading=14, textColor=NAVY,
)

cover_meta_value_style = ParagraphStyle(
    "CoverMetaValue", parent=base["Normal"],
    fontName="Helvetica", fontSize=9.5, leading=14, textColor=INK,
)

toc_heading_style = ParagraphStyle(
    "TOCHeading", parent=base["Heading1"],
    fontName="Helvetica-Bold", fontSize=15, textColor=NAVY, spaceAfter=10,
)

toc_entry_style = ParagraphStyle(
    "TOCEntry", parent=base["Normal"],
    fontName="Helvetica", fontSize=10.5, leading=16, textColor=INK,
    leftIndent=4,
)

section_style = ParagraphStyle(
    "SectionTitle", parent=base["Heading1"],
    fontName="Helvetica-Bold", fontSize=16, leading=20,
    textColor=NAVY, spaceBefore=0, spaceAfter=4, keepWithNext=True,
)

section_caption_style = ParagraphStyle(
    "SectionCaption", parent=base["Normal"],
    fontName="Helvetica-Oblique", fontSize=9.7, leading=14,
    textColor=SLATE, alignment=TA_JUSTIFY, spaceBefore=4, spaceAfter=12,
)

panel_label_style = ParagraphStyle(
    "PanelLabel", parent=base["Normal"],
    fontName="Helvetica-Bold", fontSize=11.5, leading=14,
    textColor=TEAL, spaceBefore=10, spaceAfter=4, keepWithNext=True,
)

panel_caption_style = ParagraphStyle(
    "PanelCaption", parent=base["Normal"],
    fontName="Helvetica", fontSize=9.3, leading=13.5,
    textColor=SLATE, alignment=TA_JUSTIFY, leftIndent=10,
    spaceBefore=4, spaceAfter=14,
)

# =========================================================
# Footer / header (page X of Y)
# =========================================================

class NumberedCanvas(canvas_module.Canvas):
    """Canvas that renders 'Page X of Y' by buffering pages."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        total_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self._draw_furniture(total_pages)
            super().showPage()
        super().save()

    def _draw_furniture(self, total_pages):
        page = self._pageNumber
        width, _ = A4

        # top rule + running title (skip on the cover page)
        if page > 1:
            self.setStrokeColor(LINE)
            self.setLineWidth(0.6)
            self.line(0.55 * inch, A4[1] - 0.5 * inch, width - 0.55 * inch, A4[1] - 0.5 * inch)
            self.setFont("Helvetica", 8)
            self.setFillColor(SLATE)
            self.drawString(0.55 * inch, A4[1] - 0.42 * inch, REPORT_TITLE)

        # bottom rule + page number
        self.setStrokeColor(LINE)
        self.setLineWidth(0.6)
        self.line(0.55 * inch, 0.55 * inch, width - 0.55 * inch, 0.55 * inch)
        self.setFont("Helvetica", 8.5)
        self.setFillColor(SLATE)
        self.drawCentredString(width / 2, 0.35 * inch, f"Page {page} of {total_pages}")


# =========================================================
# Figure grouping
# =========================================================
# Each section = (main_key, [(item_key, image_path, letter_or_None), ...])
# A single-figure section has one item whose key == main_key and letter=None.

SECTIONS = [
    ("figure1", [
        ("figure1a", FIG / "Fig1" / "Figure1A_Length_Heatmap.png", "A"),
        ("figure1b", FIG / "Fig1" / "Figure1B_Coverage_Heatmap.png", "B"),
        ("figure1c", FIG / "Fig1" / "Figure1C_CountDistribution.png", "C"),
        ("figure1d", FIG / "Fig1" / "Figure1D_AminoAcidComposition.png", "D"),
        ("figure1e", FIG / "Fig1" / "Figure1E_Entropy_Boxplot.png", "E"),
    ]),
    ("figure2", [
        ("figure2a", FIG / "Fig2" / "Fig2A_PeptideMotifs.png", "A"),
        ("figure2b", FIG / "Fig2" / "Fig2B_Entropy.png", "B"),
        ("figure2c", FIG / "Fig2" / "Fig2C_Purity.png", "C"),
    ]),
    ("figure3", [
        ("figure3", FIG / "Fig3" / "Fig3_Purity.png", None),
    ]),
    ("figure4", [
        ("figure4", FIG / "Fig4" / "Fig4_Jaccard.png", None),
    ]),
    ("figure5", [
        ("figure5", FIG / "Fig5" / "Fig5_Combined.png", None),
    ]),
    ("figure6", [
        ("figure6a", FIG / "Fig6" / "Fig6A_GeneLength.png", "A"),
        ("figure6b", FIG / "Fig6" / "Fig6B_LCRCount.png", "B"),
        ("figure6c", FIG / "Fig6" / "Fig6C_Coverage.png", "C"),
        ("figure6d", FIG / "Fig6" / "Fig6D_EntropyRatio.png", "D"),
    ]),
    ("figure7", [
        ("figure7", FIG / "Fig7" / "Fig7_summary.png", None),
    ]),
]

MAX_IMG_WIDTH = 6.8 * inch
MAX_IMG_HEIGHT = 6.6 * inch
MAX_PANEL_HEIGHT = 3.5 * inch  # cap per-panel height when a figure has several sub-panels


def scaled_image(path, max_width, max_height):
    img = Image(str(path))
    scale = min(max_width / img.drawWidth, max_height / img.drawHeight)
    img.drawWidth *= scale
    img.drawHeight *= scale
    img.hAlign = "CENTER"
    return img


def fig_number_label(main_key):
    """'figure1' -> 'Figure 1', 'figure6' -> 'Figure 6'"""
    digits = "".join(ch for ch in main_key if ch.isdigit())
    return f"Figure {digits}"


def section_title(main_key):
    meta = CAP.get(main_key, {})
    label = fig_number_label(main_key)
    title = meta.get("title")
    return f"{label}. {title}" if title else label


# =========================================================
# Build story
# =========================================================

story = []

# ---- Cover page -----------------------------------------------------

story.append(Spacer(1, 1.4 * inch))
story.append(Paragraph(REPORT_TITLE, cover_title_style))
story.append(Spacer(1, 0.12 * inch))
story.append(Paragraph("Benchmark Report &nbsp;|&nbsp; All 7 Figures", cover_subtitle_style))
story.append(Spacer(1, 0.55 * inch))
story.append(HRFlowable(width="60%", thickness=1, color=TEAL, hAlign="CENTER"))
story.append(Spacer(1, 0.5 * inch))

today = datetime.date.today().strftime("%d %B %Y")
info_rows = [
    [Paragraph("Dataset", cover_meta_label_style), Paragraph(BASE.name, cover_meta_value_style)],
    [Paragraph("Project Folder", cover_meta_label_style), Paragraph(str(BASE), cover_meta_value_style)],
    [Paragraph("Generated", cover_meta_label_style), Paragraph(today, cover_meta_value_style)],
]
info_table = Table(info_rows, colWidths=[1.6 * inch, 4.6 * inch], hAlign="CENTER")
info_table.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, -1), PALE),
    ("BOX", (0, 0), (-1, -1), 0.6, LINE),
    ("INNERGRID", (0, 0), (-1, -1), 0.4, LINE),
    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ("LEFTPADDING", (0, 0), (-1, -1), 10),
    ("RIGHTPADDING", (0, 0), (-1, -1), 10),
    ("TOPPADDING", (0, 0), (-1, -1), 8),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
]))
story.append(info_table)

story.append(PageBreak())

# ---- Contents / list of figures --------------------------------------

story.append(Paragraph("Contents", toc_heading_style))
story.append(HRFlowable(width="100%", thickness=0.8, color=LINE, spaceAfter=10))

for main_key, items in SECTIONS:
    label = section_title(main_key)
    story.append(Paragraph(label, toc_entry_style))

story.append(PageBreak())

# ---- Figures ----------------------------------------------------------

for s_idx, (main_key, items) in enumerate(SECTIONS):
    meta = CAP.get(main_key, {})

    story.append(Paragraph(section_title(main_key), section_style))
    story.append(HRFlowable(width="100%", thickness=1.1, color=TEAL, spaceAfter=6))

    if meta.get("caption"):
        story.append(Paragraph(meta["caption"], section_caption_style))

    multi_panel = len(items) > 1
    panel_max_h = MAX_PANEL_HEIGHT if multi_panel else MAX_IMG_HEIGHT

    for item_key, img_path, letter in items:
        if not img_path.exists():
            print(f"Skipping missing image: {img_path}")
            continue

        block = []
        if letter:
            block.append(Paragraph(f"({letter})", panel_label_style))

        block.append(scaled_image(img_path, MAX_IMG_WIDTH, panel_max_h))

        panel_caption = CAP.get(item_key, {}).get("caption", "")
        if panel_caption:
            block.append(Paragraph(panel_caption, panel_caption_style))
        else:
            block.append(Spacer(1, 0.15 * inch))

        story.append(KeepTogether(block))

    if s_idx != len(SECTIONS) - 1:
        story.append(PageBreak())

# =========================================================
# Build PDF
# =========================================================

doc = SimpleDocTemplate(
    str(OUTPUT),
    pagesize=A4,
    rightMargin=0.55 * inch,
    leftMargin=0.55 * inch,
    topMargin=0.75 * inch,
    bottomMargin=0.7 * inch,
    title=REPORT_TITLE,
)

doc.build(story, canvasmaker=NumberedCanvas)

print("\n====================================")
print("Benchmark report generated")
print(OUTPUT)
print("====================================")