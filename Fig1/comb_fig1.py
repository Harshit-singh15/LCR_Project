# combine_figure1.py

from PIL import Image

A = Image.open("Celegans\\Fig_outputs\\Figure1A_Length_Heatmap.png")
B = Image.open("Celegans\\Fig_outputs\\Figure1B_Coverage_Heatmap.png")
C = Image.open("Celegans\\Fig_outputs\\Figure1C_CountDistribution.png")
D = Image.open("Celegans\\Fig_outputs\\Figure1D_AminoAcidComposition.png")
E = Image.open("Celegans\\Fig_outputs\\Figure1E_Entropy_Boxplot.png")

# resize all top row images to same height
top_height = 700

A = A.resize((int(A.width * top_height / A.height), top_height))
B = B.resize((int(B.width * top_height / B.height), top_height))
C = C.resize((int(C.width * top_height / C.height), top_height))

# bottom row
bottom_height = 700

D = D.resize((1200, bottom_height))
E = E.resize((800, bottom_height))

width = A.width + B.width + C.width
height = top_height + bottom_height

canvas = Image.new("RGB", (width, height), "white")

canvas.paste(A, (0, 0))
canvas.paste(B, (A.width, 0))
canvas.paste(C, (A.width + B.width, 0))

canvas.paste(D, (0, top_height))
canvas.paste(E, (D.width, top_height))

canvas.save("Celegans\\Fig_outputs\\Figure1_Final.png")
print("Saved Figure1_Final.png")