from pathlib import Path

project = Path("analysis")

folders = [
    "lcrbytools",
    "bed",
    "data/Fig1",
    "data/Fig2",
    "figures/Fig1",
    "figures/Fig2",
    "report",
    "logs",
]

for folder in folders:
    (project / folder).mkdir(parents=True, exist_ok=True)