"""
run part or all of the pipeline:
    python3 main.py --run halper
    python3 main.py --run bed_g homer
    python3 main.py --run all
"""

import argparse
from pathlib import Path
import sys

project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root / "pipeline"))

import pipeline.halper as halper
import pipeline.bed_genome as bed_genome
import pipeline.bed_promoter_enhancer as bed_promoter_enhancer
import pipeline.homer as homer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run",
        nargs="+",
        choices=["halper", "bed", "homer", "all"],
        required=True
    )
    args = parser.parse_args()

    steps = args.run
    if "all" in steps:
        steps = ["halper", "bed", "homer"]

    if "halper" in steps:
        halper.main()
    if "bed_g" in steps:
        bed_genome.main()
    if "bed_pe" in steps:
        bed_promoter_enhancer.main()
    if "homer" in steps:
        homer.main()


if __name__ == "__main__":
    main()