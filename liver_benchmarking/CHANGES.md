# Changes Log

This file lists what I changed in the Green team's scripts to run them on our setup. Analysis logic is unchanged throughout.

---

## Common changes across all SLURM scripts

To avoid repeating the same notes for every file, these apply to all .sh scripts below:

- Log output paths updated to my own logs directory.
- Memory config standardized to mem-per-cpu=2000M.
- Removed the original authors' email notification lines.
- All hardcoded paths pointing to other students' directories were updated to point to my own directory.

---

## scripts/run_halper_mapping.sh

This is the file with the biggest changes. I used the Green team's HALPER pipeline to run our pancreas data instead of their liver data, to see how their pipeline behaves on a different tissue and compare with our own pipeline's results.

- Removed conda activate hal, switched to our apptainer-wrapped halLiftover. Updated PATH and PYTHONPATH accordingly.

- Input data: switched from liver narrowPeak files to our pancreas narrowPeak files. After copying them in, I renamed them to human_liver.narrowPeak / mouse_liver.narrowPeak so the downstream scripts (which hardcode these filenames) can find them.

- Added a Mouse to Human run (the original only does Human to Mouse).

- Passed --halPath explicitly to HALPER.

---

## scripts/run_motif_analysis.sh

- module load bedtools changed to module load bedtools/2.30.0 (pinned version).

- Genome path updated to our local mm10.fa location.

- JASPAR database changed from JASPAR2026_vertebrates_combined.meme to JASPAR2024_vertebrates.meme because we don't have the 2026 version locally. This may affect motif results.

- Removed a few echo statements at the end.

---

## scripts/run_pe_classification.sh

- Only changed ROOT to a hardcoded path pointing to my own directory.

---

## scripts/run_plots.sh

- Added set -euo pipefail (fail-fast on errors).

- R environment: switched from source activate rgreat_env (conda) to export R_LIBS_USER, because that conda env doesn't exist on our setup. I installed the required R packages into a user library instead.

---

## scripts/run_rgreat.sh

- Same R environment change as run_plots.sh.

- Added set -euo pipefail.

- Bug fix: the original cd target was a .R file instead of a directory, so it would fail at runtime. Changed it to cd into the project directory.

---

## rGREAT_Analysis/scripts/rgreat_analysis.R

- Bug fix: original code had a typo Analyis in the output directory path. Fixed to Analysis to match the actual folder name.

- Bug fix: added recursive = TRUE to dir.create() so it doesn't fail when parent directories don't exist.

- Minor formatting: 2-space indentation, spaces around equals sign. No logic changes.

---

## rGREAT_Analysis/scripts/rgreat_plots.R

- Updated indir and outdir paths only.

— Ziyi
