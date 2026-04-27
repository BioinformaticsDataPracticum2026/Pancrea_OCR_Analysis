#!/bin/bash
#SBATCH -J trace_halper
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 8:00:00
#SBATCH -o /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/halper_%j.out
#SBATCH -e /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/halper_%j.err

# =========================================================
# TRACE HALPER mapping - zhuang21 (pancreas data test run)
# =========================================================

# Load Anaconda (provides python)
module load anaconda3

# Make the apptainer-wrapped halLiftover visible on PATH
export PATH=/ocean/projects/bio230007p/zhuang21/Pancrea_OCR_Analysis/pipeline/tools/hal_bin:$PATH

# Let Python find orthologFind module (inside halLiftover-postprocessing)
export PYTHONPATH=/ocean/projects/bio230007p/zhuang21/Pancrea_OCR_Analysis/pipeline/tools/halLiftover-postprocessing:$PYTHONPATH

# =========================================================
# Parameters
# =========================================================
HUMAN_PEAKS="/ocean/projects/bio230007p/zhuang21/output/peaks/human_pancreas.narrowPeak"
MOUSE_PEAKS="/ocean/projects/bio230007p/zhuang21/output/peaks/mouse_pancreas.narrowPeak"
OUTPUT_DIR="/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/Mapping/outputs"
HAL_FILE="/ocean/projects/bio230007p/zhuang21/data/Alignments/10plusway-master.hal"
HALPER_SCRIPT="/ocean/projects/bio230007p/zhuang21/Pancrea_OCR_Analysis/pipeline/tools/halLiftover-postprocessing/halper_map_peak_orthologs.sh"
HAL_BIN="/ocean/projects/bio230007p/zhuang21/Pancrea_OCR_Analysis/pipeline/tools/hal_bin/halLiftover"

# =========================================================
# Setup
# =========================================================
mkdir -p "$OUTPUT_DIR"

# Copy pancreas peaks, rename to "liver" so downstream TRACE scripts find them
cp "$HUMAN_PEAKS" "$OUTPUT_DIR/human_liver.narrowPeak"
cp "$MOUSE_PEAKS" "$OUTPUT_DIR/mouse_liver.narrowPeak"

# =========================================================
# Run HALPER: Human -> Mouse
# =========================================================
echo "=============================="
echo "HALPER: Human -> Mouse"
echo "Start: $(date)"
echo "=============================="

bash "$HALPER_SCRIPT" \
    -b "$OUTPUT_DIR/human_liver.narrowPeak" \
    -o "$OUTPUT_DIR" \
    -s Human \
    -t Mouse \
    -c "$HAL_FILE" \
    --halPath "$HAL_BIN"

# =========================================================
# Run HALPER: Mouse -> Human
# =========================================================
echo ""
echo "=============================="
echo "HALPER: Mouse -> Human"
echo "Start: $(date)"
echo "=============================="

bash "$HALPER_SCRIPT" \
    -b "$OUTPUT_DIR/mouse_liver.narrowPeak" \
    -o "$OUTPUT_DIR" \
    -s Mouse \
    -t Human \
    -c "$HAL_FILE" \
    --halPath "$HAL_BIN"

echo ""
echo "=============================="
echo "All HALPER jobs finished."
echo "End: $(date)"
echo "=============================="
