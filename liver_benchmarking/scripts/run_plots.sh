#!/bin/bash
#SBATCH --job-name=rgreat_plots
#SBATCH --output=/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/rgreat_plots_%j.out
#SBATCH --error=/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/rgreat_plots_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=bio230007p
#SBATCH --partition=RM-shared

set -euo pipefail

cd /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR
export R_LIBS_USER=/jet/home/zhuang21/R/library

Rscript rGREAT_Analysis/scripts/rgreat_plots.R
