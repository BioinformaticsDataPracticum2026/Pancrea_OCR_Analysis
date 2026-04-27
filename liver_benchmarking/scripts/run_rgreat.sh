#!/bin/bash
#SBATCH --job-name=rgreat
#SBATCH --output=/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/rgreat_%j.out
#SBATCH --error=/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/rgreat_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=bio230007p
#SBATCH --partition=RM-shared

set -euo pipefail

cd /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR

# 使用系统 R，加载用户装的 R 包
export R_LIBS_USER=/jet/home/zhuang21/R/library

Rscript rGREAT_Analysis/scripts/rgreat_analysis.R
