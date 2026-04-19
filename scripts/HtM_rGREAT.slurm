#!/bin/bash
#SBATCH -A bio230007p
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH -J rgreat
#SBATCH -o /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/HtM_rGREAT.out
#SBATCH -e /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/HtM_rGREAT.err

source ~/.bashrc
conda activate rgreat_env
cd /ocean/projects/bio230007p/kfang5/rgreat_analysis

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

/usr/bin/time -v Rscript scripts/HtM_rGREAT.R > logs/run_rgreat.log 2>&1
