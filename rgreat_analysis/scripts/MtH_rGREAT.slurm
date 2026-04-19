#!/bin/bash
#SBATCH -A bio230007p
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH -t 02:00:00
#SBATCH -J rgreat

# ===== not_open =====
#SBATCH -o /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/MtH_not_open.out
#SBATCH -e /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/MtH_not_open.err

# ===== shared =====
##SBATCH -o /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/MtH_shared.out
##SBATCH -e /ocean/projects/bio230007p/kfang5/rgreat_analysis/logs/MtH_shared.err

source ~/.bashrc
conda activate rgreat_env
cd /ocean/projects/bio230007p/kfang5/rgreat_analysis

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
#/usr/bin/time -v Rscript scripts/MtH_rGREAT.R > logs/run_MtH_shared.log 2>&1
/usr/bin/time -v Rscript scripts/MtH_rGREAT.R > logs/run_MtH_not_open.log 2>&1