# rGREAT Installation on PSC Bridges-2

## Step 1: Create conda environment
conda create -n rgreat_env r-base=4.3 -y
conda activate rgreat_env

## Step 2: Install R packages
R
install.packages("BiocManager")
BiocManager::install("rGREAT")

## Step 3: Test
library(rGREAT)
