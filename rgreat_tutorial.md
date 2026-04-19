# rGREAT Installation on PSC Bridges-2

> This guide describes how to install and test **rGREAT** on PSC Bridges-2 for downstream ATAC-seq analysis.

---

## 🚀 1. Create conda environment

    conda create -n rgreat_env r-base=4.3 -y
    conda activate rgreat_env

---

## 📦 2. Install required R packages

Start R:

    R

Then run:

    install.packages("BiocManager")
    BiocManager::install("rGREAT")

---

## ✅ 3. Test installation

    library(rGREAT)

If no error appears, the installation is successful.

---

## ⚠️ Notes

- Always activate your conda environment before launching R
- Some dependencies may be required on Bridges-2

---

## 🧬 Usage in this project

- Functional enrichment analysis of ATAC-seq peaks 
- Interpretation of shared vs species-specific OCRs  
- Linking regulatory regions to biological pathways  

---

## 📚 Reference

- https://bioconductor.org/packages/release/bioc/html/rGREAT.html