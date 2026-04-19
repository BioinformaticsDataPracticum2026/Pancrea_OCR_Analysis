# rGREAT Analysis for Shared, Species-Specific, and Species-Level OCRs

> This guide describes how to run **rGREAT** analysis for functional enrichment of OCR sets from the human–mouse pancreas comparison, including directional mapped OCRs (**HtM** and **MtH**) as well as full **Human OCR** and **Mouse OCR** sets.

---

## 1. Overview

We analyzed **six OCR sets** in total:

- **HtM shared OCRs**
- **HtM not_open OCRs**
- **MtH shared OCRs**
- **MtH not_open OCRs**
- **Human OCRs**
- **Mouse OCRs**

For the directional analyses:

- **shared OCRs** represent conserved accessible regions between species  
- **not_open OCRs** represent species-specific OCRs whose orthologous regions are not accessible in the other species  

In addition, we ran rGREAT on the full **Human OCR** and **Mouse OCR** peak sets to provide species-level baseline functional profiles.

---

## 2. Set run mode

Edit in either `HtM_rGREAT.R` or `MtH_rGREAT.R`:

```r
run_mode <- "shared"   # or "not_open"
```

Run each mode separately to avoid memory issues.

---

## 3. Submit job

After chaning your run mode in R scripts, run this command for directional OCR analyses:

```bash
cd /ocean/projects/bio230007p/kfang5/rgreat_analysis/scripts/
# HtM
sbatch HtM_rGREAT.slurm
# MtH
sbatch MtH_rGREAT.slurm
```

For species-level OCR analyses:

```bash
cd /ocean/projects/bio230007p/kfang5/rgreat_analysis/scripts/
# Human OCR
sbatch human_OCR_rGREAT.slurm
# Mouse OCR
sbatch mouse_OCR_rGREAT.slurm
```

---

## 4. Input files

### For HtM (human → mouse)

- shared OCRs  

```bash
/ocean/projects/bio230007p/jji5/output/bed/HtM/shared_ocrs.bed
```

- not-open OCRs

```bash
/ocean/projects/bio230007p/jji5/output/bed/HtM/mapped_not_open_in_mouse.bed
```

- background

```bash
shared(human pancreas OCR peak set): /ocean/projects/bio230007p/jji5/data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
not_open(mouse pancreas OCR peak set): /ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
```

### For MtH (mouse → human)

- shared OCRs  

```bash
/ocean/projects/bio230007p/jji5/output/bed/MtH/shared_ocrs.bed
```

- not-open OCRs

```bash
/ocean/projects/bio230007p/jji5/output/bed/MtH/mapped_not_open_in_human.bed
```

- background

```bash
shared(human pancreas OCR peak set): /ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
not_open(mouse pancreas OCR peak set): /ocean/projects/bio230007p/jji5/data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
```

- Human OCR (full human pancreas OCR peak set)

```bash
/ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz
```

- Mouse OCR (full mouse pancreas OCR peak set)

```bash
/ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz
```

---

## 5. Output

Results are saved under:

```bash
results/HtM/
├── shared/
│   ├── shared_GO_BP.tsv
│   ├── shared_great_object.rds
├── not_open/
│   ├── not_open_GO_BP.tsv
│   ├── not_open_great_object.rds

results/MtH/
├── shared/
│   ├── shared_GO_BP.tsv
│   ├── shared_great_object.rds
├── not_open/
│   ├── not_open_GO_BP.tsv
│   ├── not_open_great_object.rds

results/Human_OCR_conservative/
├── human_OCR_GO_BP.tsv
├── human_OCR_great_object.rds

results/Mouse_OCR_conservative/
├── mouse_OCR_GO_BP.tsv
├── mouse_OCR_great_object.rds
```

---

## 6. Plotting

- Plotting is handled in separate scripts:

```bash
plot_rGREAT_nonOCR.R for HtM and MtH shared / not_open OCR sets
plot_rGREAT_OCR.R for Human OCR and Mouse OCR sets
```


## 7. Notes

- Analyses are run separately for memory stability
- Input regions and background regions must be in the same coordinate system
- For directional analyses, genome build depends on the OCR set being analyzed:
- HtM shared → `mm10`
- HtM not_open → `hg38`
- MtH shared → `hg38`
- MtH not_open → `mm10`

- Human OCR and Mouse OCR analyses are run directly on species-level OCR peak sets
- Human OCR and Mouse OCR analyses do not use an explicit background file in the current scripts