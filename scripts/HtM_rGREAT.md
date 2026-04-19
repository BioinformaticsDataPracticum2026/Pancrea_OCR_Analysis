# rGREAT Analysis for HtM OCRs

> This guide describes how to run **rGREAT** analysis on mapped OCRs from human to mouse (HtM) for downstream functional enrichment.

---

## 🧬 1. Overview

We analyze two sets of regions:

- shared OCRs: conserved accessible regions between species  
- not-open OCRs: regions mapped to mouse but not accessible  

Both are compared against a common background of mapped OCRs.

---

## ⚙️ 2. Set run mode

Edit in `HtM_rGREAT.R`:

```r
run_mode <- "shared"   # or "not_open"
```

Run each mode separately to avoid memory issues.

---

## 🚀 3. Submit job

```bash
cd /ocean/projects/bio230007p/kfang5/rgreat_analysis/scripts/
sbatch HtM_rGREAT.slurm
```

---

## 📂 4. Input files

- shared OCRs  

```bash
shared_ocrs.bed
```

- not-open OCRs  

```bash
mapped_not_open_in_mouse.bed
```

- background  

```bash
mapped_mouse_halper.bed
```

---

## 📊 5. Output

Results are saved under:

```bash
results/HtM/
├── shared/
│   ├── shared_GO_BP.tsv
│   ├── shared_great_object.rds
├── not_open/
│   ├── not_open_GO_BP.tsv
│   ├── not_open_great_object.rds
```

---

## ⚠️ 6. Notes

- Analyses are run separately for memory stability  
- Background is consistent across runs for fair comparison  
- Plotting is handled in a separate script  
