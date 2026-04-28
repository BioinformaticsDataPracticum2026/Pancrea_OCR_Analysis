# ATAC-seq Project README

## Overview
This project analyzes ATAC-seq open chromatin regions (OCRs) across two species to identify conserved and species-specific regulatory elements. The workflow focuses on:

- mapping OCR peak files between species using halLiftover and HALPER
- identifying shared and species-specific OCRs with BEDTools
- separating OCRs into likely promoter-like and enhancer-like regions
- testing for enriched biological processes using rGREAT
- testing for enriched sequence motifs using HOMER

The overall goal is to assess how chromatin accessibility and regulatory architecture are conserved or diverge across species.

---

## Dependencies:
- Python version 3.6 or 3.7
- Python libraries: `yaml`, `matplotlib`, `numpy`
- R version >= 4.0.0
- R libraries: `rGREAT`, `rtracklayer`, `yaml`, `ggplot2`
- Conda environment
- Slurm cluster environment

## Prerequisites and Installation:
This project is being developed on Bridges-2, a Linux HPC cluster. Therefore, it is best to run the pipeline on a cluster with a similar level of computational power and software environment.

If the required software is not installed on the cluster or device you are using, install them following the instructions below （also see [Tools](#tools)）:

- HAL / halLiftover: https://github.com/ComparativeGenomicsToolkit/hal
- HALPER postprocessing: https://github.com/pfenninglab/halLiftover-postprocessing
- BEDTools: https://bedtools.readthedocs.io/en/latest/content/installation.html
- rGREAT: https://github.com/jokergoo/rGREAT
- HOMER: http://homer.ucsd.edu/homer/

---

## Table of Contents

- [Project Directory Structure and Setup Guide](#project-directory-structure-and-setup-guide)
- [Conda environment setup](#conda-environment-setup)
- [Input Data](#input-data)
- [Output Folders](#output-folders)
- [Config File Requirements](#config-file-requirements)
- [Notes on Naming Consistency](#notes-on-naming-consistency)
- [How to Run the Pipeline](#how-to-run-the-pipeline)
  - [1. Peak preparation using halLiftover and HALPER](#1-peak-preparation-using-halliftover-and-halper)
  - [2. Identify shared and species-specific mapped peaks using BEDTools](#2-identify-shared-and-species-specific-mapped-peaks-using-bedtools)
  - [3. Functional annotation of OCR sets using rGREAT](#3-functional-annotation-of-ocr-sets-using-rgreat)
  - [4. Separate peaks into likely enhancers and likely promoters](#4-separate-peaks-into-likely-enhancers-and-likely-promoters)
  - [5. Finding enriched sequence motifs using HOMER](#5-finding-enriched-sequence-motifs-using-homer)
- [Notes on Slurm Jobs](#notes-on-slurm-jobs)
- [To Cite this Repository](#to-cite-this-repository)
- [Citations for the software we used.](#citations-for-the-software-we-used)
- [Gen AI usage](#gen-ai-usage)

## Project Directory Structure and Setup Guide

### Directory Structure Tree
Below is the ideal directory structure. The config file assume that the folders are arranged this way. Therefore, it is recommended to set your directory like this (otherwise you need to change the config file):

```text
Pancrea_OCR_Analysis/
├── ATAC_env/
├── data/
├── output/
│   ├── bed/
│   ├── bed_promoter_enhancer/
│   ├── halper/
│   ├── homer/
│   └── rgreat/
├── pipeline/
├── tools/
├── config.yaml
├── environment.yml
├── main.py
└── README.md
```
Everything under `output/` can be created by the pipeline.

After installing the repo, verify that the `pipeline/` folder contains the Python and R scripts used by the pipeline:

```text
pipeline/
├── utils.py
├── halper.py
├── bed_genome.py
├── bed_promoter_enhancer.py
├── homer.py
└── r_scripts/
    ├── rGREAT.R
    └── plot_rGREAT.R
```
Download all necessary tools/software under the tools folder:

```text
tools/
├── bedtools2/
├── hal/
├── halLiftover-postprocessing/
├── homer/
└── sonLib/
```
[Back to Table of Contents](#table-of-contents)

## Conda environment setup

This project uses a Conda environment defined in `environment.yml`.


### Load Conda on the cluster

On Bridges-2, Conda may not be available by default in a new shell. First load the Anaconda module:

```bash
module spider anaconda
module load anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
```

### Option 1: Create the environment by name

```bash
conda env create -f environment.yml
conda activate ATAC_env
```

This creates an environment named `ATAC_env` in Conda's default environment location.

### Option 2 (Recommended): Create the environment inside the project folder

```bash
conda env create -p ./ATAC_env -f environment.yml
conda activate ./ATAC_env
```
(Received feedback from multiple users that Option2 is better.)

This creates the environment inside a local `conda_envs` folder in the repository, which matches the folder structure used in this project. 

### Verify the installation

After activation, check that the required Python packages are available:

```bash
python -c "import numpy, matplotlib, yaml, pandas, scipy; print('Environment setup successful')"
```
[Back to Table of Contents](#table-of-contents)

## Input Data

The user must provide the main biological input files in the `data/` directory.

Required alignment file:

```text
data/Alignments/10plusway-master.hal
```

Required ATAC-seq peak files (both optimal and conserved ok, but must be decompressed! .gz file will not work.):

```text
data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak
```

Required genome FASTA files:

```text
data/HumanGenomeInfo/hg38.fa
data/MouseGenomeInfo/mm10.fa
```

Required TSS annotation files:

```text
data/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed
data/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed
```

Optional motif database or reference data:

```text
data/CIS-BP_2.00/
```
[Back to Table of Contents](#table-of-contents)

## Output Folders

The `output/` directory stores files generated by the pipeline. Most of these folders can be created automatically by the scripts.

Recommended output structure:

```text
output/
├── halper/
├── bed/
├── bed_promoter_enhancer/
├── bed_pe/
├── homer/
└── rgreat/
```

Each pipeline step has its own temporary folder:

```text
output/halper/temp/
output/bed/temp/
output/bed_pe/temp/
output/homer/temp/
output/rgreat/temp/
```

This prevents intermediate files from different scripts from being mixed together.

[Back to Table of Contents](#table-of-contents)

## Config File Requirements

The `config.yaml` file should define `project_root` and then use relative paths from that project root.

Example:

```yaml
project_root: "/ocean/projects/bio230007p/project_folder"
conda_env: "ATAC_env"

species_1: "Human"
species_2: "Mouse"
organ: "Pancreas"

halper_repo: "tools/halLiftover-postprocessing"
hal_file: "data/Alignments/10plusway-master.hal"

species_1_peak_file: "data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak"
species_2_peak_file: "data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak"

halper_output_dir_htm: "output/halper/HtM"
halper_output_dir_mth: "output/halper/MtH"
halper_temp_dir: "output/halper/temp"

mapped_htm_file: "output/halper/HtM/idr.optimal_peak.HumanToMouse.HALPER.narrowPeak"
mapped_mth_file: "output/halper/MtH/idr.optimal_peak.MouseToHuman.HALPER.narrowPeak"

bed_output_dir_htm: "output/bed/HtM"
bed_output_dir_mth: "output/bed/MtH"
bed_g_temp_dir: "output/bed/temp"

mouse_tss_file: "data/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed"
human_tss_file: "data/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed"

bed_pe_output_dir_htm: "output/bed_promoter_enhancer/HtM"
bed_pe_output_dir_mth: "output/bed_promoter_enhancer/MtH"
bed_pe_temp_dir: "output/bed_pe/temp"

promoter_threshold: 2000

homer_output_dir: "output/homer"
homer_temp_dir: "output/homer/temp"

rgreat_output_dir: "output/rgreat"
rgreat_temp_dir: "output/rgreat/temp"

```
[Back to Table of Contents](#table-of-contents)

## Notes on Naming Consistency

Use consistent folder names across the config and scripts.

Recommended:

```text
output/rgreat/
```

Avoid mixing this with:

```text
output/rGreat/
```

The current config uses lowercase:

```yaml
rgreat_output_dir: "output/rgreat"
rgreat_temp_dir: "output/rgreat/temp"
```

Similarly, avoid using a single shared temp folder such as:

```text
output/temp/
```

The current pipeline uses script-specific temp folders:

```text
output/halper/temp/
output/bed/temp/
output/bed_pe/temp/
output/homer/temp/
output/rgreat/temp/
```
[Back to Table of Contents](#table-of-contents)

# How to Run the Pipeline:
There are multiple parts in this pipeline. They can be run all at once or separately. First navigate to the project top directory from terminal.

```bash
cd path_to_project_top_dir
```

To run all steps:

```bash
python3 main.py --run all
```

To run only HALPER:

```bash
python3 main.py --run halper
```

To run multiple selected steps:

```bash
python3 main.py --run bed_g bed_pe homer
```

For other functions, replace the keyword after `--run` with the following:

| Keyword | Function |
| -------- | -------- |
| halper | Run halLiftover and HALPER to map peaks from source to target species in both directions |
| bed_g | Run BEDTools to find shared and species-specific OCRs |
| bed_pe | Run BEDTools to separate OCRs into likely enhancers and likely promoters |
| homer | Run HOMER motif enrichment analysis |
| rgreat all | Run all rGREAT tasks |
| rgreat Human_full | Run rGREAT only for full human OCRs |
| rgreat Mouse_full | Run rGREAT only for full mouse OCRs |
| rgreat HtM_shared | Run rGREAT only for human-to-mouse shared OCRs |
| rgreat HtM_human_specific | Run rGREAT only for human-specific OCRs in the human-to-mouse direction |
| rgreat MtH_shared | Run rGREAT only for mouse-to-human shared OCRs |
| rgreat MtH_mouse_specific | Run rGREAT only for mouse-specific OCRs in the mouse-to-human direction |
| plot_rgreat all | Plot all rGREAT results |
| plot_rgreat HtM_shared | Plot one selected rGREAT result |

Examples:

```bash
python3 main.py --run rgreat all
python3 main.py --run rgreat HtM_shared
python3 main.py --run rgreat all plot_rgreat all
python3 main.py --run plot_rgreat MtH_mouse_specific
```
[Back to Table of Contents](#table-of-contents)

---


### 1. Peak preparation using `halLiftover` and `HALPER`

#### Input:
- Human `idr.optimal_peak.narrowPeak`
- Mouse `idr.optimal_peak.narrowPeak`
- HAL multiple-genome alignment file
- HALPER postprocessing repository

#### What the pipeline does:
The pipeline maps reproducible ATAC-seq peaks between species using `halLiftover` and `HALPER`.

Starting from the `idr.optimal_peak.narrowPeak` files, peaks are lifted over in both directions:

- Human to Mouse
- Mouse to Human

`halLiftover` uses a multiple-species genome alignment in HAL format to project each peak’s genomic coordinates from one species onto the other genome. Because raw liftover can produce fragmented or partial mappings, `HALPER` is then used to postprocess these lifted regions and reconstruct orthologous regulatory elements.

Running the mapping in both directions helps reduce directional bias and provides two coordinate-specific sets for downstream comparison.

#### Output:
For Human to Mouse:

```text
output/halper/HtM/idr.optimal_peak.HumanToMouse.HALPER.narrowPeak
```

For Mouse to Human:

```text
output/halper/MtH/idr.optimal_peak.MouseToHuman.HALPER.narrowPeak
```

The HALPER temp directory stores temporary BED and liftover files:

```text
output/halper/temp/
```

---

### 2. Identify shared and species-specific mapped peaks using `BEDTools`

#### Input:
- Human-to-Mouse HALPER output:
  - `output/halper/HtM/idr.optimal_peak.HumanToMouse.HALPER.narrowPeak`
- Mouse-to-Human HALPER output:
  - `output/halper/MtH/idr.optimal_peak.MouseToHuman.HALPER.narrowPeak`
- Native target species peak file:
  - Mouse peaks for the Human-to-Mouse direction
  - Human peaks for the Mouse-to-Human direction

#### What the pipeline does:
The pipeline uses `bedtools intersect` to compare mapped OCRs against native OCRs in the target species.

For each direction:

- `-u` reports mapped OCRs that overlap at least one native target OCR. These are classified as shared OCRs.
- `-v` reports mapped OCRs that do not overlap native target OCRs. These are classified as species-specific OCRs.

For Human to Mouse:

- shared OCRs are mapped human OCRs that overlap native mouse OCRs
- human-specific OCRs are mapped human OCRs that do not overlap native mouse OCRs

For Mouse to Human:

- shared OCRs are mapped mouse OCRs that overlap native human OCRs
- mouse-specific OCRs are mapped mouse OCRs that do not overlap native human OCRs

#### Output:
For Human to Mouse:

```text
output/bed/HtM/shared_ocrs.bed
output/bed/HtM/human_specific_ocrs.bed
output/bed/HtM/summary.txt
```

For Mouse to Human:

```text
output/bed/MtH/shared_ocrs.bed
output/bed/MtH/mouse_specific_ocrs.bed
output/bed/MtH/summary.txt
```

Intermediate files are stored in:

```text
output/bed/temp/
```

---

### 3. Functional annotation of OCR sets using rGREAT

#### Input:
For full species analysis:

- Human full OCRs:
  - `species_1_peak_file`
- Mouse full OCRs:
  - `species_2_peak_file`

For cross-species comparison:

- Human to Mouse:
  - `output/bed/HtM/shared_ocrs.bed`
  - `output/bed/HtM/human_specific_ocrs.bed`

- Mouse to Human:
  - `output/bed/MtH/shared_ocrs.bed`
  - `output/bed/MtH/mouse_specific_ocrs.bed`

#### What the pipeline does:
The rGREAT pipeline submits each OCR set to GREAT and retrieves enriched Gene Ontology Biological Process terms.

Genome builds are chosen based on the coordinate system of the input file:

- Human full OCRs use `hg38`
- Mouse full OCRs use `mm10`
- Human-to-Mouse mapped regions use `mm10`
- Mouse-to-Human mapped regions use `hg38`

The rGREAT tasks can be run all at once or individually through `main.py`.

#### Output:
Example outputs:

```text
output/rgreat/Human_full/human_full_GO_BP.tsv
output/rgreat/Mouse_full/mouse_full_GO_BP.tsv
output/rgreat/HtM_shared/htm_shared_GO_BP.tsv
output/rgreat/HtM_human_specific/htm_human_specific_GO_BP.tsv
output/rgreat/MtH_shared/mth_shared_GO_BP.tsv
output/rgreat/MtH_mouse_specific/mth_mouse_specific_GO_BP.tsv
```

The plotting script generates PDF bubble plots in:

```text
output/rgreat/plots/
```

---

### 4. Separate peaks into likely enhancers and likely promoters

#### Input:
This step uses:

- HALPER mapped peak file:
  - `idr.optimal_peak.SourceToTarget.HALPER.narrowPeak`
- Native target species peak file:
  - `idr.optimal_peak.narrowPeak`
- Target species TSS annotation file
- Promoter distance threshold, default:
  - `2000 bp`

For Human to Mouse, the mapped OCRs are in mouse coordinates, so the mouse TSS file is used.

For Mouse to Human, the mapped OCRs are in human coordinates, so the human TSS file is used.

#### What the pipeline does:
The pipeline uses `bedtools closest` with the `-d` option to find the nearest transcription start site for each OCR.

OCRs are classified as:

- promoter-like if the nearest TSS is within the promoter threshold
- enhancer-like if the nearest TSS is farther than the promoter threshold

After classifying mapped OCRs and native target OCRs into promoter-like and enhancer-like groups, the pipeline compares the corresponding classes across species.

For each class:

- `bedtools intersect -u` identifies shared promoter-like or enhancer-like OCRs
- `bedtools intersect -v` identifies source-species-specific promoter-like or enhancer-like OCRs

#### Output:
For Human to Mouse:

```text
output/bed_promoter_enhancer/HtM/shared_promoter_ocrs.bed
output/bed_promoter_enhancer/HtM/shared_enhancer_ocrs.bed
output/bed_promoter_enhancer/HtM/human_specific_promoter_ocrs.bed
output/bed_promoter_enhancer/HtM/human_specific_enhancer_ocrs.bed
output/bed_promoter_enhancer/HtM/summary.txt
```

For Mouse to Human:

```text
output/bed_promoter_enhancer/MtH/shared_promoter_ocrs.bed
output/bed_promoter_enhancer/MtH/shared_enhancer_ocrs.bed
output/bed_promoter_enhancer/MtH/mouse_specific_promoter_ocrs.bed
output/bed_promoter_enhancer/MtH/mouse_specific_enhancer_ocrs.bed
output/bed_promoter_enhancer/MtH/summary.txt
```

Intermediate files are stored in:

```text
output/bed_pe/temp/
```

---

### 5. Finding enriched sequence motifs using HOMER

#### Input:
HOMER uses BED files from the promoter/enhancer step and genome FASTA files.

Example input sets:

- shared enhancer OCRs
- human-specific enhancer OCRs
- mouse-specific enhancer OCRs
- promoter OCRs
- enhancer OCRs

Genome FASTA files:

- Human genome FASTA, for example `hg38.fa`
- Mouse genome FASTA, for example `mm10.fa`

#### What the pipeline does:
The pipeline uses `findMotifsGenome.pl` from HOMER to identify enriched de novo and known DNA motifs within genomic regions.

`findMotifsGenome.pl` is used because the inputs are genomic intervals in BED format.

For Human-to-Mouse outputs, the regions are in mouse coordinates, so the mouse genome is used.

For Mouse-to-Human outputs, the regions are in human coordinates, so the human genome is used.

#### Output:
Each HOMER run produces one output directory containing HTML reports and motif result files.

Example output folders:

```text
output/homer/shared_enhancer_htm/
output/homer/human_specific_enhancer_htm/
output/homer/shared_enhancer_mth/
output/homer/mouse_specific_enhancer_mth/
```

HOMER job scripts and logs are stored in:

```text
output/homer/temp/
```
[Back to Table of Contents](#table-of-contents)

---

## Notes on Slurm Jobs
Most pipeline steps generate Slurm job scripts and submit them using `sbatch`.

You do not need to manually enter an interactive node with `interact` for normal pipeline runs.

Use `interact` only for debugging commands manually on a compute node.

---

## To Cite this Repository
Ji A, Fang K, Huang Z. Pancrea_OCR_Analysis. GitHub repository, branch SF. 2026. Available at: https://github.com/BioinformaticsDataPracticum2026/Pancrea_OCR_Analysis

## Citations for the software we used.

1. Hickey G, Paten B, Earl D, et al. HAL: a hierarchical format for storing and analyzing multiple genome alignments. *Bioinformatics*. 2013;29(10):1341-1342. https://academic.oup.com/bioinformatics/article/29/10/1341/256598

2. Zhang M, Vicario DS, Rivas MV, et al. HALPER facilitates the identification of regulatory element orthologs across species. *Bioinformatics*. 2020;36(15):4339-4340.

3. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. https://academic.oup.com/bioinformatics/article/26/6/841/244688

4. Heinz S, Benner C, Spann N, et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. *Molecular Cell*. 2010;38(4):576-589. https://pubmed.ncbi.nlm.nih.gov/20513432/

5. Gu Z, Hübschmann D. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics*. 2023;39(1):btac745. https://academic.oup.com/bioinformatics/article/39/1/btac745/6832038

[Back to Table of Contents](#table-of-contents)

## Gen AI Usage:
the `main.py`, `utils.py` and textwrapper are created by Chatgpt. All script text in `make_job_script` function in each script are created manually. 

All the directory structure trees in this readme are also created by AI. The rest of the contents in readme are manually created by curated by AI for better format. 

