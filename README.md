# ATAC-seq Project README



## Overview
This project analyzes ATAC-seq open chromatin regions (OCRs) across two species to identify conserved and species-specific regulatory elements. The workflow focuses on:
- processing mapped OCR peak files for each species
- separating OCRs into likely promoter-like and enhancer-like regions
- comparing OCR classes across species to define shared and species-specific elements
- testing for enriched sequence patterns and transcription factor motifs using HOMER

The overall goal is to assess how chromatin accessibility and regulatory architecture are conserved or diverge across species.

---

## Dependencies:
- Python version 3.6 or 3.7 (https://www.python.org/downloads/release/python-371/)
- Python libraries `matplotlib` and `numpy`
- R version >= 4.0.0 (https://www.r-project.org/)
- Conda environment

## Prerequisits and Installation:
This project is being developed on Bridges-2, a Linux HPC cluster. Therefre, it is best that the pipeline is ran on a cluster with similar level of computational power and settings.

If the software is not installed on cluster or on the device you are using to run this pipeline, install them following the instructions below:

- Halper and HalLiftover: https://github.com/ComparativeGenomicsToolkit/hal
- Bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html
- rGreat: https://github.com/jokergoo/rgreat
- Homer: http://homer.ucsd.edu/homer/


## How to Run the pipeline:
There are multiple parts in this pipeline. They can be ran all at once or separately. First navigate to the project's top directory from terminal.
```Bash
cd path_to_project_top_dir
```
To run all steps at once:
```Bash
python3 main.py --run all
```

However, if you wish to run steps separately, the first step is still mandatory for the downstream analysis:
```Bash
python3 main.py --run halper
```

For other functions, replace the keyword after `--run` with the following:

| Keyword | Function|
| -------- | -------- |
| halper | Run halLiftover+halper to map peaks from source to target species in both direction|
| bed_pe   | Run bed pipeline to separate peaks into likely enhancers and likely promoters |
| bed_g | Run bed pipeline to find OCRs that are species specific or shared |
| homer | Run homer for finding sequence pattern |




## Workflow summary

### 1\. Peak preparation using `HalLiftover` and `Halper`

#### Input: 
- `idr.optimal_peak.narrowPeak` or `idr.conservative_peak`.narrowPeak from both human and mouse. These are all the 

#### What the pipeline does:

- The pipeline first maps reproducible ATAC-seq peaks between species using `halLiftover` and `HALPER`. Starting from the `idr.optimal_peak.narrowPeak` files, peaks are lifted over in both directions: from human to mouse and from mouse to human.

- `halLiftover` uses a multiple-species genome alignment in HAL format to project each peak’s genomic coordinates from one species onto the other genome. Because raw liftover can produce fragmented or partial mappings, `HALPER` is then used to post-process these lifted regions, reconstructing orthologous regulatory elements in a more consistent and interpretable way.

- Running the mapping in both directions helps reduce directional bias and provides a more complete set of orthologous open chromatin regions for downstream comparison.

#### Output:

- Mapped OCR peak sets in the other species’ genome coordinates, including human-to-mouse and mouse-to-human orthologous peak regions, ready for overlap and conservation analysis.

### 2\. Identify shared and species-specific mapped peaks using `BEDTools`

#### Input:

- Mapped OCR peak sets `idr.optimal_peak.SourcetoTarget.HALPER.narrowPeak` generated from step 1. There should be one for Human to Mouse mapping and one from the other direction.

#### What the pipline does:

- The pipeline uses the `bed_intersect` function from `BEDTools` to compare the mapped peak sets and determine which OCRs are shared between species and which are species-specific. Peaks that overlap between the compared mapped sets are classified as shared peaks. Peaks that do not overlap a corresponding region in the other species are classified as species-specific peaks. The `-v` flag is used to find species specific peaks, and the `-u` flag is used to find overlapping but unique peeaks.

#### Output:

If ran in both direction, there should be a set of these output for each direction:
- shared_ocrs.bed
- mapped_not_open_in_target.bed

For example if we run Human (source) to Mouse (target), we will get one bed file for shared ocr, and one for Human specific peaks that are not open in mouse.


### 3\. Functional annotation of OCR sets using rGREAT
#### Input: 

For both directions from source_species to target_species:
- shared OCRs: `shared_ocrs.bed`
- Source Species Specific OCRs: `mapped_not_open_in_target_species.bed`
- background information: `mapped_source_species_halper.bed`

#### What the pipeline does:
The rGreat pipeline will analyze two sets of regions
- shared OCRs: conserved accessible regions between species
- species specific OCRs: regions mapped to target_species but not accessible in source_species.

The pipeline look at these regions and analyze their biological functions. 

#### Output:
For shared region: `shared_GO_BP.tsv` and `shared_great_object.rds`
For specific region: `specific_GO_BP.tsv` and `specific_great_object.rds`
The pipeline will also take the tsv file and plot bubble plot for the top genes. 

---

### 4\. Separate Peaks into likely-Enhancer or likely-Promoters
#### Input:
This step comes right after running halper, so we need two files:
mapped peak: `idr.optimal_peak.SourceToTarget.HALPER.narrowPeak`
original peak of Target species: `idr.optimal_peak.narrowPeak`

You also need to define a promoter distance threshold (`default = 2000bp`)

To distinguish promoter from enhancer in species specific peaks, we also need information about the locations of Transcription Start Sties for that particular species.
- for human specific: `gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed`
- for mice specific: `gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed`

#### What the pipeline does:
The pipeline first use `bedtools_closest ` to classify peaks into likely promoters and likely enhancer based on the distance threshold.

After dividing the mapped OCRs into promoter-like and enhancer-like sets in each species, it uses the `-u` flag in `bedtools_intersect` to compare the corresponding classes across species and determine which promoter-like and enhancer-like OCRs were shared. 

Finally, to identify species_specific OCRs within each class, it uses the `-v` flag in `bedtools_intersect` to subtract the shared OCRs from the full mapped OCR set for that class. 

#### Output:
For both SourceToTarget direction, the result will look different but the outputted files are the same:
- `bed_pe_source_pancrea_to_target.out.txt` contains summary of the result.
- `species_enhancer_ocrs.bed` and `species_promoter_ocrs.bed`
- `species_specific_enhancer_ocrs.bed` and `species_specific_promoter_ocrs.bed`
- `shared_enhancer_ocrs.bed` and `shared_promoter_ocrs.bed`

### 5\. Finding frequent sequence pattern:
#### Input:
- all or part of the output files from step 4, depending on your analysis goal. 
- A FASTA genome reference file: Human (hg18, hg19, hg38), Mouse (mm8, mm9, mm10)

#### What the piprline does:
Using the `findMotifsGenome.pl` function from Homer, identifies enriched de novo and known DNA motifs within genomic regions (BED files). It extracts sequences from a specified genome, accounts for background sequence composition (GC content), and produces HTML reports of motif enrichment. 

#### Output:
- One `.html` files for each of the input. 

## To Cite this repository
Ji A, Fang K, Huang Z. Pancrea_OCR_Analysis. GitHub repository, branch SF. 2026. Available at: https://github.com/BioinformaticsDataPracticum2026/Pancrea_OCR_Analysis

