# ATAC-seq Project README

## Overview
This project analyzes ATAC-seq open chromatin regions (OCRs) across multiple species to identify conserved and species-specific regulatory elements. The workflow focuses on:

- processing mapped OCR peak files for each species
- separating OCRs into likely promoter-like and enhancer-like regions
- comparing OCR classes across species to define shared and species-specific elements
- testing for enriched sequence patterns and transcription factor motifs using HOMER

The overall goal is to assess how chromatin accessibility and regulatory architecture are conserved or diverge across species.

---

## Project goals

1. Identify OCRs from ATAC-seq peak sets for each species.
2. Classify OCRs into likely promoters and likely enhancers.
3. Compare promoter-like and enhancer-like OCRs across species.
4. Define shared versus species-specific OCRs within each regulatory class.
5. Perform motif enrichment analysis to identify sequence patterns that may underlie regulatory conservation or divergence.

---

## Input data

The main inputs are mapped OCR peak files for each species.

Typical inputs include:

- peak BED files for each species
- mapped/orthologous OCR coordinates across species
- genome annotation files for assigning OCRs relative to transcription start sites (TSS)
- reference genome files required for downstream motif analysis

If applicable, peak sets were selected from the `idr.optimal` output rather than more conservative peak sets in order to balance reproducibility with sensitivity. This is often appropriate when the downstream analysis depends on capturing a comprehensive set of biologically meaningful accessible regions, especially for comparative analyses where overly conservative filtering may remove true OCRs that are weaker but reproducible.

---

## Workflow summary

### 1. Peak preparation
For each species, begin with the mapped OCR peak file. These peaks represent candidate open chromatin regions identified from ATAC-seq data and mapped into a comparable coordinate framework for cross-species analysis.

### 2. Divide OCRs into likely promoter-like and enhancer-like regions
OCRs are classified based on their distance to annotated transcription start sites.

General rationale:

- **promoter-like OCRs** are located near TSSs
- **enhancer-like OCRs** are located farther from TSSs

This distinction is useful because promoters and enhancers differ in regulatory role, sequence composition, and degree of evolutionary conservation.

### 3. Compare promoter-like and enhancer-like OCRs across species
After dividing mapped OCRs into promoter-like and enhancer-like sets in each species, overlap analysis is performed again within each class to determine:

- promoter-like OCRs shared across species
- promoter-like OCRs specific to a given species
- enhancer-like OCRs shared across species
- enhancer-like OCRs specific to a given species

Shared OCRs are identified by intersecting OCR sets across species.
Species-specific OCRs within a class are defined by subtracting shared OCRs from the mapped OCRs for that class.

This allows direct comparison of conservation patterns between promoter-associated and enhancer-associated accessible regions.

### 4. Motif enrichment / sequence pattern analysis
After identifying enhancer and promoter sets, motif enrichment analysis is performed with HOMER to determine whether the observed sequence patterns are enriched beyond chance expectation.

HOMER is run on the following peak sets:

1. enhancers for each species
2. promoters for each species
3. enhancers shared across species
4. enhancers specific to each species

These analyses help identify transcription factor binding motifs and other enriched sequence patterns that may explain conserved regulatory programs or species-specific regulatory evolution.

---

## Suggested directory structure

```text
project/
├── data/
│   ├── peaks/
│   ├── mapped_ocr/
│   ├── annotations/
│   └── genomes/
├── results/
│   ├── promoter_enhancer_split/
│   ├── shared_vs_specific/
│   ├── homer/
│   └── figures/
├── scripts/
│   ├── classify_ocr/
│   ├── intersect/
│   ├── subtract/
│   └── motif_analysis/
└── README.md
```

---

## Example analysis logic

### Promoter/enhancer classification
A typical strategy is to define promoter-like OCRs as peaks within a specified distance of a TSS, and enhancer-like OCRs as peaks outside that window. The exact threshold should be reported in the Methods section and justified with literature when possible.

### Shared versus species-specific OCRs
For each OCR class:

- use intersection across species to define shared regions
- subtract shared regions from mapped OCRs to obtain species-specific regions

### Motif analysis
For each peak set of interest:

- provide HOMER with the target BED file
- define an appropriate background set
- compare motif enrichment across shared and species-specific regulatory elements

---

## Interpretation framework

This project is built on the expectation that promoter-like OCRs may be more conserved across species than enhancer-like OCRs, because promoters are often more tightly constrained by gene regulation requirements. Enhancer-like OCRs may show more lineage- or species-specific accessibility, reflecting regulatory innovation or adaptation.

Motif enrichment results can then be used to infer:

- which regulatory programs are broadly conserved
- which transcription factor families may drive species-specific accessibility patterns
- whether shared OCRs are associated with core regulatory functions

---

## Software and tools

Common tools used in this workflow may include:

- BEDTools for genomic intersections and subtraction
- HOMER for motif enrichment analysis
- genome annotation utilities for assigning OCRs relative to TSS
- standard shell, Python, or R scripts for workflow organization and plotting

---

## Notes and assumptions

- This README is a project-level draft based on prior discussion context and may need to be adjusted to match the exact scripts, file names, thresholds, and species used in the actual project.
- The exact promoter distance threshold, species list, reference genomes, and HOMER parameters should be filled in explicitly.
- If you used `idr.optimal` peak sets, document that choice clearly in the Methods section together with the reason for preferring it over a conservative peak set.

---

## Items to fill in

Before finalizing this README, replace the placeholders below with project-specific details:

- species analyzed:
- genome builds used:
- promoter distance threshold:
- exact peak file names:
- exact overlap criteria:
- HOMER command or parameter settings:
- script names and execution order:
- output figure/table descriptions:

---

## Contact / maintainer

Add project owner, lab, and contact information here.


