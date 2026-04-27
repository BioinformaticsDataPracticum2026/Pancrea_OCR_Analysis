# Liver Benchmarking

This folder contains code from the Green team's repo. I copied it over for two reasons: (1) to run their pipeline alongside our pancreas results as a comparison, and (2) to help test whether their code actually runs cleanly in a different environment or on someone else's setup, in case there are any bugs.

## Source

- Original repo: https://github.com/BioinformaticsDataPracticum2026/liver-ATAC-OCR
- Original authors: Green team (Bioinformatics Data Practicum 2026)
- Pulled on: April 2026

## Why

The Green team does basically the same kind of analysis we do (ATAC-seq OCR, HALPER mapping, motif analysis, rGREAT), just on a different tissue. Running their code on our side accomplishes a few things:

- Helps debug their pipeline — code that runs on the original author's machine doesn't always run elsewhere
- Sanity-checks our own methodology
- Lets us compare results across tissues
- Confirms our compute environment is set up correctly

## What's in here

- `scripts/` — 5 shell scripts: HALPER mapping, motif analysis, PE classification, plotting, rGREAT
- `r_scripts/` — 2 R scripts for rGREAT analysis and plotting

## What I changed

To get the code running under our environment (PSC Bridges-2, our directory layout), I had to change a few things. The specifics are in `CHANGES.md`. The analysis logic itself is unchanged.

Any bugs or issues I ran into during testing will be written up in the course report feedback to share with the Green team.

— Ziyi
