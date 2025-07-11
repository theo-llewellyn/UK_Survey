# UK_Survey
A Reference ITS2 DNA Barcode Library for UK Fungi associated with Bark and Ambrosia Beetles

## 1. Sequence Denoising
1. `./vsearch.sh` This script denoises reads, merges FWD and REV files, detects consensus chimaeras and clusters reads into 0% ASVs

## 2. ASV and OTU identification
1. `./dnabarcoder.sh` BLASTs ASVs against UNITE ITS2 and assigns taxonomy based on dnabarcorder taxon-specific thresholds
2. `./dynamic_clustering.R` an edited version of the dynamic clustering approach of Florence et al. (2025) to work on our data files. Requires all the prerequisites of the orginial scripts which can be found here: https://github.com/LukeLikesDirt/AusMycobiome

## 3. Plotting data
1. `./alpha_diversity.R` calculates rarefied richness per sample and visualises diversity stats across various metadata variables
2. `./taxonomy_plots.R` compares the results of dnabarcoder identifications to TBAS and RDP
3. `./plot_map.R` plots site map shown in Figure 1
