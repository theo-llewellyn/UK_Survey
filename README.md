# UK_Survey
This repository contains the code associated to the paper:
<br/>

A Reference DNA Barcode Library for UK Fungi associated with Bark and Ambrosia Beetles
<br/>

**Authors**:

Angelina Ceballos-Escalera<sup>1,2,#</sup>, Theo Llewellyn<sup>1,2,#,</sup>*, John Richards<sup>1,2</sup>, Daegan Inward<sup>3</sup>, Alfried Vogler<sup>1,2</sup>
<br/>

**Affilitions**<br/>
1. Leverhulme Centre for the Holobiont, Department of Life Sciences, Imperial College London, Silwood Park Campus, Ascot, Berkshire, SL5 7PY, UK
2. Department of Life Sciences, Natural History Museum, Cromwell Road, London, SW6 7BD UK
3. Forest Research, Alice Holt Research Station, Farnham, Surrey, GU10 4LH, UK
<br/>
<sup>#</sup>These authors contributed equally

*Corresponding author: t.llewellyn19@imperial.ac.uk


## 1. Sequence Denoising
1. `./vsearch.sh` This script denoises reads, merges FWD and REV files, detects consensus chimaeras and clusters reads into 0% ASVs

## 2. ASV and OTU identification
1. `./dnabarcoder.sh` BLASTs ASVs against UNITE ITS2 and assigns taxonomy based on dnabarcorder taxon-specific thresholds
2. `./dynamic_clustering.R` an edited version of the dynamic clustering approach of Florence et al. (2025) to work on our data files. Requires all the prerequisites of the orginial scripts which can be found here: https://github.com/LukeLikesDirt/AusMycobiome

## 3. Plotting data
1. `./alpha_diversity.R` calculates rarefied richness per sample and visualises diversity stats across various metadata variables
2. `./taxonomy_plots.R` compares the results of dnabarcoder identifications to TBAS and RDP
3. `./plot_map.R` plots site map shown in Figure 1
