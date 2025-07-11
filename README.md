# A Reference DNA Barcode Library for UK Fungi associated with Bark and Ambrosia Beetles
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

## Data Records

The processed sequence data and metadata are available on the public NCBI SRA under BioProject accession XXX.  All other data files are available at Figshare https://doi.org/10.xxxxx. Figshare contains:
1. The OTU x Sample table
2. Full Taxonomic Identifications for each OTU as determined by dnabarcoder and BLAST against UNITE and subsequent dynamic clustering
3. An interactive KRONA plot showing taxonomy and abundance of OTUs
4. The ASV x Sample table
5. Representative ITS2 sequences for all ASVs in FASTA format
6. Representative ITS2 sequences for all OTUs in FASTA format
7. Sample metadata
8. site metadata
9. The OTU identifications of T-BAS and RDP
10. the OTU alpha diversity statistics calculated for each beetle sample

## Analysis scripts
The following scripts contain all the code to produce the ASVs and OTUs and calculate diversity statistics.

### 1. Sequence Denoising
1. `./vsearch.sh` This script denoises reads, merges FWD and REV files, detects consensus chimaeras and clusters reads into 0% ASVs

### 2. ASV and OTU identification
1. `./dnabarcoder.sh` BLASTs ASVs against UNITE ITS2 and assigns taxonomy based on dnabarcorder taxon-specific thresholds
2. `./dynamic_clustering.R` an edited version of the dynamic clustering approach of Florence et al. (2025) to work on our data files. Requires all the prerequisites of the orginial scripts which can be found here: https://github.com/LukeLikesDirt/AusMycobiome

### 3. Plotting data
1. `./alpha_diversity.R` calculates rarefied richness per sample and visualises diversity stats across various metadata variables
2. `./taxonomy_plots.R` compares the results of dnabarcoder identifications to TBAS and RDP
3. `./plot_map.R` plots site map shown in Figure 1
