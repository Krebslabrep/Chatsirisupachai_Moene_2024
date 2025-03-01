# Chatsirisupachai, Moene et al., 2025 - Mouse promoters are characterised by low occupancy and high turnover of RNA polymerase II

Kasit Chatsirisupachai<sup>#</sup>, Christina J.I. Moene<sup>#</sup>, Rozemarijn Kleinendorst, Elisa Kreibich, Nacho Molina*, and Arnaud Krebs*\
This GitHub repository contains custom scripts for data analysis and figures from the manuscript
"Mouse promoters are characterised by low occupancy and high turnover of RNA polymerase II". The manuscript has recently been accepted for publication in _Molecular Systems Biology_

BioRxiv: https://www.biorxiv.org/content/10.1101/2024.09.23.614464v1

## Abstract
The general transcription machinery and its occupancy at promoters are highly conserved across metazoans. This contrasts with the kinetics of mRNA production that considerably differ between model species such as _Drosophila_ and mouse. The molecular basis for these kinetic differences is currently unknown. Here, we used Single Molecule Footprinting to measure RNA Polymerase II (Pol II) occupancy, the fraction of DNA molecules bound, at promoters in mouse and _Drosophila_ cell lines. Single molecule data reveals that Pol II occupancy is on average 3-5 times more frequent at transcriptionally active _Drosophila_ promoters than active mouse promoters. Kinetic modelling of the occupancy states suggests that these differences in Pol II occupancy are determined by the ratio between the transcription initiation and Pol II turnover rates. We used chemical perturbation of transcription initiation to determine Pol II turnover rate in both species. Integration of these data into the model shows that infrequent Pol II occupancy in mouse is explained by the combination of high Pol II turnover and low transcription initiation rates.

## Description
data: Contains processed data used in and obtained from the analyses. Note that raw data generated in this study have been deposited in ArrayExpress under accession numbers E-MTAB-14461 and E-MTAB-14462. Publicly available datasets used in this study are listed in Table EV2 and can be found in the online version of the manuscript.

scripts: Contains all scripts used to perform the analyses in this project. The scripts are organised into subfolders based on the analysis type. Each subfolder includes a README.md file that describes the corresponding script files.

## Requirements
Most analyses were done using R (version 4.2.2). Only the modelling analysis was performed in Matlab.
### R packages
* `AnnotationDbi (1.60.2)`
* `Biostrings (2.66.0)`
* `BSgenome.Dmelanogaster.UCSC.dm6 (1.4.1)`
* `BSgenome.Mmusculus.UCSC.mm10 (1.4.3)`
* `caTools (1.18.2)`
* `cowplot (1.1.1)`
* `dplyr (1.1.4)`
* `GenomicFeatures (1.50.4)`
* `GenomicRanges (1.50.2)`
* `GGally (1.50.2)`
* `ggExtra (0.10.0)`
* `ggplot2 (3.5.1)`
* `ggpointdensity (0.1.0)`
* `ggpubr (0.6.0)`
* `gplots (3.1.3)`
* `liftOver (1.22.0)`
* `parallel (4.2.2)`
* `pheatmap (1.0.12)`
* `plyranges (1.18.0)`
* `QuasR (1.38.0)`
* `RColorBrewer (1.1.3)`
* `Rsamtools (2.14.0)`
* `readr (2.1.4)`
* `reshape2 (1.4.4)`
* `rtracklayer (1.58.0)`
* `SingleMoleculeFootprinting ('promoter' branch)`
* `stringi (1.7.12)`
* `tibble (3.2.1)`
* `tidyr (1.3.0)`
* `tidyverse (2.0.0)`
* `viridis (0.6.2)`

Note that the package `SingleMoleculeFootprinting` used in this manuscript come from the `promoter` branch, which can be installed using the following command.
```r
remotes::install_github(repo = "https://github.com/Krebslabrep/SingleMoleculeFootprinting.git", ref = "promoter", build_vignettes = FALSE)
```
