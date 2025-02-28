# Data preparation
This folder contains scripts to prepare all data used in the analyses.
If you are only interested in the analyses, you can skip the scripts in this folder and directly go to the `scripts/01_SMF_bait_capture` or other scripts folders.

### Order of scripts
1. [01_MM_TSSs.R](/scripts/00_data_prep/01_MM_TSSs.R)
* Annotate transcription start sites (TSSs) of mouse Refseq transcripts using CAGE-data and define TATA-containing promoters.
2. 02_MM_PRO-seq.R 
* Prepare mouse PRO-seq data from Kreibich et al., 2023 for the analysis in Figure 1E-F.
3. 03_MM_ChIP_MNase.R
* Prepare ChIP-seq data from Langer et al., 2016, this is used to define top 5% promoters by Pol II ChIP-seq signal.
* Prepare MNase-seq from Barisic et al., 2019 for the analysis in Figure 1E-F.
