# Data Preparation

This folder contains scripts to prepare all data used in the analyses.
If you are only interested in the analyses, you can skip the scripts in this folder and directly go to the `scripts/01_SMF_bait_capture` or other script folders.

## Order of Scripts

1. [01_MM_TSSs.R](/scripts/00_data_prep/01_MM_TSSs.R)  
    * Annotate transcription start sites (TSSs) of mouse RefSeq transcripts using CAGE data and define TATA-containing promoters.

2. [02_MM_PRO-seq.R](/scripts/00_data_prep/02_MM_PRO-seq.R)  
    * Calculate mouse promoter PRO-seq signal from Kreibich et al., 2023, for the analysis in Figure 1E-F.

3. [03_MM_ChIP_MNase.R](/scripts/00_data_prep/03_MM_ChIP_MNase.R)  
    * Calculate mouse promoter ChIP-seq signal from Langer et al., 2016; this is used to define the top 5% promoters by Pol II ChIP-seq signal.
    * Calculate mouse promoter MNase-seq signal from Barisic et al., 2019, for the analysis in Figure 1E-F.

4. [04_DM_TSSs.R](/scripts/00_data_prep/04_DM_TSSs.R)  
    * Annotate transcription start sites (TSSs) of _Drosophila_ RefSeq transcripts using CAGE data and define TATA-containing promoters.

5. [05_DM_PRO-seq.R](/scripts/00_data_prep/05_DM_PRO-seq.R)  
    * Calculate _Drosophila_ promoter PRO-seq signal from Kwak et al., 2013, for the analysis in Figure 1E-F.

6. [06_DM_ChIP_MNase.R](/scripts/00_data_prep/06_DM_ChIP_MNase.R)  
    * Calculate _Drosophila_ promoter ChIP-seq signal from Tettey et al., 2019; this is used to define the top 5% promoters by Pol II ChIP-seq signal.
    * Calculate _Drosophila_ promoter MNase-seq signal from Gilchrist et al., 2010, for the analysis in Figure 1E-F.

7. [07_DM_SMF_S2.R](/scripts/00_data_prep/07_DM_SMF_S2.R)  
    * Call context methylation of Krebs et al., 2017, SMF data from _Drosophila_ S2 cells. This average methylation level is used for Figure 1E-F.

8. [08_DM_SMF_OSC.R](/scripts/00_data_prep/08_DM_SMF_OSC.R)  
    * Call context methylation of Krebs et al., 2017, SMF data from _Drosophila_ OSC cells. This average methylation level is used for Figure EV1C-D.
