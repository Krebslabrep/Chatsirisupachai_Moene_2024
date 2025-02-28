# Bait-Capture SMF Data Analysis

This folder contains scripts to analyse and plot bait-capture SMF (mouse; Sönmezer et al., 2021) and whole-genome SMF (_Drosophila_; Krebs et al., 2017).

## Order of Scripts

1. [01_Composite_plots_Fig1B-D_FigEV1.R](/scripts/01_SMF_bait_capture/01_Composite_plots_Fig1B-D_FigEV1.R)  
   * Plots the composite SMF signal from mESCs (Sönmezer et al., 2021), other mouse cell lines (Kreibich et al., 2023), and _Drosophila_ S2 and OSC cells (Krebs et al., 2017). These plots correspond to Figures 1B-D and Figure EV1.

2. [02_Composite_plots_MM_DM_Fig1E-F.R](/scripts/01_SMF_bait_capture/02_Composite_plots_MM_DM_Fig1E-F.R)  
   * Compares promoter SMF, PRO-seq, and MNase-seq signals from mouse and _Drosophila_. These plots correspond to Figures 1E-F.

3. [03_MM_SMF_bait-capture_SM_sorting.R](/scripts/01_SMF_bait_capture/03_MM_SMF_bait-capture_SM_sorting.R)  
   * Performs promoter-state sorting on mouse SMF data (Sönmezer et al., 2021) as described in the Methods section of the manuscript. This results in a 4-bit vector classifying the state of every read among 2^4^ = 16 theoretical possibilities.

4. [04_MM_SMF_bait-capture_promoter_states.R](/scripts/01_SMF_bait_capture/04_MM_SMF_bait-capture_promoter_states.R)  
   * Assigns promoter state frequency to each mouse promoter, along with the biological interpretation (Unassigned, Nucleosome, Unbound, PIC, PIC + Pol II, and Pol II), as described in the Methods section of the manuscript and Figure EV2B.

5. [05_MM_SMF_bait-capture_corr_replicates_FigEV2C.R](/scripts/01_SMF_bait_capture/05_MM_SMF_bait-capture_corr_replicates_FigEV2C.R)  
   * Checks the correlation of promoter state frequency from mouse SMF data between replicates. This results in Figure EV2C.

6. [06_MM_corr_plot_SMF_others_FigEV2D.R](/scripts/01_SMF_bait_capture/06_MM_corr_plot_SMF_others_FigEV2D.R)  
   * Checks the correlation of promoter state frequency from mouse SMF data with other omics data. This results in Figure EV2D.

7. [07_DM_SMF_genome-wide_SM_sorting.R](/scripts/01_SMF_bait_capture/07_DM_SMF_genome-wide_SM_sorting.R)  
   * Performs promoter-state sorting on _Drosophila_ SMF data (Krebs et al., 2017) as described in the Methods section of the manuscript. This results in a 4-bit vector classifying the state of every read among 2^4^ = 16 theoretical possibilities.

8. [08_DM_SMF_genome-wide_promoter_states.R](/scripts/01_SMF_bait_capture/08_DM_SMF_genome-wide_promoter_states.R)  
   * Assigns promoter state frequency to each _Drosophila_ promoter, along with the biological interpretation (Unassigned, Nucleosome, Unbound, PIC, PIC + Pol II, and Pol II), as described in the Methods section of the manuscript and Figure EV2B.

9. [09_State_freq_vs_ChIP_seq_Fig2A_Fig2C.R](/scripts/01_SMF_bait_capture/09_State_freq_vs_ChIP_seq_Fig2A_Fig2C.R)  
   * Compares SMF-derived promoter state frequency with Pol II ChIP-seq data. This corresponds to Figures 2A and 2C.

10. [10_Comparison_PIC_PolII_occupancy_Fig2E-F.R](/scripts/01_SMF_bait_capture/10_Comparison_PIC_PolII_occupancy_Fig2E-F.R)  
    * Compares PIC and Pol II occupancy between mouse and _Drosophila_. This corresponds to Figures 2E-F.

11. [11_Single-molecule_plotting_Fig2B_Fig2D.R](/scripts/01_SMF_bait_capture/11_Single-molecule_plotting_Fig2B_Fig2D.R)  
    * Generates single-molecule plots for Figures 2B and 2D.

12. [12_CpG_density_and_promoter_accessibility_Appendix_FigS2.R](/scripts/01_SMF_bait_capture/12_CpG_density_and_promoter_accessibility_Appendix_FigS2.R)  
    * Analyses CpG density and promoter accessibility. This results in Appendix Figure S2.