# Amplicon SMF Data Analysis

This folder contains scripts to analyse and plot Amplicon SMF for mouse (this study) and _Drosophila_ (Krebs et al., 2017).

## Order of Scripts

1. [01_MM_SMF_Amplicon_TRP_SM_sorting.R](/scripts/02_SMF_amplicon/01_MM_SMF_Amplicon_TRP_SM_sorting.R) 
   * Performs promoter-state sorting on mouse amplicon SMF data as described in the Methods section of the manuscript. This results in a 4-bit vector classifying the state of every read among 2<sup>4</sup> = 16 theoretical possibilities.

2. [02_MM_SMF_Amplicon_TRP_promoter_states.R](/scripts/02_SMF_amplicon/02_MM_SMF_Amplicon_TRP_promoter_states.R) 
   * Assigns promoter state frequency to each mouse promoter in amplicon regions, along with the biological interpretation (Unassigned, Nucleosome, Unbound, PIC, PIC + Pol II, and Pol II), as described in the Methods section of the manuscript and Figure EV2B.

3. [03_MM_corr_plot_SMF_others_FigEV3B.R](/scripts/02_SMF_amplicon/03_MM_corr_plot_SMF_others_FigEV3B.R)  
   * Checks the correlation of promoter state frequency from mouse SMF data with other omics data. This results in Figure EV3B.

4. [04_MM_SMF_Amplicon_TRP_single-molecule_plotting_Amd1_Fig3A.R](/scripts/02_SMF_amplicon/04_MM_SMF_Amplicon_TRP_single-molecule_plotting_Amd1_Fig3A.R) 
   * Generates single-molecule plots for Figure 3A.

5. [05_MM_SMF_Amplicon_TRP_single-molecule_plotting_Rsrp1_FigEV4A.R](/scripts/02_SMF_amplicon/05_MM_SMF_Amplicon_TRP_single-molecule_plotting_Rsrp1_FigEV4A.R) 
   * Generates single-molecule plots for Figure EV4A.

6. [06_DM_SMF_Amplicon_TRP_SM_sorting.R](/scripts/02_SMF_amplicon/06_DM_SMF_Amplicon_TRP_SM_sorting.R) 
   * Performs promoter-state sorting on _Drosophila_ amplicon SMF data (Krebs et al., 2017) as described in the Methods section of the manuscript. This results in a 4-bit vector classifying the state of every read among 2<sup>4</sup> = 16 theoretical possibilities.

7. [07_DM_SMF_Amplicon_TRP_promoter_states.R](/scripts/02_SMF_amplicon/07_DM_SMF_Amplicon_TRP_promoter_states.R) 
   * Assigns promoter state frequency to each _Drosophila_ promoter in amplicon regions, along with the biological interpretation (Unassigned, Nucleosome, Unbound, PIC, PIC + Pol II, and Pol II), as described in the Methods section of the manuscript and Figure EV2B.

8. [08_MM_DM_SMF_Amplicon_TRP_promoter_states_plot_Fig3B.R](/scripts/02_SMF_amplicon/08_MM_DM_SMF_Amplicon_TRP_promoter_states_plot_Fig3B.R) 
   * Compares Pol II occupancy before and after TRP treatment. This corresponds to Figure 3B.

9. [09_MM_Amplicon_TRP_promoter_states_plot_by_TATA_FigEV4B.R](/scripts/02_SMF_amplicon/09_MM_Amplicon_TRP_promoter_states_plot_by_TATA_FigEV4B.R) 
   * Compares Pol II occupancy before and after TRP treatment, seperated by TATA-containing or TATA-less promoters. This corresponds to Figure EV4B.
