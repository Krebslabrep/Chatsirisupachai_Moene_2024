# Amplicon SMF Data Analysis

This folder contains scripts to analyse and plot PRO-seq data for mouse and _Drosophila_.

## Order of Scripts

1. [01_MM_PRO-seq_TRP_RPM_norm_Appendix_FigS3A.R](/scripts/03_PRO-seq/01_MM_PRO-seq_TRP_RPM_norm_Appendix_FigS3A.R) 
   * Counts PRO-seq signal at mouse promoters, performs RPM normalisation, and plots correlation between two replicates (Appendix Figure S3A).

2. [02_DM_PRO-seq_TRP_RPM_norm_Appendix_FigS3B.R](/scripts/03_PRO-seq/02_DM_PRO-seq_TRP_RPM_norm_Appendix_FigS3B.R) 
   * Counts PRO-seq signal at _Drosophila_ promoters, performs RPM normalisation, and plots correlation between two replicates (Appendix Figure S3B).

3. [03_MM_PRO-seq_TRP_SpikeIn_norm_Appendix_FigS3C.R](/scripts/03_PRO-seq/03_MM_PRO-seq_TRP_SpikeIn_norm_Appendix_FigS3C.R)  
   * Counts PRO-seq signal at mouse promoters, performs spike-in normalisation, calculates fold change of PRO-seq signal upon TRP treatment, and plots correlation between fold changes of the two replicates (Appendix Figure S3C).

4. [04_DM_PRO-seq_TRP_SpikeIn_norm_Appendix_FigS3D.R](/scripts/03_PRO-seq/04_DM_PRO-seq_TRP_SpikeIn_norm_Appendix_FigS3D.R) 
   * Counts PRO-seq signal at _Drosophila_ promoters, performs spike-in normalisation, calculates fold change of PRO-seq signal upon TRP treatment, and plots correlation between fold changes of the two replicates (Appendix Figure S3D).

5. [05_Comparison_of_PolII_turnover_MM_DM_Fig5A.R](/scripts/03_PRO-seq/05_Comparison_of_PolII_turnover_MM_DM_Fig5A.R) 
   * Compares Pol II turnover between mouse and _Drosophila_. This corresponds to Figure 5A.

6. [06_single-site_examples_MM_DM_Fig5B-C.R](/scripts/03_PRO-seq/06_single-site_examples_MM_DM_Fig5B-C.R) 
   * Plots single-site genomic tracks for PRO-seq data in mouse and _Drosophila_. This corresponds to Figure 5B-C.

7. [07_k_means_clustering_of_PolII_turnover_Fig5D-E.R](/scripts/03_PRO-seq/07_k_means_clustering_of_PolII_turnover_Fig5D-E.R) 
   * Performs k-means clustering of Pol II turnover and plots heatmaps. This corresponds to Figure 5D-E.

8. [08_Pausing_index_calculation.R](/scripts/03_PRO-seq/08_Pausing_index_calculation.R) 
   * Calculates pausing index using PRO-seq data without TRP treatment.

9. [09_Pausing_index_distribution_Appendix_FigS1.R](/scripts/03_PRO-seq/09_Pausing_index_distribution_Appendix_FigS1.R) 
   * Compares pausing index between mouse and _Drosophila_ promoters. This corresponds to Appendix Figure S1.
