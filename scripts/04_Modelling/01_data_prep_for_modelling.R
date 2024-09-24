########## Prepare data for modelling ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 24.09.2024

library(tibble)

### Top 5% promoters
prom_list_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")
prom_list_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")

### promoter state frequencies from SMF for top 5% promoters
SMF_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix_avg.rds")
SMF_DM_top5 <- SMF_DM[row.names(SMF_DM) %in% prom_list_DM$top5,]
SMF_DM_top5 <- rownames_to_column(SMF_DM_top5, var = "gene")
write.table(SMF_DM_top5, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/Drosophila_state_frequency.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

SMF_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix_avg.rds")
SMF_MM_top5 <- SMF_MM[row.names(SMF_MM) %in% prom_list_MM$top5,]
SMF_MM_top5 <- rownames_to_column(SMF_MM_top5, var = "gene")
write.table(SMF_MM_top5, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/Mouse_state_frequency.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

### PRO-seq of top 5% promoters
DM_PROseq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm.rds")
DM_PROseq_top5 <- DM_PROseq[row.names(DM_PROseq) %in% prom_list_DM$top5,]
DM_PROseq_top5 <- rownames_to_column(DM_PROseq_top5, var = "gene")
write.table(DM_PROseq_top5, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm_top5.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

MM_PROseq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/mouse_TKO_TRP_promoter_counts_SpikeIn_norm.rds")
MM_PROseq_top5 <- MM_PROseq[row.names(MM_PROseq) %in% prom_list_MM$top5,]
MM_PROseq_top5 <- MM_PROseq_top5[!(rownames(MM_PROseq_top5) == "NM_001163521.1"),] # remove outlier
MM_PROseq_top5 <- rownames_to_column(MM_PROseq_top5, var = "gene")
write.table(MM_PROseq_top5, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/mouse_TKO_TRP_promoter_counts_SpikeIn_norm_top5.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

### promoter cluster mapping from k-means clustering of PRO-seq data
DM_cluster <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/DM_K_means_clustering_TSS_cluster_mapping.rds")
DM_cluster <- rownames_to_column(DM_cluster, var = "gene")
write.table(DM_cluster, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/Drosophila_K_means_clustering_TSS_cluster_mapping.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

MM_cluster <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/MM_K_means_clustering_TSS_cluster_mapping.rds")
MM_cluster <- rownames_to_column(MM_cluster, var = "gene")
write.table(MM_cluster, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/mouse_K_means_clustering_TSS_cluster_mapping.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

