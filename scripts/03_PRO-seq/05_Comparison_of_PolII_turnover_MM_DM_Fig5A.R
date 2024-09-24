########## PRO-seq: Comparison of Pol II turnover between mouse and Drosophila ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 24.09.2024

library(ggplot2)
library(reshape2)

### TSS
TSSsc_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Top 5% promoters
prom_list_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")
prom_list_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")

### Load Spike-in normalised data
mouse_norm <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/mouse_TKO_TRP_promoter_counts_SpikeIn_norm.rds")
droso_norm <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm.rds")

### Average reps
mouse_norm_avg <- data.frame(TRP_0_min = (mouse_norm$TKO_0min_R1 + mouse_norm$TKO_0min_R2)/2,
                             TRP_2.5_min = mouse_norm$TKO_2.5min_R1,
                             TRP_5_min = (mouse_norm$TKO_5min_R1 + mouse_norm$TKO_5min_R2)/2,
                             TRP_10_min = (mouse_norm$TKO_10min_R1 + mouse_norm$TKO_10min_R2)/2,
                             TRP_20_min = (mouse_norm$TKO_20min_R1 + mouse_norm$TKO_20min_R2)/2)
row.names(mouse_norm_avg) <- row.names(mouse_norm)

droso_norm_avg <- data.frame(TRP_0_min = (droso_norm$S2_0min_R1 + droso_norm$S2_0min_R2)/2,
                             TRP_2.5_min = droso_norm$S2_2.5min_R2,
                             TRP_5_min = (droso_norm$S2_5min_R1 + droso_norm$S2_5min_R2)/2,
                             TRP_10_min = (droso_norm$S2_10min_R1 + droso_norm$S2_10min_R2)/2,
                             TRP_20_min = (droso_norm$S2_20min_R1 + droso_norm$S2_20min_R2)/2)
row.names(droso_norm_avg) <- row.names(droso_norm)

### top 5%
mouse_norm_avg <- data.frame(mouse_norm_avg[rownames(mouse_norm_avg) %in% prom_list_MM$top5,])
droso_norm_avg <- data.frame(droso_norm_avg[rownames(droso_norm_avg) %in% prom_list_DM$top5,])

### Fold changes
mouse_norm_avg_FC <- mouse_norm_avg / mouse_norm_avg$TRP_0_min
droso_norm_avg_FC <- droso_norm_avg / droso_norm_avg$TRP_0_min

# remove outlier in mouse
mouse_norm_avg_FC <- mouse_norm_avg_FC[!(rownames(mouse_norm_avg_FC) == "NM_001163521.1"),]   # FC at 20 min goes to 8.48

saveRDS(mouse_norm_avg_FC, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/mouse_TKO_TRP_promoter_counts_SpikeIn_norm_avg_top5_FC.rds")
saveRDS(droso_norm_avg_FC, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm_avg_top5_FC.rds")


### function to extract time point
extract_time <- function(x){
  return(unlist(strsplit(x, split = "_"))[2])
}

# mouse
mouse_df <- melt(mouse_norm_avg_FC)
mouse_df$species <- rep("MM", nrow(mouse_df))
mouse_df$variable <- as.character(mouse_df$variable)
mouse_df$time_point <- unlist(lapply(mouse_df$variable, extract_time))

# droso
droso_df <- melt(droso_norm_avg_FC)
droso_df$species <- rep("DM", nrow(droso_df))
droso_df$variable <- as.character(droso_df$variable)
droso_df$time_point <- unlist(lapply(droso_df$variable, extract_time))

# merge
df <- rbind(mouse_df, droso_df)
df$time_point <- factor(df$time_point, levels = c(0, 2.5, 5, 10, 20))


##### *** Figure 5A *** #####
p <- ggplot(df, aes(x = time_point, y = value)) + 
  geom_boxplot(aes(fill = species), outlier.shape = NA) +
  ggtitle("Pol II turnover between species: \nTop 5% highly active promoters") +
  xlab("TRP treatment (min)") +
  ylab("Fold change of qPRO-seq signal") +
  ylim(c(0, 1.6)) +
  scale_fill_manual(values = c(alpha("darkblue",alpha=0.4), alpha("darkred",alpha=0.4)), 
                    labels = c("Drosophila", "Mouse")) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)  
dev.off()

