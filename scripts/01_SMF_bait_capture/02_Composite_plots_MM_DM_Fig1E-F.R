########## Composite plot for SMF, PRO-seq, and MNase-seq in TKO mES and S2 Drosophila cells ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(ggplot2)
library(cowplot)
library(reshape2)
library(GenomicRanges)

### Load data MM
TSSsc_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TATA_promoters_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")
Top_promoters_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")
AllC_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_avg_meth_all_cell_types.rds")
PRO_seq_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_Kreibich2023_PRO_seq_profile_RPM_smoothed.rds")
MNase_seq_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_MNase_profile_RPM_smoothed.rds")


### Load data DM
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TATA_promoters_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_TATA_promoters.rds")
Top_promoters_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")
ContextMeth_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_context_methylation.rds")
PRO_seq_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_PRO_seq_profile_RPM_smoothed.rds")
MNase_seq_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_MNase_profile_RPM_smoothed.rds")

av_met_DM <- as.matrix((ContextMeth_DM$SMF_DM_DE_S2_R1_MethRate + ContextMeth_DM$SMF_DM_DE_S2_R2_MethRate)/2)
colnames(av_met_DM) <- c("S2")

# add avg methylation
ContextMeth_DM$S2_avg <- av_met_DM


########## Top 5% promoters ##########
########## MM data ##########
# subset cytosines
plotC <- as.matrix(findOverlaps(resize(TSSsc_MM, fix = "center", width = 500)[Top_promoters_MM$top5], AllC_MM, ignore.strand = TRUE))

# calculate distance between Cs and corresponding TSSsc
Neg.i <- strand(TSSsc_MM[Top_promoters_MM$top5][plotC[,1]]) == "-"
orientation <- rep(1, length(TSSsc_MM[Top_promoters_MM$top5][plotC[,1]]))
orientation[as.vector(Neg.i)] <- -1

dist_from_TSS <- (start(AllC_MM[plotC[,2]]) - start(TSSsc_MM[Top_promoters_MM$top5][plotC[,1]])) * orientation

SMF_df_MM <- data.frame(distance = dist_from_TSS, SMF = (1 - AllC_MM@elementMetadata[["TKO_meth"]][plotC[,2]]) * 100)
SMF_df_MM <- SMF_df_MM[complete.cases(SMF_df_MM),]

SMF_df_avg_MM <- aggregate(.~distance, data=SMF_df_MM, mean)

# create a trend line
avg_profile_MM <- caTools::runmean(SMF_df_avg_MM$SMF, 20, endrule = "constant")
avg_profile_MM <- data.frame(distance = SMF_df_avg_MM$distance, SMF = avg_profile_MM)

### Combine PRO-seq and MNase-seq
PRO_seq_top5_MM <- PRO_seq_MM[rownames(PRO_seq_MM) %in% Top_promoters_MM$top5,]
PRO_seq_top5_MM <- apply(PRO_seq_top5_MM, 2, mean)
PRO_seq_profile_MM <- data.frame(distance = names(PRO_seq_top5_MM), PRO_seq = PRO_seq_top5_MM, row.names = NULL)

MNase_seq_top5_MM <- MNase_seq_MM[rownames(MNase_seq_MM) %in% Top_promoters_MM$top5,]
MNase_seq_top5_MM <- apply(MNase_seq_top5_MM, 2, mean)
MNase_seq_profile_MM <- data.frame(distance = names(MNase_seq_top5_MM), MNase_seq = MNase_seq_top5_MM, row.names = NULL)

avg_profile_MM <- merge(avg_profile_MM, PRO_seq_profile_MM, by = "distance")
avg_profile_MM <- merge(avg_profile_MM, MNase_seq_profile_MM, by = "distance")

avg_profile_MM$species <- rep("Mouse", nrow(avg_profile_MM))
head(avg_profile_MM)

avg_profile_MM$PRO_seq_z <- (avg_profile_MM$PRO_seq - mean(avg_profile_MM$PRO_seq))/sd(avg_profile_MM$PRO_seq)
avg_profile_MM$MNase_seq_z <- (avg_profile_MM$MNase_seq - mean(avg_profile_MM$MNase_seq))/sd(avg_profile_MM$MNase_seq)


########## DM data ##########
# subset cytosines
plotC <- as.matrix(findOverlaps(resize(TSSsc_DM, fix = "center", width = 500)[Top_promoters_DM$top5], ContextMeth_DM, ignore.strand = TRUE))

# calculate distance between Cs and corresponding TSSsc
Neg.i <- strand(TSSsc_DM[Top_promoters_DM$top5][plotC[,1]]) == "-"
orientation <- rep(1, length(TSSsc_DM[Top_promoters_DM$top5][plotC[,1]]))
orientation[as.vector(Neg.i)] <- -1

dist_from_TSS <- (start(ContextMeth_DM[plotC[,2]]) - start(TSSsc_DM[Top_promoters_DM$top5][plotC[,1]])) * orientation

SMF_df_DM <- data.frame(distance = dist_from_TSS, SMF = (1 - ContextMeth_DM@elementMetadata[["S2_avg"]][plotC[,2]]) * 100)
SMF_df_DM <- SMF_df_DM[complete.cases(SMF_df_DM),]

SMF_df_avg_DM <- aggregate(.~distance, data=SMF_df_DM, mean)

# create a trend line
avg_profile_DM <- caTools::runmean(SMF_df_avg_DM$SMF, 20, endrule = "constant")
avg_profile_DM <- data.frame(distance = SMF_df_avg_DM$distance, SMF = avg_profile_DM)

### Combine PRO-seq and MNase-seq
PRO_seq_top5_DM <- PRO_seq_DM[rownames(PRO_seq_DM) %in% Top_promoters_DM$top5,]
PRO_seq_top5_DM <- apply(PRO_seq_top5_DM, 2, mean)
PRO_seq_profile_DM <- data.frame(distance = names(PRO_seq_top5_DM), PRO_seq = PRO_seq_top5_DM, row.names = NULL)

MNase_seq_top5_DM <- MNase_seq_DM[rownames(MNase_seq_DM) %in% Top_promoters_DM$top5,]
MNase_seq_top5_DM <- apply(MNase_seq_top5_DM, 2, mean)
MNase_seq_profile_DM <- data.frame(distance = names(MNase_seq_top5_DM), MNase_seq = MNase_seq_top5_DM, row.names = NULL)

avg_profile_DM <- merge(avg_profile_DM, PRO_seq_profile_DM, by = "distance")
avg_profile_DM <- merge(avg_profile_DM, MNase_seq_profile_DM, by = "distance")

avg_profile_DM$species <- rep("Drosophila", nrow(avg_profile_DM))
head(avg_profile_DM)

avg_profile_DM$PRO_seq_z <- (avg_profile_DM$PRO_seq - mean(avg_profile_DM$PRO_seq))/sd(avg_profile_DM$PRO_seq)
avg_profile_DM$MNase_seq_z <- (avg_profile_DM$MNase_seq - mean(avg_profile_DM$MNase_seq))/sd(avg_profile_DM$MNase_seq)


########## merge df from both species
plot_df <- rbind(avg_profile_MM, avg_profile_DM)

### Plotting
# plot
#pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/composite_plots/", save_file_name),
#    height = 4.5, width = 4.5, useDingbats = FALSE)
p1 <- ggplot(plot_df, aes(x = distance, y = SMF, group = species)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  #xlab("Position relative to TSS") +
  ylab("SMF (%)") +
  ggtitle("Top 5% TSSs") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p1)
dev.off()

p2 <- ggplot(plot_df, aes(x = distance, y = MNase_seq_z)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  #xlab("Position relative to TSS") +
  ylab("MNase-seq\n(scaled)") +
  theme(plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p2)
dev.off()

p3 <- ggplot(plot_df, aes(x = distance, y = PRO_seq_z)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  xlab("Position relative to TSS") +
  ylab("PRO-seq\n(scaled)") +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p3)
dev.off()

### SMF, MNase-seq, PRO-seq *** Figure 1E ***
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_DM_comparison/SMF_MNase-seq_PRO-seq_MM_and_DM_top5_percent_TSSs.pdf",
#    width = 5, height = 8, useDingbats = FALSE)
p_arranged <- plot_grid(p1, p2, p3, 
                        ncol = 1, 
                        align = "v",
                        rel_heights = c(2, 1, 1.5))
print(p_arranged)
dev.off()



########## Top 5% TATA promoters ##########
########## MM data ##########
Top_promoters_MM_TATA <- Top_promoters_MM$top5[Top_promoters_MM$top5 %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id]

# subset cytosines
plotC <- as.matrix(findOverlaps(resize(TSSsc_MM, fix = "center", width = 500)[Top_promoters_MM_TATA], AllC_MM, ignore.strand = TRUE))

# calculate distance between Cs and corresponding TSSsc
Neg.i <- strand(TSSsc_MM[Top_promoters_MM_TATA][plotC[,1]]) == "-"
orientation <- rep(1, length(TSSsc_MM[Top_promoters_MM_TATA][plotC[,1]]))
orientation[as.vector(Neg.i)] <- -1

dist_from_TSS <- (start(AllC_MM[plotC[,2]]) - start(TSSsc_MM[Top_promoters_MM_TATA][plotC[,1]])) * orientation

SMF_df_MM <- data.frame(distance = dist_from_TSS, SMF = (1 - AllC_MM@elementMetadata[["TKO_meth"]][plotC[,2]]) * 100)
SMF_df_MM <- SMF_df_MM[complete.cases(SMF_df_MM),]

SMF_df_avg_MM <- aggregate(.~distance, data=SMF_df_MM, mean)

# create a trend line
avg_profile_MM <- caTools::runmean(SMF_df_avg_MM$SMF, 20, endrule = "constant")
avg_profile_MM <- data.frame(distance = SMF_df_avg_MM$distance, SMF = avg_profile_MM)

### Combine PRO-seq and MNase-seq
PRO_seq_top5_MM <- PRO_seq_MM[rownames(PRO_seq_MM) %in% Top_promoters_MM_TATA,]
PRO_seq_top5_MM <- apply(PRO_seq_top5_MM, 2, mean)
PRO_seq_profile_MM <- data.frame(distance = names(PRO_seq_top5_MM), PRO_seq = PRO_seq_top5_MM, row.names = NULL)

MNase_seq_top5_MM <- MNase_seq_MM[rownames(MNase_seq_MM) %in% Top_promoters_MM_TATA,]
MNase_seq_top5_MM <- apply(MNase_seq_top5_MM, 2, mean)
MNase_seq_profile_MM <- data.frame(distance = names(MNase_seq_top5_MM), MNase_seq = MNase_seq_top5_MM, row.names = NULL)

avg_profile_MM <- merge(avg_profile_MM, PRO_seq_profile_MM, by = "distance")
avg_profile_MM <- merge(avg_profile_MM, MNase_seq_profile_MM, by = "distance")

avg_profile_MM$species <- rep("Mouse", nrow(avg_profile_MM))
head(avg_profile_MM)

avg_profile_MM$PRO_seq_z <- (avg_profile_MM$PRO_seq - mean(avg_profile_MM$PRO_seq))/sd(avg_profile_MM$PRO_seq)
avg_profile_MM$MNase_seq_z <- (avg_profile_MM$MNase_seq - mean(avg_profile_MM$MNase_seq))/sd(avg_profile_MM$MNase_seq)


########## DM data ##########
Top_promoters_DM_TATA <- Top_promoters_DM$top5[Top_promoters_DM$top5 %in% TATA_promoters_DM[TATA_promoters_DM$TATA == TRUE,]$gene_id]

# subset cytosines
plotC <- as.matrix(findOverlaps(resize(TSSsc_DM, fix = "center", width = 500)[Top_promoters_DM_TATA], ContextMeth_DM, ignore.strand = TRUE))

# calculate distance between Cs and corresponding TSSsc
Neg.i <- strand(TSSsc_DM[Top_promoters_DM_TATA][plotC[,1]]) == "-"
orientation <- rep(1, length(TSSsc_DM[Top_promoters_DM_TATA][plotC[,1]]))
orientation[as.vector(Neg.i)] <- -1

dist_from_TSS <- (start(ContextMeth_DM[plotC[,2]]) - start(TSSsc_DM[Top_promoters_DM_TATA][plotC[,1]])) * orientation

SMF_df_DM <- data.frame(distance = dist_from_TSS, SMF = (1 - ContextMeth_DM@elementMetadata[["S2_avg"]][plotC[,2]]) * 100)
SMF_df_DM <- SMF_df_DM[complete.cases(SMF_df_DM),]

SMF_df_avg_DM <- aggregate(.~distance, data=SMF_df_DM, mean)

# create a trend line
avg_profile_DM <- caTools::runmean(SMF_df_avg_DM$SMF, 20, endrule = "constant")
avg_profile_DM <- data.frame(distance = SMF_df_avg_DM$distance, SMF = avg_profile_DM)

### Combine PRO-seq and MNase-seq
PRO_seq_top5_DM <- PRO_seq_DM[rownames(PRO_seq_DM) %in% Top_promoters_DM_TATA,]
PRO_seq_top5_DM <- apply(PRO_seq_top5_DM, 2, mean)
PRO_seq_profile_DM <- data.frame(distance = names(PRO_seq_top5_DM), PRO_seq = PRO_seq_top5_DM, row.names = NULL)

MNase_seq_top5_DM <- MNase_seq_DM[rownames(MNase_seq_DM) %in% Top_promoters_DM_TATA,]
MNase_seq_top5_DM <- apply(MNase_seq_top5_DM, 2, mean)
MNase_seq_profile_DM <- data.frame(distance = names(MNase_seq_top5_DM), MNase_seq = MNase_seq_top5_DM, row.names = NULL)

avg_profile_DM <- merge(avg_profile_DM, PRO_seq_profile_DM, by = "distance")
avg_profile_DM <- merge(avg_profile_DM, MNase_seq_profile_DM, by = "distance")

avg_profile_DM$species <- rep("Drosophila", nrow(avg_profile_DM))
head(avg_profile_DM)

avg_profile_DM$PRO_seq_z <- (avg_profile_DM$PRO_seq - mean(avg_profile_DM$PRO_seq))/sd(avg_profile_DM$PRO_seq)
avg_profile_DM$MNase_seq_z <- (avg_profile_DM$MNase_seq - mean(avg_profile_DM$MNase_seq))/sd(avg_profile_DM$MNase_seq)


########## merge df from both species
plot_df <- rbind(avg_profile_MM, avg_profile_DM)

### Plotting
# plot
#pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/composite_plots/", save_file_name),
#    height = 4.5, width = 4.5, useDingbats = FALSE)
p1 <- ggplot(plot_df, aes(x = distance, y = SMF, group = species)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  #xlab("Position relative to TSS") +
  ylab("SMF (%)") +
  ggtitle("Top 5% TSSs - TATA-containing") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p1)
dev.off()

p2 <- ggplot(plot_df, aes(x = distance, y = MNase_seq_z)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  #xlab("Position relative to TSS") +
  ylab("MNase-seq\n(scaled)") +
  theme(plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p2)
dev.off()

p3 <- ggplot(plot_df, aes(x = distance, y = PRO_seq_z)) +
  geom_line(linewidth = 1.5, aes(color = species)) +
  scale_color_manual(values=c("darkblue", "darkred")) +
  xlab("Position relative to TSS") +
  ylab("PRO-seq\n(scaled)") +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p3)
dev.off()

### SMF, MNase-seq, PRO-seq  *** Figure 1F ***
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_DM_comparison/SMF_MNase-seq_PRO-seq_MM_and_DM_top5_percent_TSSs_TATA.pdf",
#    width = 5, height = 8, useDingbats = FALSE)
p_arranged <- plot_grid(p1, p2, p3, 
                        ncol = 1, 
                        align = "v",
                        rel_heights = c(2, 1, 1.5))
print(p_arranged)
dev.off()
