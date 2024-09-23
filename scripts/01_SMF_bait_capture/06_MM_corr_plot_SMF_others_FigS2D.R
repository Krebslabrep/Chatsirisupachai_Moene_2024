########## Correlation heatmap between SMF states and other datasets ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gplots)

### Load data
# TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

# SMF state frequencies
SMF <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix_avg.rds")
SMF %>% rownames_to_column() -> SMF

# RNA-seq (from Domcke et al., 2015)
RNA_seq <- readRDS('/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_geneRPKM.rds')
RNA_seq %>% as.data.frame() %>% filter(row.names(RNA_seq) %in% SMF$rowname) %>% 
  rownames_to_column() %>% gather(name, value, -rowname) %>% 
  group_by(rowname) %>% summarize(mean = mean(value)) %>% as.data.frame() %>% 
  rename(RNA_seq = mean) -> RNA_seq

# PolII ChIP-seq
ChIP_seq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_RPM.rds")
ChIP_seq %>% filter(row.names(ChIP_seq) %in% SMF$rowname) %>% rownames_to_column() %>%
  dplyr::select(c("rowname", "PolII_enrichment")) %>% rename(ChIP_seq = PolII_enrichment) -> ChIP_seq

# MNase-seq
MNase_seq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_MNase_RPM.rds")
MNase_seq %>% filter(row.names(MNase_seq) %in% SMF$rowname) %>% rownames_to_column() %>%
  dplyr::select(c("rowname", "merged_reps")) %>% rename(MNase_seq = merged_reps) -> MNase_seq

# PRO-seq
PRO_seq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_Kreibich2023_PRO-seq_RPM.rds")
PRO_seq %>% filter(row.names(PRO_seq) %in% SMF$rowname) %>% rownames_to_column() %>%
  dplyr::select(c("rowname", "TKO_avg")) %>% rename(PRO_seq = TKO_avg) -> PRO_seq

### Merge all
list_df <- list(SMF = SMF, RNA_seq = RNA_seq, ChIP_seq = ChIP_seq, MNase_seq = MNase_seq, PRO_seq = PRO_seq)
list_df %>% reduce(inner_join, by = "rowname") %>% column_to_rownames("rowname") -> merged_df
head(merged_df)

colnames(merged_df) <- c("unassigned", "nucleosome", "unbound", "PIC", "PIC + PolII", "PolII", 
                         "RNA-seq", "PolII ChIP-seq", "MNase-seq", "PRO-seq")

### Correlation matrices (*** Figure S2D ***)
spearman_corr <- cor(merged_df, method = "spearman")

brk <- seq(-1,1,0.1)
RdBu <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(length(brk) - 1))

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/corr_plot/corr_SMF_RNA_ChIP_MNase_PRO.pdf",
#    width = 6, height = 6, useDingbats = FALSE)
p <- heatmap.2(spearman_corr, 
               trace = "none", 
               dendrogram = "none",
               margins = c(10, 10), 
               breaks = brk, 
               col = RdBu,
               density.info = "none", 
               main = "Spearman correlation",
               #lwid = lwid,
               #lhei = lhei,
               key.title = "correlation")
print(p)
dev.off()

