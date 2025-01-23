########## Promoters' CG content analysis ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.01.2025

library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Dmelanogaster.UCSC.dm6)

### Load data
# TSSs
TSSs_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSs_MM <- resize(TSSs_MM, width = 500, fix = "center")
TSSs_MM <- keepStandardChromosomes(TSSs_MM)

TSSs_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSs_DM <- resize(TSSs_DM, width = 500, fix = "center")
TSSs_DM <- keepStandardChromosomes(TSSs_DM)


### Load SMF data
MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix_avg.rds")
MM %>% 
  rownames_to_column() %>%
  dplyr::rename(gene_id = rowname) %>%
  mutate(non_nucleosome = unbound + PIC + PIC.polII + polII,
         PolII_all = PIC.polII + polII) -> MM

DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix_avg.rds")
DM %>% 
  rownames_to_column() %>%
  dplyr::rename(gene_id = rowname) %>%
  mutate(non_nucleosome = unbound + PIC + PIC.polII + polII,
         PolII_all = PIC.polII + polII) -> DM


### calculate number of CpG within promoter regions
# Get sequences around the TSSs [-250:250] for the promoters that used in SMF analysis
MM_seq <- getSeq(Mmusculus, TSSs_MM[unlist(TSSs_MM$gene_id) %in% MM$gene_id])
DM_seq <- getSeq(Dmelanogaster, TSSs_DM[unlist(TSSs_DM$gene_id) %in% DM$gene_id])

# Count num CpG
MM_CpG <- data.frame(gene_id = names(MM_seq), CpG = vcountPattern("CG", MM_seq))
DM_CpG <- data.frame(gene_id = names(DM_seq), CpG = vcountPattern("CG", DM_seq))

# plot
plot_df <- data.frame(species = c(rep("Mouse", nrow(MM_CpG)), rep("Drosophila", nrow(DM_CpG))),
                      CpG = c(MM_CpG$CpG, DM_CpG$CpG))

p <- ggplot(plot_df, aes(x = species, y = CpG)) +
  geom_boxplot(aes(fill = species)) +
  scale_fill_manual(values = c(alpha("darkblue", alpha = 0.4), alpha("darkred", alpha = 0.4))) +
  ylim(0, 100) +
  labs(x = "species", y = "Number of CGs at promoter [-250:250]") +
  stat_compare_means(size = 5,
                     label = "p.format",
                     method = "wilcox.test") +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_DM_comparison/Revision_num_CGs_at_promoters.pdf",
#    width = 4.5, height = 4.5, useDingbats = FALSE)
print(p)
#dev.off()


### Subset promoters by CpG density
# merge CpG with SMF data
MM %>%
  left_join(MM_CpG, by = "gene_id") -> MM

DM %>%
  left_join(DM_CpG, by = "gene_id") -> DM

# calculate quartile by CpG
MM %>%
  mutate(CpG_quartile = cut(CpG,
                            breaks = quantile(CpG, probs = seq(0, 1, 0.25), na.rm = TRUE),
                            include.lowest = TRUE,
                            labels = c("Q1", "Q2", "Q3", "Q4"))) -> MM

DM %>%
  mutate(CpG_quartile = cut(CpG,
                            breaks = quantile(CpG, probs = seq(0, 1, 0.25), na.rm = TRUE),
                            include.lowest = TRUE,
                            labels = c("Q1", "Q2", "Q3", "Q4"))) -> DM

### Unbound by CpG quartile
p <- ggplot(MM, aes(x = CpG_quartile, y = unbound)) +
  geom_boxplot(aes(fill = CpG_quartile)) +
  scale_fill_manual(values = c("#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"),
                    name = "CpG quartile") +
  labs(x = "CpG abundance (low -> high)", y = "unbound (%)") +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size=15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = list(c("Q1", "Q2"), 
                                        c("Q2", "Q3"),
                                        c("Q3", "Q4"),
                                        c("Q1", "Q3"),
                                        c("Q2", "Q4"),
                                        c("Q1", "Q4")),
                     size = 4,
                     label = "p.adjust",
                     method = "wilcox.test",
                     p.adjust.method = "BH")

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/Revision_unbound_by_CGs_quartile.pdf",
#    width = 6, height = 4.5, useDingbats = FALSE)
print(p)
#dev.off()

