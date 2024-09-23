########## PIC and Pol II state abundance in mouse vs Drosophila ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

### TATA
TATA_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")
TATA_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_TATA_promoters.rds")

### Top 5% promoters
prom_list_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")
prom_list_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")

### Load data
MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix_avg.rds")
DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix_avg.rds")

### create PIC_all (PIC + PIC.polII) and polII_all (polII + PIC.polII)
MM$PIC_all <- MM$PIC + MM$PIC.polII
MM$polII_all <- MM$polII + MM$PIC.polII

DM$PIC_all <- DM$PIC + DM$PIC.polII
DM$polII_all <- DM$polII + DM$PIC.polII


### join binding frequency and TATA promoter dataframes
MM <- merge(TATA_MM, MM, by.x = "gene_id", by.y = "row.names")
MM$species <- rep("Mouse", nrow(MM))

DM <- merge(TATA_DM, DM, by.x = "gene_id", by.y = "row.names")
DM$species <- rep("Drosophila", nrow(DM))


### PIC binding frequencies (Top 5% active promoters)
MM %>% 
  filter(gene_id %in% prom_list_MM$top5) %>%
  dplyr::select(TATA, PIC_all, species) -> df1

DM %>% 
  filter(gene_id %in% prom_list_DM$top5) %>%
  dplyr::select(TATA, PIC_all, species) -> df2

df <- rbind(df1, df2)
df$group <- do.call(paste, c(df[c("species", "TATA")], sep = "_"))
df$group <- as.factor(df$group)

##### PIC occupancy (*** Figure 2E ***)
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_DM_comparison/PIC_state_abundance_DM_MM_downstream_bins.pdf",
#    height = 6, width = 4.5, useDingbats = FALSE)
p <- ggboxplot(df, x = "group", y = "PIC_all", 
               fill = "species", palette = c(alpha("darkblue", alpha = 0.4), alpha("darkred", alpha = 0.4)),
               x.order = c("Drosophila_FALSE", "Drosophila_TRUE", "Mouse_FALSE", "Mouse_TRUE"),
               xlab = "TATA",
               ylab = "state abundance (%)",
               title = "PIC") +
  scale_x_discrete(limits = c("Drosophila_FALSE", "Drosophila_TRUE", "Mouse_FALSE", "Mouse_TRUE"),
                   labels = c("TATA-less", "TATA", "TATA-less", "TATA")) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = list(c("Drosophila_FALSE", "Drosophila_TRUE"), 
                                        c("Mouse_FALSE", "Mouse_TRUE"),
                                        c("Drosophila_FALSE", "Mouse_FALSE"),
                                        c("Drosophila_TRUE", "Mouse_TRUE")), 
                     size = 4,
                     label = "p.adjust",
                     method = "wilcox.test")

print(p)
dev.off()


### PolII binding frequencies (Top 5% active promoters)
MM %>% 
  filter(gene_id %in% prom_list_MM$top5) %>%
  dplyr::select(TATA, polII_all, species) -> df1

DM %>% 
  filter(gene_id %in% prom_list_DM$top5) %>%
  dplyr::select(TATA, polII_all, species) -> df2

df <- rbind(df1, df2)
df$group <- do.call(paste, c(df[c("species", "TATA")], sep = "_"))
df$group <- as.factor(df$group)

##### Pol II occupancy (*** Figure 2F ***)
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_DM_comparison/polII_state_abundance_DM_MM_downstream_bins.pdf",
#    height = 6, width = 4.5, useDingbats = FALSE)
p <- ggboxplot(df, x = "group", y = "polII_all", 
               fill = "species", palette = c(alpha("darkblue", alpha = 0.4), alpha("darkred", alpha = 0.4)),
               x.order = c("Drosophila_FALSE", "Drosophila_TRUE", "Mouse_FALSE", "Mouse_TRUE"),
               xlab = "TATA",
               ylab = "state abundance (%)",
               title = "Pol II") +
  scale_x_discrete(limits = c("Drosophila_FALSE", "Drosophila_TRUE", "Mouse_FALSE", "Mouse_TRUE"),
                   labels = c("TATA-less", "TATA", "TATA-less", "TATA")) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = list(c("Drosophila_FALSE", "Drosophila_TRUE"), 
                                        c("Mouse_FALSE", "Mouse_TRUE"),
                                        c("Drosophila_FALSE", "Mouse_FALSE"),
                                        c("Drosophila_TRUE", "Mouse_TRUE")), 
                     size = 4,
                     label = "p.adjust",
                     method = "wilcox.test")

print(p)
dev.off()


