########## Pausing Index distribution ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.01.2025

library(tidyverse)
library(ggplot2)

### TATA promoters
TATA_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")
TATA_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_TATA_promoters.rds")

### Top 5% promoters
prom_list_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")
prom_list_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")

### PI
PI_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_qPRO-seq_PI_GB_TSS300-TSS600.rds")
PI_MM %>% 
  mutate(log2PI = log2(PI),
         species = "Mouse") %>%
  drop_na() %>% 
  filter(is.finite(log2PI)) -> PI_MM_cleaned
dim(PI_MM_cleaned)    # 13459 genes

PI_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_qPRO-seq_PI_GB_TSS300-TSS600.rds")
PI_DM %>% 
  mutate(log2PI = log2(PI),
         species = "Drosophila") %>%
  drop_na() %>% 
  filter(is.finite(log2PI)) -> PI_DM_cleaned
dim(PI_DM_cleaned)    # 8466 genes

# merge dfs
df <- rbind(PI_MM_cleaned, PI_DM_cleaned)

##### Appendix Fig. S1A #####
p <- ggplot(df, aes(x = log2PI, color = species, fill = species)) +
  #geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.4, binwidth = 0.25) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c(alpha("darkblue", alpha = 0.4), alpha("darkred", alpha = 0.4))) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  labs(x = "log2(Pausing index)", y = "Density") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "right",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/qPRO_seq_061123/PI/Revision_distribution_of_PI_MM_DM.pdf",
#    width = 6, height = 4.5, useDingbats = FALSE)
print(p)
dev.off()


### Top 5%
PI_MM %>% 
  mutate(log2PI = log2(PI)) %>%
  drop_na() %>% 
  filter(is.finite(log2(PI))) %>%
  filter(genes %in% prom_list_MM$top5) %>%
  dplyr::rename(gene_id = genes) %>%
  left_join(TATA_MM, by = "gene_id") %>%
  mutate(species = "Mouse") -> PI_MM_top5

dim(PI_MM_top5)   # 1126


# Drosophila
PI_DM %>% 
  mutate(log2PI = log2(PI)) %>%
  drop_na() %>% 
  filter(is.finite(log2(PI))) %>%
  filter(genes %in% prom_list_DM$top5) %>%
  dplyr::rename(gene_id = genes) %>%
  left_join(TATA_DM, by = "gene_id") %>%
  mutate(species = "Drosophila") -> PI_DM_top5

dim(PI_DM_top5)   # 786


### merge species
df <- rbind(PI_MM_top5, PI_DM_top5)
df$species <- as.factor(df$species)
df$group <- do.call(paste, c(df[c("species", "TATA")], sep = "_"))
df$group <- as.factor(df$group)

##### Appendix Fig. S1B #####
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/qPRO_seq_061123/PI/PI_GB_TSS300-TSS600_comparison_DM_MM_Top5percent.pdf", 
#    width = 4.5, height = 4.5, useDingbats = FALSE)
p <- ggplot(df, aes(x = group, y = log2PI)) + 
  geom_boxplot(aes(fill = species)) +
  ggtitle("Pausing index") +
  ylab("log2(Pausing index)") +
  scale_x_discrete(limits = c("Drosophila_FALSE", "Drosophila_TRUE", "Mouse_FALSE", "Mouse_TRUE"),
                   labels = c("TATA-less", "TATA", "TATA-less", "TATA")) +
  scale_fill_manual(values = c(alpha("darkblue", alpha = 0.4), alpha("darkred", alpha = 0.4))) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 15, angle = 45, hjust=1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        #legend.position = "none",
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

