########## Pol II occupancy in mouse following TRP treatment (separated by TATA, TATA-less) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.01.2025

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)

### Load data
TATA <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")
merged_freq_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_TKO_DE_amplicon_TRP_promoter_states_freq_matrix_collapsed.rds")

merged_freq_MM$DMSO_merged %>% 
  rownames_to_column() %>%
  select(rowname, polII) %>%
  dplyr::rename(gene_id = rowname,
                PolII_DMSO = polII) %>%
  left_join(TATA, by = "gene_id") %>%
  mutate(Treatment = "DMSO") -> df_DMSO

merged_freq_MM$TRP_merged %>% 
  rownames_to_column() %>%
  select(rowname, polII) %>%
  dplyr::rename(gene_id = rowname,
                PolII_TRP = polII) %>%
  left_join(TATA, by = "gene_id") %>%
  mutate(Treatment = "TRP") -> df_TRP

df_DMSO %>%
  select(gene_id, TATA, PolII_DMSO) %>%
  cbind(df_TRP %>% select(PolII_TRP)) -> df

### plot
data_long <- melt(df, id.vars = c("gene_id", "TATA"), 
                  variable.name = "Condition", value.name = "Pol II occupancy")

data_long$TATA <- ifelse(data_long$TATA, "TATA", "TATA-less")
data_long$Condition <- ifelse(data_long$Condition == "PolII_DMSO", "DMSO", "TRP")


# Create a paired boxplot grouped by TATA
p <- ggpaired(data_long, x = "Condition", y = "Pol II occupancy", id = "gene_id", color = "Condition",
         line.color = "gray", line.size = 0.5) +
  scale_color_manual(values = c("#525252", "#a50f15")) +
  facet_wrap(~ TATA) +
  xlab("Condition") +
  ylab("Pol II occupancy (%)") +
  labs(x = "Condition", y = "Pol II occupancy (%)") +
  stat_compare_means(paired = TRUE) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #legend.title = element_blank(),
        #legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/TRP_analysis/Revision_Boxplot_MM_PolII_binding_freq_TRP_treatment_TATA.pdf", 
#    width = 6, height = 4.5, useDingbats = FALSE)
print(p)
dev.off()
