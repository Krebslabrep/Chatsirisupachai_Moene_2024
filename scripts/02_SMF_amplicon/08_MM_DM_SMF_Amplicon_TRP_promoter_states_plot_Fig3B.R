########## Plot state frequencies of TRP analysis (merged rep) for MM and DM ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

merged_freq_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_TKO_DE_amplicon_TRP_promoter_states_freq_matrix_collapsed.rds")
merged_freq_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/amplicon_SMF/DM_S2_DE_amplicon_TRP_promoter_states_freq_matrix_collapsed.rds")

# MM
plot_df <- data.frame(cbind(merged_freq_MM$DMSO_merged$polII, merged_freq_MM$TRP_merged$polII))
colnames(plot_df) <- c("DMSO", "TRP")
rownames(plot_df) <- rownames(merged_freq_MM$DMSO_merged)
plot_df <- melt(plot_df)
colnames(plot_df) <- c("Condition", "PolII")

plot_df$Species <- rep("Mouse", nrow(plot_df))

# DM
plot_df_DM <- data.frame(cbind(merged_freq_DM$DMSO_merged$polII, merged_freq_DM$TRP_merged$polII))
colnames(plot_df_DM) <- c("DMSO", "TRP")
rownames(plot_df_DM) <- rownames(merged_freq_DM$DMSO_merged)
plot_df_DM <- melt(plot_df_DM)
colnames(plot_df_DM) <- c("Condition", "PolII")

plot_df_DM$Species <- rep("Drosophila", nrow(plot_df_DM))

# merge df
plot_df <- rbind(plot_df, plot_df_DM)
plot_df$Species_Condition <- paste0(plot_df$Species, "-", plot_df$Condition)

### Boxplot (*** Figure 3B ***)
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/TRP_analysis/merged_rep/Boxplot_MM_DM_S2_PolII_binding_freq_TRP_treatment_merged_rep_new.pdf",
#    width = 6, height = 4.5, useDingbats = FALSE)
p <- ggboxplot(plot_df, x = "Condition", y = "PolII", color = "Species_Condition",
               add = "jitter") +
  scale_color_manual(values = c("#525252", "#08519c", "#525252", "#a50f15")) +
  labs(x = "Condition", y = "Pol II occupancy (%)") +
  stat_compare_means(paired = TRUE) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #legend.title = element_blank(),
        #legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
p <- facet(p, facet.by = "Species")
print(p)
dev.off()



