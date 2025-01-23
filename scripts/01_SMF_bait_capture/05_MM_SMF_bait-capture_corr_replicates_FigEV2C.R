########## SMF bait-capture - Correlation of state frequencies between replicates ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(ggplot2)
library(GGally)
library(ggpointdensity)
library(viridis)


### Qinput
# this file points toward BAM files of bait-capture SMF in TKO mESCs from Sonmezer et al., 2021
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_SMF_bait_capture_MM.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])

### Load data
grouped_freq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix.rds")

### Combine PIC.polII and polII
combine_PICpolII <- function(state_matrix){
  tmp <- as.data.frame(state_matrix)
  tmp$polII <- tmp$polII + tmp$PIC.polII
  tmp$PIC.polII <- NULL
  return(tmp)
}

grouped_freq_list <- lapply(grouped_freq, combine_PICpolII)
names(grouped_freq_list) <- MySample


##### Plot correlation between replicates (*** Figure S2C ***)

scatter_corr <- function(state, grouped_freq_list){
  
  state_list <- lapply(grouped_freq_list, "[", state)
  names(state_list) <- names(grouped_freq_list)
  state_df <- data.frame(state_list)
  
  corr <- cor.test(state_df[,1], state_df[,2])
  coefficient <- round(corr$estimate, 2)
  
  #pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/genome_wide_sorting/scatterplot_replicates_MM_TKO_DE_bait_capture_state_freq_", state, ".pdf"),
  #    height = 4.5, width = 5.5, useDingbats = FALSE)
  p <- ggplot(state_df, aes(x = state_df[,1], y = state_df[,2])) +
    geom_pointdensity() +
    scale_color_viridis(name = "num promoter") +
    ggtitle(state) +
    xlim(0, 100) +
    ylim(0, 100) +
    xlab("State frequency: Rep 1") +
    ylab("State frequency: Rep 2") +
    annotate("text", x=10, y=90, label= paste0("r = ", coefficient), size = 5) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14))
  print(p)
  dev.off()
}

lapply(c("unassigned", "nucleosome",  "unbound",  "PIC",  "polII"), scatter_corr, grouped_freq_list)

