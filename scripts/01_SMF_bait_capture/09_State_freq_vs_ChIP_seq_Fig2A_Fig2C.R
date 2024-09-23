########## Redistribution plot of ChIP-seq signal and SMF (Figure 2A,C) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(tidyverse)
library(dplyr)
library(ggplot2)

##### Mouse #####
### Load data
# TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

# SMF state frequencies
SMF <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_promoter_states_freq_matrix_avg.rds")

# PolII ChIP-seq
ChIP_seq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_counts.rds")

ChIPdata <- ChIP_seq %>%
  filter(row.names(ChIP_seq) %in% row.names(SMF)) %>%
  dplyr::select(PolII_mESC_SRX1089844) %>%
  rename(ChIP_seq = PolII_mESC_SRX1089844) %>%
  dplyr::select(ChIP_seq) %>%
  arrange(ChIP_seq)


brk <- c(-10, seq(0, 9, 0.1), 10) #log scale to follow extreme value are collapsed


Freqmat <- SMF[rownames(ChIPdata),]
stateList <- colnames(SMF)

sp <- names(Freqmat)
mat <- Freqmat[!is.na(rowSums(Freqmat)),]
mat <- cbind(Freqmat, ChIPdata)

mat$ChIP_enrich.log2 <- log2(mat$ChIP_seq + 1)

ci <- cut(mat$ChIP_enrich.log2 , brk)

statesv <- lapply(stateList, function(j){
  us <- mat[,j]
  nus <- us/100   #quantile(us,seq(0,1,0.01),na.rm=T)[90]
  x <- tapply(nus, ci,function(x){
    if(sum(!is.na(x))>5){
      median(x, na.rm=T)
    } 
    else{
      NA
    }
  })
  return(x)
})

mt <- (do.call(rbind, statesv) * 100) + 1

smw <- length(brk)/5

mts <- apply(mt, 1, function(x){caTools::runmean(x, smw, endrule='constant')})
Q2 <- mts/rowSums(mts)
colnames(Q2) <- stateList

### Plotting
colour_set <- colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,2,3,1,5,4)]#[c(2,3,4,9)]

plot_df <- reshape2::melt(Q2)
plot_df$value <- plot_df$value * 100
colnames(plot_df) <- c("Var1", "State", "value")

##### *** Figure 2C #####
#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/composite_plots/SMF_State_freq_vs_ChIP_seq_MM.pdf",
#    width = 6, height = 4.5, useDingbats = FALSE)
p <- ggplot(plot_df, aes(x = Var1, y = value, fill = forcats::fct_rev(State))) +
  geom_bar(stat = "identity", width = 1) + 
  scale_fill_manual(values = rev(colour_set), name = "State") +
  ggtitle("Mouse") +
  ylab("State Frequency (%)") +
  xlab("Pol II ChIP-seq Rank") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.position = "right",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()



##### Drosophila #####
### Load data
# TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

# SMF state frequencies
SMF <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix_avg.rds")

# PolII ChIP-seq
ChIP_seq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_counts.rds")
ChIP_seq %>%
  mutate(PolII_avg = (PolII_S2_SRX3981728 + PolII_S2_SRX3981727)/2) -> ChIP_seq

ChIPdata <- ChIP_seq %>%
  filter(row.names(ChIP_seq) %in% row.names(SMF)) %>%
  dplyr::select(PolII_avg) %>%
  rename(ChIP_seq = PolII_avg) %>%
  dplyr::select(ChIP_seq) %>%
  arrange(ChIP_seq)


brk <- c(-10, seq(0, 16, 0.1), 17) #log scale to follow extreme value are collapsed


Freqmat <- SMF[rownames(ChIPdata),]
stateList <- colnames(SMF)

sp <- names(Freqmat)
mat <- Freqmat[!is.na(rowSums(Freqmat)),]
mat <- cbind(Freqmat, ChIPdata)

mat$ChIP_enrich.log2 <- log2(mat$ChIP_seq + 1)

ci <- cut(mat$ChIP_enrich.log2 , brk)

statesv <- lapply(stateList, function(j){
  us <- mat[,j]
  nus <- us/100   #quantile(us,seq(0,1,0.01),na.rm=T)[90]
  x <- tapply(nus, ci,function(x){
    if(sum(!is.na(x))>5){
      median(x, na.rm=T)
    } 
    else{
      NA
    }
  })
  return(x)
})

mt <- (do.call(rbind, statesv) * 100) + 1

smw <- length(brk)/5

mts <- apply(mt, 1, function(x){caTools::runmean(x, smw, endrule='constant')})
Q2 <- mts/rowSums(mts)
colnames(Q2) <- stateList

### Plotting
colour_set <- colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,2,3,1,5,4)]#[c(2,3,4,9)]

plot_df <- reshape2::melt(Q2)
plot_df$value <- plot_df$value * 100
colnames(plot_df) <- c("Var1", "State", "value")
plot_df <- plot_df[complete.cases(plot_df),]


##### *** Figure 2A #####
# pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/DM_analysis/composite_plots/SMF_State_freq_vs_ChIP_seq_DM.pdf",
#    width = 6, height = 4.5, useDingbats = FALSE)
p <- ggplot(plot_df, aes(x = Var1, y = value, fill = forcats::fct_rev(State))) +
  geom_bar(stat = "identity", width = 1) + 
  scale_fill_manual(values = rev(colour_set), name = "State") +
  ggtitle("Drosophila") +
  ylab("State Frequency (%)") +
  xlab("Pol II ChIP-seq Rank") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.position = "right",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))
print(p)
dev.off()

