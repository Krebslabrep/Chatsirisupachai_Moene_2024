##### Single-site examples for Pol II turnover
# Author: Kasit Chatsirisupachai
# LastUpdate: 24.09.2024

library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(dplyr)
library(tidyverse)
library(Rsamtools)
library(QuasR)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome.Mmusculus.UCSC.mm10)

########## Mouse example: NM_001291068.1 (Polr2a) ##########
TSSsc_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Total read counts of the BAM files
# this file points to the BAM files of PRO-seq in TKO mESC cells (dm6 Spike-in) following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_SpikeIn <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_MM_TRP_PRO-seq_dm6SpikeIn.txt"
SpikeIn_BAMs <- read.delim(Qinput_SpikeIn)$FileName

# this file points to the BAM files of PRO-seq in TKO mESC cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_MM_TRP_PRO-seq.txt"
sample_BAMs <- read.delim(Qinput_sample)$FileName

# prepare countBam parameters
param <- ScanBamParam(mapqFilter = 10)
bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)

# counts mapped reads
lib_size_samples <- unlist(lapply(sample_BAMs, function(BAM){
  tmp <- countBam(BAM, param = param)
  return(tmp$records)
}))

lib_size_SpikeIn <- unlist(lapply(SpikeIn_BAMs, function(BAM){
  tmp <- countBam(BAM, param = param)
  return(tmp$records)
}))

# strsplit to get sample names
BAMs_shorten <- unlist(lapply(sample_BAMs, function(x){
  tmp <- tail(unlist(strsplit(x, split = "/")), n = 1)
  tmp <- gsub(pattern = "_no_rRNAs[.]bam", replacement = "", tmp)
  return(tmp)
}))

SpikeIn_BAMs_shorten <- unlist(lapply(SpikeIn_BAMs, function(x){
  tmp <- tail(unlist(strsplit(x, split = "/")), n = 1)
  tmp <- gsub(pattern = "_no_rRNAs[.]bam", replacement = "", tmp)
  return(tmp)
}))


lib_size_samples_df <- data.frame(BAMfile = BAMs_shorten, mapped_reads = lib_size_samples)
lib_size_SpikeIn_df <- data.frame(BAMfile = BAMs_shorten, mapped_reads = lib_size_SpikeIn)

### Calculating normalising factors
lib_size_samples_df$norm_factor <- min(lib_size_samples_df$mapped_reads) / lib_size_samples_df$mapped_reads
lib_size_SpikeIn_df$norm_factor <- min(lib_size_SpikeIn_df$mapped_reads) / lib_size_SpikeIn_df$mapped_reads


### QuasR project
genome <- "BSgenome.Mmusculus.UCSC.mm10"
paired <- "no"

proj <- qAlign(Qinput_sample, genome, paired = paired, checkOnly = TRUE)

cluObj <- makeCluster(12)

### Count reads profile at the promoters
promCounts <- qProfile(proj, query = TSSsc_MM[unlist(TSSsc_MM$gene_id == "NM_001291068.1")],
                       upstream = 500, downstream = 500,
                       orientation = "opposite", clObj = cluObj)
promCounts$coverage <- NULL
names(promCounts)

### normalise by norm_factor from SpikeIn
norm_spikeIn <- function(sample, promCounts, lib_size_SpikeIn_df){
  
  # select sample counts
  sample_counts <- promCounts[[sample]]
  sample_counts <- melt(sample_counts)
  sample_counts$Var1 <- NULL
  sample_counts$sample <- rep(sample, nrow(sample_counts))
  
  colnames(sample_counts) <- c("position", "counts", "sample")
  sample_counts <- sample_counts[,c("sample", "position", "counts")]
  
  # select norm factor
  norm_factor <- lib_size_SpikeIn_df[lib_size_SpikeIn_df$BAMfile == sample,]$norm_factor
  
  # normalise read counts
  sample_counts$norm_counts <- sample_counts$counts * norm_factor
  
  # smooth the signal
  sample_counts$smooth <- caTools::runmean(sample_counts$norm_counts, 25, endrule = "constant")
  
  return(sample_counts)
}

normCounts <- lapply(names(promCounts), norm_spikeIn, promCounts = promCounts, lib_size_SpikeIn_df = lib_size_SpikeIn_df)
names(normCounts) <- names(promCounts)

### average between replicates
normCounts_avg <- list(data.frame(sample = rep("0_min", nrow(normCounts$TKO_0min_R1)), 
                                  position = normCounts$TKO_0min_R1$position,
                                  value = (normCounts$TKO_0min_R1$smooth + normCounts$TKO_0min_R2$smooth)/2),
                       data.frame(sample = rep("2.5_min", nrow(normCounts$TKO_2.5min_R1)), 
                                  position = normCounts$TKO_2.5min_R1$position,
                                  value = normCounts$TKO_2.5min_R1$smooth),
                       data.frame(sample = rep("5_min", nrow(normCounts$TKO_5min_R1)), 
                                  position = normCounts$TKO_5min_R1$position,
                                  value = (normCounts$TKO_5min_R1$smooth + normCounts$TKO_5min_R2$smooth)/2),
                       data.frame(sample = rep("10_min", nrow(normCounts$TKO_10min_R1)), 
                                  position = normCounts$TKO_10min_R1$position,
                                  value = (normCounts$TKO_10min_R1$smooth + normCounts$TKO_10min_R2$smooth)/2),
                       data.frame(sample = rep("20_min", nrow(normCounts$TKO_20min_R1)), 
                                  position = normCounts$TKO_20min_R1$position,
                                  value = (normCounts$TKO_20min_R1$smooth + normCounts$TKO_20min_R2$smooth)/2)
)
normCounts_avg <- do.call(rbind, normCounts_avg)
normCounts_avg$sample <- factor(normCounts_avg$sample, levels = c("0_min", "2.5_min", "5_min", "10_min", "20_min"))


##### *** Figure 5C *** #####
normCounts_avg %>%
  ggplot(aes(x = position, y = value), ) +
  geom_area(fill = "darkred", alpha = 0.5) +
  ggtitle("Polr2a") +
  xlab("Position relative to TSS") +
  ylab("Normalised PRO-seq counts") +
  facet_wrap(~sample, ncol = 1, strip.position = "right") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "right",
        panel.background = element_blank(),
        #strip.background = element_blank(),
        strip.text = element_text(size=14),
        axis.line = element_line(colour = "black")) -> p

# Adding annotation to only the first facet
p + geom_rect(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
              aes(xmin = 0, xmax = 500, ymin = -3, ymax = -2),
              fill = "grey80", color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  geom_segment(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
               aes(x = 0, xend = 0, y = -2, yend = -1),
               color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  geom_segment(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
               aes(x = 0, xend = 50, y = -1, yend = -1),
               color = "black", linewidth = 0.2, inherit.aes = FALSE,
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) -> p1

print(p1)
dev.off()



########## DM example: NM_001276025.1 (Fur1) ##########
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Total read counts of the BAM files
# this file points to the BAM files of PRO-seq in S2 cells (mm10 Spike-in) following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_SpikeIn <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq_mm10SpikeIn.txt"
SpikeIn_BAMs <- read.delim(Qinput_SpikeIn)$FileName

# this file points to the BAM files of PRO-seq in S2 cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq.txt"
sample_BAMs <- read.delim(Qinput_sample)$FileName

# prepare countBam parameters
param <- ScanBamParam(mapqFilter = 10)
bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)

# counts mapped reads
lib_size_samples <- unlist(lapply(sample_BAMs, function(BAM){
  tmp <- countBam(BAM, param = param)
  return(tmp$records)
}))

lib_size_SpikeIn <- unlist(lapply(SpikeIn_BAMs, function(BAM){
  tmp <- countBam(BAM, param = param)
  return(tmp$records)
}))

# strsplit to get sample names
BAMs_shorten <- unlist(lapply(sample_BAMs, function(x){
  tmp <- tail(unlist(strsplit(x, split = "/")), n = 1)
  tmp <- gsub(pattern = "_no_rRNAs[.]bam", replacement = "", tmp)
  return(tmp)
}))

SpikeIn_BAMs_shorten <- unlist(lapply(SpikeIn_BAMs, function(x){
  tmp <- tail(unlist(strsplit(x, split = "/")), n = 1)
  tmp <- gsub(pattern = "_no_rRNAs[.]bam", replacement = "", tmp)
  return(tmp)
}))


lib_size_samples_df <- data.frame(BAMfile = BAMs_shorten, mapped_reads = lib_size_samples)
lib_size_SpikeIn_df <- data.frame(BAMfile = BAMs_shorten, mapped_reads = lib_size_SpikeIn)


### Calculating normalising factors
lib_size_samples_df$norm_factor <- min(lib_size_samples_df$mapped_reads) / lib_size_samples_df$mapped_reads
lib_size_SpikeIn_df$norm_factor <- min(lib_size_SpikeIn_df$mapped_reads) / lib_size_SpikeIn_df$mapped_reads


### QuasR project
genome <- "BSgenome.Dmelanogaster.UCSC.dm6"
paired <- "no"

proj <- qAlign(Qinput_sample, genome, paired = paired, checkOnly = TRUE)

cluObj <- makeCluster(12)


### Count reads at the promoters
promCounts <- qProfile(proj, query = TSSsc_DM[unlist(TSSsc_DM$gene_id == "NM_001276025.1")],
                       upstream = 500, downstream = 500,
                       orientation = "opposite", clObj = cluObj)
promCounts$coverage <- NULL
names(promCounts)

### normalise by norm_factor from SpikeIn
norm_spikeIn <- function(sample, promCounts, lib_size_SpikeIn_df){
  
  # select sample counts
  sample_counts <- promCounts[[sample]]
  sample_counts <- melt(sample_counts)
  sample_counts$Var1 <- NULL
  sample_counts$sample <- rep(sample, nrow(sample_counts))
  
  colnames(sample_counts) <- c("position", "counts", "sample")
  sample_counts <- sample_counts[,c("sample", "position", "counts")]
  
  # select norm factor
  norm_factor <- lib_size_SpikeIn_df[lib_size_SpikeIn_df$BAMfile == sample,]$norm_factor
  
  # normalise read counts
  sample_counts$norm_counts <- sample_counts$counts * norm_factor
  
  # smooth the signal
  sample_counts$smooth <- caTools::runmean(sample_counts$norm_counts, 25, endrule = "constant")
  
  return(sample_counts)
}

normCounts <- lapply(names(promCounts), norm_spikeIn, promCounts = promCounts, lib_size_SpikeIn_df = lib_size_SpikeIn_df)
names(normCounts) <- names(promCounts)

### average between replicates
normCounts_avg <- list(data.frame(sample = rep("0_min", nrow(normCounts$S2_0min_R1)), 
                                  position = normCounts$S2_0min_R1$position,
                                  value = (normCounts$S2_0min_R1$smooth + normCounts$S2_0min_R2$smooth)/2),
                       data.frame(sample = rep("2.5_min", nrow(normCounts$S2_2.5min_R2)), 
                                  position = normCounts$S2_2.5min_R2$position,
                                  value = normCounts$S2_2.5min_R2$smooth),
                       data.frame(sample = rep("5_min", nrow(normCounts$S2_5min_R1)), 
                                  position = normCounts$S2_5min_R1$position,
                                  value = (normCounts$S2_5min_R1$smooth + normCounts$S2_5min_R2$smooth)/2),
                       data.frame(sample = rep("10_min", nrow(normCounts$S2_10min_R1)), 
                                  position = normCounts$S2_10min_R1$position,
                                  value = (normCounts$S2_10min_R1$smooth + normCounts$S2_10min_R2$smooth)/2),
                       data.frame(sample = rep("20_min", nrow(normCounts$S2_20min_R1)), 
                                  position = normCounts$S2_20min_R1$position,
                                  value = (normCounts$S2_20min_R1$smooth + normCounts$S2_20min_R2$smooth)/2)
)
normCounts_avg <- do.call(rbind, normCounts_avg)
normCounts_avg$sample <- factor(normCounts_avg$sample, levels = c("0_min", "2.5_min", "5_min", "10_min", "20_min"))


##### *** Figure 5B *** #####
normCounts_avg %>%
  ggplot(aes(x = position, y = value), ) +
  geom_area(fill = "darkblue", alpha = 0.5) +
  ggtitle("Fur1") +
  xlab("Position relative to TSS") +
  ylab("Normalised PRO-seq counts") +
  facet_wrap(~sample, ncol = 1, strip.position = "right") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "right",
        panel.background = element_blank(),
        #strip.background = element_blank(),
        strip.text = element_text(size=14),
        axis.line = element_line(colour = "black")) -> p


# Adding annotation to only the first facet
p + geom_rect(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
              aes(xmin = 0, xmax = 500, ymin = -5, ymax = -3),
              fill = "grey80", color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  geom_segment(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
               aes(x = 0, xend = 0, y = -3, yend = -1),
               color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  geom_segment(data = subset(normCounts_avg, sample == levels(normCounts_avg$sample)[5]),
               aes(x = 0, xend = 50, y = -1, yend = -1),
               color = "black", linewidth = 0.2, inherit.aes = FALSE,
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) -> p1

print(p1)
dev.off()

