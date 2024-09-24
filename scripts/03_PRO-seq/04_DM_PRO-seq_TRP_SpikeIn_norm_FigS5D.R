########## PRO-seq S2 cells: Spike-in normalisation ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 24.09.2024

library(QuasR)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(Rsamtools)
library(viridis)

### TSS
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

promReg <- promoters(TSSsc_DM, upstream = 100, downstream = 200)

### Total read counts of the BAM files
# list all files
# this file points to the BAM files of PRO-seq in S2 cells (mm10 Spike-in) following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_SpikeIn <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq_mm10SpikeIn.txt"
SpikeIn_BAMs <- read.delim(Qinput_SpikeIn)$FileName
SpikeIn_BAMs
# this file points to the BAM files of PRO-seq in S2 cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq.txt"
sample_BAMs <- read.delim(Qinput_sample)$FileName
sample_BAMs


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


### Load read counts at the promoters
promCounts <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts.rds")

### Normalisation by Spike-In
norm_mapped_reads <- function(sample){
  return(promCounts[sample] * lib_size_SpikeIn_df[grepl(sample, lib_size_SpikeIn_df$BAMfile), ]$norm_factor)
}

promCounts_norm <- lapply(colnames(promCounts), norm_mapped_reads)
promCounts_norm <- do.call(cbind, promCounts_norm)

promCounts_norm <- promCounts_norm[,c("S2_0min_R1", "S2_0min_R2",
                                      "S2_2.5min_R2", 
                                      "S2_5min_R1", "S2_5min_R2",
                                      "S2_10min_R1", "S2_10min_R2",
                                      "S2_20min_R1", "S2_20min_R2")]

saveRDS(promCounts_norm, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm.rds")


##### calculate fold-change for each timepoint and replicate comparing with 0 timepoint
FC_each_rep <- data.frame(TKO_5min_R1 = promCounts_norm$S2_5min_R1/promCounts_norm$S2_0min_R1,
                          TKO_5min_R2 = promCounts_norm$S2_5min_R2/promCounts_norm$S2_0min_R2,
                          TKO_10min_R1 = promCounts_norm$S2_10min_R1/promCounts_norm$S2_0min_R1,
                          TKO_10min_R2 = promCounts_norm$S2_10min_R2/promCounts_norm$S2_0min_R2,
                          TKO_20min_R1 = promCounts_norm$S2_20min_R1/promCounts_norm$S2_0min_R1,
                          TKO_20min_R2 = promCounts_norm$S2_20min_R2/promCounts_norm$S2_0min_R2)
row.names(FC_each_rep) <- row.names(promCounts_norm)


##### *** Figure S5D *** #####
time_point_cor <- function(timepoint, df){
  # select time point of interest
  df %>% dplyr::select(contains(timepoint)) -> tmp
  colnames(tmp) <- c("Rep1", "Rep2")
  
  # cleaning
  tmp <- tmp[complete.cases(tmp),]
  tmp <- tmp[is.finite(unlist(tmp[1])),]
  tmp <- tmp[is.finite(unlist(tmp[2])),]
  tmp <- tmp[tmp[1] > 0 & tmp[2] > 0,]
  
  # log2 transform
  tmp <- log2(tmp)
  
  # corr
  corr <- as.numeric(round(cor(tmp[1], tmp[2]), 2))
  
  # plotting
  #pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/qPRO_seq_061123/DM/QC/Fold_change_comparing_with_0min_between_reps", time_point, ".pdf"), 
  #    width = 4.5, height = 4.5, useDingbats=FALSE)
  p <- ggplot(tmp, aes(x = Rep1, y = Rep2)) +
    geom_pointdensity() +
    scale_color_viridis() +
    ggtitle(paste0("Fold change S2: ", gsub(pattern = "_", replacement = "", timepoint))) +
    xlab("Rep 1") +
    ylab("Rep 2") +
    xlim(c(-7,4)) +
    ylim(c(-7,4)) +
    annotate("text", x=-6, y=3, label= paste0("r = ", corr), size = 5) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey90"),
          axis.line = element_line(colour = "black"))
  print(p)
  dev.off()
}


lapply(c("_5min", "_10min", "_20min"), time_point_cor, df = FC_each_rep)


