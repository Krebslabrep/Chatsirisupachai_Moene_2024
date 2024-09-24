########## PRO-seq S2 cells: RPM normalisation ##########
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

### QuasR project
# Qinput
# this file points to the BAM files of PRO-seq in TKO mESC cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq.txt"
genome <- "BSgenome.Dmelanogaster.UCSC.dm6"
paired <- "no"

proj <- qAlign(Qinput_sample, genome, paired = paired, checkOnly = TRUE)

cluObj <- makeCluster(12)

### Count reads at the promoters
promCounts <- data.frame(qCount(proj, query = promReg, orientation = "opposite", clObj = cluObj))
promCounts$width <- NULL

saveRDS(promCounts, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts.rds")


### Normalisation by RPM
# function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

align_stats <- as.data.frame(alignmentStats(proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

res <- lapply(colnames(promCounts), function(x){
  return(unlist(lapply(promCounts[[x]], RPM, all_mapped = align_stats[rownames(align_stats) == x,]$mapped)))
})

res_df <- as.data.frame(do.call(cbind, res))
colnames(res_df) <- colnames(promCounts)
rownames(res_df) <- rownames(promCounts)

promCounts_RPM <- res_df[,c("S2_0min_R1", "S2_0min_R2", "S2_2.5min_R2", "S2_5min_R1", "S2_5min_R2",
                            "S2_10min_R1", "S2_10min_R2", "S2_20min_R1", "S2_20min_R2")]

saveRDS(promCounts_RPM, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_RPM.rds")



##### *** Figure S5B *** #####
genes <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected_whole_gene.rds")
# keep only genes that are longer than 600 bases 
# as initially we check also the correlation between PRO-seq signal at the gene body, which also highly correlated. But we didn't put this gene body correlation in the manuscript.
genes <- genes[width(genes) >= 600, ]

promCounts_RPM <- promCounts_RPM[row.names(promCounts_RPM) %in% genes$gene_id,]

### Correlation between replicates
plot_cor <- function(timepoint, df){
  
  # select timepoint of interest
  df %>% dplyr::select(contains(timepoint)) -> tmp
  
  rep_corr <- as.numeric(cor(tmp[1], tmp[2]))
  rep_corr <- round(rep_corr, 2)
  
  # log 10
  tmp <- log10(tmp)
  
  #pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/qPRO_seq_061123/MM/QC/scatter_between_reps", timepoint, "Promoter.pdf"),
  #    width = 4.5, height = 4.5, useDingbats = FALSE)
  p <- ggplot(tmp, aes(x = .data[[colnames(tmp)[1]]], y = .data[[colnames(tmp)[2]]])) +
    geom_pointdensity() +
    scale_color_viridis() +
    ggtitle(paste0("TKO ", (gsub("_", "", timepoint)))) +
    xlab("log10(RPM) Rep 1") +
    ylab("log10(RPM) Rep 2") +
    annotate("text", x=-1, y=3, label= paste0("r = ", rep_corr), size = 5) +
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

lapply(c("_0min_", "_5min_", "_10min_", "_20min_"), plot_cor, df = promCounts_RPM)

