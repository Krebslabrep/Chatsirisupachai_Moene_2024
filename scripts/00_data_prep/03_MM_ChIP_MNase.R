########## Data preparation for mouse (ChIP-seq & MNase-seq) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(QuasR)
library(dplyr)

##### CAGE-corrected TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
# get promoter regions
promReg <- promoters(TSSsc, upstream = 200, downstream = 100)

##### function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

##### Pol II ChIP-seq for mESCs
### Qinput file
# this file points toward BAM files of Pol II ChIP-seq in mESCs from Langer et al., 2016
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_ChIP-seq_MM.txt"

ChIP_proj <- qAlign(Qinput,
                    genome = "BSgenome.Mmusculus.UCSC.mm10",
                    paired = "no")

align_stats <- as.data.frame(alignmentStats(ChIP_proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

count_promoters <- as.data.frame(qCount(ChIP_proj, promReg, selectReadPosition = "start", shift = 75, clObj = 4))
count_promoters$width <- NULL

saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_counts.rds")

# normalise count promoters
count_promoters$input_mESC_SRX032475 <- unlist(lapply(count_promoters$input_mESC_SRX032475, RPM, all_mapped = align_stats[rownames(align_stats) == "input_mESC_SRX032475",]$mapped))
count_promoters$PolII_mESC_SRX1089844 <- unlist(lapply(count_promoters$PolII_mESC_SRX1089844, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_mESC_SRX1089844",]$mapped))
count_promoters$PolII_enrichment <- count_promoters$PolII_mESC_SRX1089844 - count_promoters$input_mESC_SRX032475

saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_RPM.rds")

### Select top X% promoters
prom_QT <- quantile(count_promoters$PolII_enrichment, seq(0, 1, 0.01))
top5_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment >= prom_QT[96],])
top10_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment >= prom_QT[91],])
top20_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment >= prom_QT[81],])
top50_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment >= prom_QT[51],])
bottom10_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment <= prom_QT[10],])
bottom20_promoters <- rownames(count_promoters[count_promoters$PolII_enrichment <= prom_QT[20],])

prom_list <- list(top5 = top5_promoters, 
                  top10 = top10_promoters,
                  top20 = top20_promoters,
                  top50 = top50_promoters,
                  bottom10 = bottom10_promoters,
                  bottom20 = bottom20_promoters)

saveRDS(prom_list,"/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")


##### Pol II ChIP-seq for other cell types
### Qinput
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_ChIP_PolII_mouse_cells.txt"
PolII_proj <- qAlign(Qinput,
                     genome = "BSgenome.Mmusculus.UCSC.mm10",
                     paired = "no")

align_stats <- as.data.frame(alignmentStats(PolII_proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

count_promoters <- as.data.frame(qCount(PolII_proj, promReg, selectReadPosition = "start", shift = 75, clObj = 8))
count_promoters$width <- NULL

# normalise count promoters
count_promoters$PolII_C2C12_SRX062102 <- unlist(lapply(count_promoters$PolII_C2C12_SRX062102, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_C2C12_SRX062102",]$mapped))
count_promoters$PolII_MEL_SRX140358 <- unlist(lapply(count_promoters$PolII_MEL_SRX140358, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_MEL_SRX140358",]$mapped))
count_promoters$PolII_NP_SRX032484 <- unlist(lapply(count_promoters$PolII_NP_SRX032484, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_NP_SRX032484",]$mapped))
count_promoters$PolII_mESC_SRX1089844 <- unlist(lapply(count_promoters$PolII_mESC_SRX1089844, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_mESC_SRX1089844",]$mapped))

### Select top 5% promoters
prom_QT <- quantile(count_promoters$PolII_C2C12_SRX062102, seq(0, 1, 0.01))
top5_C2C12 <- rownames(count_promoters[count_promoters$PolII_C2C12_SRX062102 >= prom_QT[96],])

prom_QT <- quantile(count_promoters$PolII_MEL_SRX140358, seq(0, 1, 0.01))
top5_MEL <- rownames(count_promoters[count_promoters$PolII_MEL_SRX140358 >= prom_QT[96],])

prom_QT <- quantile(count_promoters$PolII_NP_SRX032484, seq(0, 1, 0.01))
top5_NP<- rownames(count_promoters[count_promoters$PolII_NP_SRX032484 >= prom_QT[96],])

prom_QT <- quantile(count_promoters$PolII_mESC_SRX1089844, seq(0, 1, 0.01))
top5_mESC<- rownames(count_promoters[count_promoters$PolII_mESC_SRX1089844 >= prom_QT[96],])

prom_list <- list(top5_C2C12 = top5_C2C12, 
                  top5_MEL = top5_MEL,
                  top5_NP = top5_NP,
                  top5_mESC = top5_mESC)

saveRDS(prom_list,"/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_top_5_list_other_cell_types.rds")



##### MNase-seq for mESCs
### Qinput
# this file points toward BAM files of MNase-seq in mESCs from Barisic et al., 2019
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_mESCs_MNase-seq.txt"

MNase_proj <- qAlign(Qinput,
                     genome = "BSgenome.Mmusculus.UCSC.mm10",
                     paired = "no")
align_stats <- as.data.frame(alignmentStats(MNase_proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

# count on specific regions
regions <- promoters(TSSsc, upstream = 40, downstream = 30)

MNase_count_TSS <- as.data.frame(qCount(MNase_proj, regions, selectReadPosition = "start", shift = 80, orientation = "any", clObj = 8))
MNase_count_TSS$width <- NULL

# normalise MNase counts
MNase_count_TSS$Mnase_mESC_WT_Rep1_SRX3827078 <- unlist(lapply(MNase_count_TSS$Mnase_mESC_WT_Rep1_SRX3827078, RPM, all_mapped = align_stats[rownames(align_stats) == "Mnase_mESC_WT_Rep1_SRX3827078", "mapped"]))
MNase_count_TSS$Mnase_mESC_WT_Rep2_SRX3827079 <- unlist(lapply(MNase_count_TSS$Mnase_mESC_WT_Rep2_SRX3827079, RPM, all_mapped = align_stats[rownames(align_stats) == "Mnase_mESC_WT_Rep2_SRX3827079", "mapped"]))

# average replicates
MNase_count_TSS$merged_reps <- (MNase_count_TSS$Mnase_mESC_WT_Rep1_SRX3827078 + MNase_count_TSS$Mnase_mESC_WT_Rep2_SRX3827079) / 2
saveRDS(MNase_count_TSS, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_MNase_RPM.rds")


# MNase profile
profile_counts <- qProfile(MNase_proj, TSSsc, upstream = 500, downstream = 500, selectReadPosition = "start", 
                           shift = 80, clObj = 8)

profile_counts_merged <- profile_counts$Mnase_mESC_WT_Rep1_SRX3827078 + profile_counts$Mnase_mESC_WT_Rep2_SRX3827079

profile_counts_merged_RPM <- RPM(profile_counts_merged, all_mapped = sum(filter(align_stats, grepl("Mnase_mESC_WT", rownames(align_stats)))$mapped))


#smoothened MNase data
smw <- 20
MNase_sm <- apply(profile_counts_merged_RPM, 1, function(x){caTools::runmean(x, smw, endrule = "constant")})
MNase_smt <- t(MNase_sm)
colnames(MNase_smt) <- -500:500

saveRDS(MNase_smt, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_MNase_profile_RPM_smoothed.rds")

