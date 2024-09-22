########## Data preparation for mouse (PRO-seq) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 22.09.2024

library(QuasR)
library(dplyr)

##### CAGE-corrected TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
# get promoter regions
promReg <- promoters(TSSsc, upstream = 200, downstream = 100)

##### Qinput file
# this file points toward BAM files of PRO-seq in TKO mESCs from Kreibich et al., 2023
input_file <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_PRO-seq_MM_Kreibich2023.txt"

##### function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

##### PRO-seq for TKO and WT mESCs
PRO_seq_MM <- qAlign(sampleFile = input_file, 
                     genome = "BSgenome.Mmusculus.UCSC.mm10", 
                     paired = "no",
                     checkOnly = T)

cluObj <- makeCluster(10)

align_stats <- as.data.frame(alignmentStats(PRO_seq_MM))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

count_promoters <- as.data.frame(qCount(PRO_seq_MM, promReg, selectReadPosition = "start", orientation = "opposite", clObj = cluObj))
count_promoters$width <- NULL

##### normalised count promoters
count_promoters$TKO_R1 <- unlist(lapply(count_promoters$TKO_R1, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_R1",]$mapped))
count_promoters$TKO_R2 <- unlist(lapply(count_promoters$TKO_R2, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_R2",]$mapped))
count_promoters$TKO_R3 <- unlist(lapply(count_promoters$TKO_R3, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_R3",]$mapped))
count_promoters$TKO_avg <- as.numeric(apply(dplyr::select(count_promoters, contains("TKO")), 1, mean))

saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_Kreibich2023_PRO-seq_RPM.rds")


##### PRO-seq profile
profile_counts <- qProfile(PRO_seq_MM, TSSsc, upstream = 500, downstream = 500, selectReadPosition = "start",
                           orientation = "opposite", clObj = cluObj)
profile_counts_merged <- profile_counts$TKO_R1 + profile_counts$TKO_R2 + profile_counts$TKO_R3

profile_counts_merged_RPM <- RPM(profile_counts_merged, all_mapped = sum(filter(align_stats, grepl("TKO", rownames(align_stats)))$mapped))

#smoothened PRO-seq data
smw <- 20
PRO_seq_sm <- apply(profile_counts_merged_RPM, 1, function(x){caTools::runmean(x, smw, endrule = "constant")})
PRO_seq_smt <- t(PRO_seq_sm)
colnames(PRO_seq_smt) <- -500:500
saveRDS(PRO_seq_smt, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_Kreibich2023_PRO_seq_profile_RPM_smoothed.rds")
