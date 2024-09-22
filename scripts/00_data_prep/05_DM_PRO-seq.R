########## Data preparation for Drosophila (PRO-seq) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(QuasR)
library(dplyr)

##### CAGE-corrected TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
# get promoter regions
promReg <- promoters(TSSsc, upstream = 200, downstream = 100)

##### Qinput file
# this file points toward BAM files of PRO-seq of S2 cells from Kwak et al., 2013
input_file <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_PRO-seq_DM.txt"

##### function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

##### PRO-seq
PRO_seq_DM <- qAlign(sampleFile = input_file, 
                     genome = "BSgenome.Dmelanogaster.UCSC.dm6", 
                     paired = "no",
                     checkOnly = T)

cluObj <- makeCluster(10)

align_stats <- as.data.frame(alignmentStats(PRO_seq_DM))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

count_promoters <- as.data.frame(qCount(PRO_seq_DM, promReg, selectReadPosition = "start", orientation = "opposite"))
count_promoters$width <- NULL

##### normalised count promoters
count_promoters$PROseq_S2_SRX203291 <- unlist(lapply(count_promoters$PROseq_S2_SRX203291, RPM, all_mapped = align_stats[rownames(align_stats) == "PROseq_S2_SRX203291",]$mapped))
count_promoters$PROseq_S2_SRX203292 <- unlist(lapply(count_promoters$PROseq_S2_SRX203292, RPM, all_mapped = align_stats[rownames(align_stats) == "PROseq_S2_SRX203292",]$mapped))
colnames(count_promoters) <- c("S2_R1", "S2_R2")
count_promoters$S2_avg <- (count_promoters$S2_R1 + count_promoters$S2_R2) / 2

saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_PRO-seq_RPM.rds")


##### PRO-seq profile
profile_counts <- qProfile(PRO_seq_DM, TSSsc, upstream = 500, downstream = 500, selectReadPosition = "start",
                           orientation = "opposite", clObj = cluObj)
profile_counts_merged <- profile_counts$PROseq_S2_SRX203291 + profile_counts$PROseq_S2_SRX203292

profile_counts_merged_RPM <- RPM(profile_counts_merged, all_mapped = sum(align_stats$mapped))

#smoothened PRO-seq data
smw <- 20
PRO_seq_sm <- apply(profile_counts_merged_RPM, 1, function(x){caTools::runmean(x, smw, endrule = "constant")})
PRO_seq_smt <- t(PRO_seq_sm)
colnames(PRO_seq_smt) <- -500:500
saveRDS(PRO_seq_smt, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_PRO_seq_profile_RPM_smoothed.rds")

