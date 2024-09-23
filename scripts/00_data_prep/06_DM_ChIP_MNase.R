########## Data preparation for Drosophila (ChIP-seq & MNase-seq) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(QuasR)
library(dplyr)

##### CAGE-corrected TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
# get promoter regions
promReg <- promoters(TSSsc, upstream = 200, downstream = 100)

##### function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

##### Pol II ChIP-seq for Drosophila S2 cells
### Qinput file
# this file points toward BAM files of Pol II ChIP-seq in S2 cells from Tettey et al., 2019
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_ChIP-seq_DM.txt"
PolII_proj <- qAlign(Qinput,
                     genome = "BSgenome.Dmelanogaster.UCSC.dm6",
                     paired = "no")
align_stats <- as.data.frame(alignmentStats(PolII_proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

count_promoters <- as.data.frame(qCount(PolII_proj, promReg, selectReadPosition = "start", shift = 75, clObj = 8))
count_promoters$width <- NULL

saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_counts.rds")

# normalise count promoters
count_promoters$PolII_S2_SRX3981728 <- unlist(lapply(count_promoters$PolII_S2_SRX3981728, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_S2_SRX3981728",]$mapped))
count_promoters$PolII_S2_SRX3981727 <- unlist(lapply(count_promoters$PolII_S2_SRX3981727, RPM, all_mapped = align_stats[rownames(align_stats) == "PolII_S2_SRX3981727",]$mapped))

count_promoters$PolII_avg <- (count_promoters$PolII_S2_SRX3981728 + count_promoters$PolII_S2_SRX3981727) / 2
saveRDS(count_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_RPM.rds")

### Select top X% promoters
prom_QT <- quantile(count_promoters$PolII_avg, seq(0, 1, 0.01))
top5_promoters <- rownames(count_promoters[count_promoters$PolII_avg >= prom_QT[96],])
top10_promoters <- rownames(count_promoters[count_promoters$PolII_avg >= prom_QT[91],])
top20_promoters <- rownames(count_promoters[count_promoters$PolII_avg >= prom_QT[81],])
top50_promoters <- rownames(count_promoters[count_promoters$PolII_avg >= prom_QT[51],])
bottom10_promoters <- rownames(count_promoters[count_promoters$PolII_avg <= prom_QT[10],])
bottom20_promoters <- rownames(count_promoters[count_promoters$PolII_avg <= prom_QT[20],])

prom_list <- list(top5 = top5_promoters, 
                  top10 = top10_promoters,
                  top20 = top20_promoters,
                  top50 = top50_promoters,
                  bottom10 = bottom10_promoters,
                  bottom20 = bottom20_promoters)

saveRDS(prom_list,"/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_ChIP-seq_quantile_list.rds")


##### MNase-seq for S2 cells
# this file points toward BAM files of MNase-seq in S2 cells from Gilchrist et al., 2010
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_MNase-seq_DM.txt"
MNase_proj <- qAlign(Qinput,
                     genome = "BSgenome.Dmelanogaster.UCSC.dm6",
                     paired = "no")
align_stats <- as.data.frame(alignmentStats(MNase_proj))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

# count on specific regions
regions <- promoters(TSSsc, upstream = 40, downstream = 30)

MNase_count_TSS <- as.data.frame(qCount(MNase_proj, regions, selectReadPosition = "start", shift = 80, orientation = "any", clObj = 8))
MNase_count_TSS$width <- NULL

# normalise MNase counts
MNase_count_TSS$MNase_S2_SRX021645 <- unlist(lapply(MNase_count_TSS$MNase_S2_SRX021645, RPM, all_mapped = align_stats[rownames(align_stats) == "MNase_S2_SRX021645", "mapped"]))
MNase_count_TSS$MNase_S2_SRX021646 <- unlist(lapply(MNase_count_TSS$MNase_S2_SRX021646, RPM, all_mapped = align_stats[rownames(align_stats) == "MNase_S2_SRX021646", "mapped"]))
MNase_count_TSS$MNase_S2_SRX021647 <- unlist(lapply(MNase_count_TSS$MNase_S2_SRX021647, RPM, all_mapped = align_stats[rownames(align_stats) == "MNase_S2_SRX021647", "mapped"]))

# average replicates
MNase_count_TSS$merged_reps <- (MNase_count_TSS$MNase_S2_SRX021645 + MNase_count_TSS$MNase_S2_SRX021646 + MNase_count_TSS$MNase_S2_SRX021647) / 3
saveRDS(MNase_count_TSS, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_MNase_RPM.rds")


# MNase profile
profile_counts <- qProfile(MNase_proj, TSSsc, upstream = 500, downstream = 500, selectReadPosition = "start", 
                           shift = 80, clObj = 8)

profile_counts_merged <- profile_counts$MNase_S2_SRX021645 + profile_counts$MNase_S2_SRX021646 + profile_counts$MNase_S2_SRX021647

profile_counts_merged_RPM <- RPM(profile_counts_merged, all_mapped = sum(align_stats$mapped))


#smoothened MNase data
smw <- 20
MNase_sm <- apply(profile_counts_merged_RPM, 1, function(x){caTools::runmean(x, smw, endrule = "constant")})
MNase_smt <- t(MNase_sm)
colnames(MNase_smt) <- -500:500

saveRDS(MNase_smt, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_promoter_MNase_profile_RPM_smoothed.rds")

