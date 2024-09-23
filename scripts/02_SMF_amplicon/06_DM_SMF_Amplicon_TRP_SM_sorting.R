########## Data preparation for Drosophila (SMF amplicon SM sorting for TRP analysis 0 min vs 20 min) merged replicates ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(QuasR)
library(SingleMoleculeFootprinting)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)
library(ggplot2)
library(stringi)

### amplicon regions (since I don't have the regions Granges object or the promoters, I try calling all TSSs and will get those that pass the filter)

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Qinput
# this file points toward BAM files of amplicon SMF in S2 cells from Krebs et al., 2017
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_amplicon_SMF_DM_TRP_experiment.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])


### Sort reads around promoter using SingleMoleculeFootprinting package
res_list <- SortReadsByPromoter_MultiSiteWrapper(sampleSheet = Qinput, 
                                                 sample = MySample, 
                                                 genome = BSgenome.Dmelanogaster.UCSC.dm6, 
                                                 coverage = 20, 
                                                 ConvRate.thr = 0.2, 
                                                 TSSsc = TSSsc, 
                                                 max_interTF_distance = 50000, 
                                                 max_window_width = 1200000, 
                                                 min_cluster_width = 600, 
                                                 sorting_coverage = 20, cores = 1)

res_list[[3]] %>% filter(complete.cases(.)) %>% distinct(TSS) %>% saveRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/amplicon_SMF/DM_amplicon_TSSs.rds")
saveRDS(res_list, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/amplicon_SMF/DM_S2_DE_amplicon_TRP_sorting_output.rds")


