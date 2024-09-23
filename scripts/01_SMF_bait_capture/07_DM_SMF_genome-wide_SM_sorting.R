########## Data preparation for Drosophila (SMF whole-genome SM sorting) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(stringi)
library(SingleMoleculeFootprinting)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)


### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Qinput
# this file points toward BAM files of whole-genome SMF in S2 cells from Krebs et al., 2017
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_whole-genome_SMF_DM.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
MySample <- unique(MySample)

### Sort reads around promoter using SingleMoleculeFootprinting package
res_list <- SortReadsByPromoter_MultiSiteWrapper(sampleSheet = Qinput, 
                                                 sample = MySample, 
                                                 genome = BSgenome.Dmelanogaster.UCSC.dm6, 
                                                 coverage = 20, 
                                                 ConvRate.thr = NULL, 
                                                 TSSsc = TSSsc, 
                                                 max_interTF_distance = 50000, 
                                                 max_window_width = 1200000, 
                                                 min_cluster_width = 600, 
                                                 sorting_coverage = 20, cores = 36)

saveRDS(res_list, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_sorting_output.rds")

