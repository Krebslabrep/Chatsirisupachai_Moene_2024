########## Data preparation for mouse (SMF bait-capture SM sorting) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(SingleMoleculeFootprinting)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(stringi)

### baits regions
baits_mm10 <- import("/g/krebs/krebs/DB/SureSelect/MouseMethyl_Bait_merged_mm10.bed")
baits_mm10_ext <- GRanges(seqnames(baits_mm10), IRanges(start(baits_mm10) - 500, end(baits_mm10) + 500))

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Qinput
# this file points toward BAM files of bait-capture SMF in TKO mESCs from Sonmezer et al., 2021
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_SMF_bait_capture_MM.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])


### Sort reads around promoter using SingleMoleculeFootprinting package
res_list <- SortReadsByPromoter_MultiSiteWrapper(sampleSheet = Qinput, 
                                                 sample = MySample, 
                                                 genome = BSgenome.Mmusculus.UCSC.mm10, 
                                                 coverage = 20, 
                                                 ConvRate.thr = NULL, 
                                                 TSSsc = TSSsc, 
                                                 max_interTF_distance = 100000, 
                                                 max_window_width = 5000000, 
                                                 min_cluster_width = 600, 
                                                 sorting_coverage = 20, cores = 20)

saveRDS(res_list, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/bait_capture_SMF/MM_TKO_DE_bait_capture_sorting_output.rds")


