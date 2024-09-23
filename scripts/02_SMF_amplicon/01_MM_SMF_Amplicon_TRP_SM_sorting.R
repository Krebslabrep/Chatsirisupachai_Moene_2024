########## Data preparation for mouse (SMF amplicon SM sorting for TRP analysis) merged replicates ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(QuasR)
library(SingleMoleculeFootprinting)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(ggplot2)
library(stringi)

### amplicon regions
amplicons <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_amplicon_regions.rds")
amplicons <- amplicons[grep("^NM|^NR", amplicons$TFBSname)]

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Qinput
# this file points toward BAM files of amplicon SMF in TKO mESCs from this study (E-MTAB-14461)
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_amplicon_SMF_MM_TRP_experiment.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
TSSsc_of_interest <- TSSsc[unlist(TSSsc$gene_id) %in% unique(amplicons$TFBSname)]


### Sort reads around promoter using SingleMoleculeFootprinting package
res_list <- SortReadsByPromoter_MultiSiteWrapper(sampleSheet = Qinput, 
                                                 sample = MySample, 
                                                 genome = BSgenome.Mmusculus.UCSC.mm10, 
                                                 coverage = 20, 
                                                 ConvRate.thr = 0.2, 
                                                 TSSsc = TSSsc_of_interest, 
                                                 max_interTF_distance = 100000, 
                                                 max_window_width = 5000000, 
                                                 min_cluster_width = 600, 
                                                 sorting_coverage = 20, cores = 20)

saveRDS(res_list, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_TKO_DE_amplicon_TRP_sorting_output.rds")


