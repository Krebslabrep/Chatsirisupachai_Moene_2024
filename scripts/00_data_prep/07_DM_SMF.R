########## Data preparation for Drosophila (SMF) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(BSgenome)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(QuasR)
library(caTools)
library(plyranges)

### Qinput
# this file points toward BAM files of genome-wide SMF in S2 cells from Krebs et al., 2017
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_whole-genome_SMF_DM.txt"
MySample <- unique(suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]]))

### Create methylation calling window
# Partition a genome into smaller regions
dm_length <- seqlengths(Dmelanogaster)
chroms <- tileGenome(dm_length, tilewidth = max(dm_length), cut.last.tile.in.chrom = TRUE)
regs <- tileGenome(dm_length, tilewidth=3000000, cut.last.tile.in.chrom = TRUE)
norm_chr <- seqnames(chroms)[1:8]

MethylationCallingWindows <- regs

### Call context methylation
# function to call context methylation in each calling window
Context_methylation <- function(i, sampleSheet = Qinput, 
                                sample = MySample, 
                                genome = BSgenome.Dmelanogaster.UCSC.dm6, 
                                coverage = 20, 
                                ConvRate.thr = NULL){
  
  # select methylation calling window
  CurrentWindow <- MethylationCallingWindows[i]
  
  # single enzyme (NO) or double enzyme (DE) experiment (in this case is the double enzyme)
  ExperimentType <- suppressMessages(SingleMoleculeFootprinting::DetectExperimentType(Samples = sample))
  
  CallContextMethylation(sampleSheet = sampleSheet,
                         sample = sample,
                         genome = genome,
                         RegionOfInterest = CurrentWindow,
                         coverage = coverage,
                         ConvRate.thr = ConvRate.thr,
                         returnSM = FALSE) -> Methylation
  return(Methylation)
}

# apply the function
ContextMeth <- parallel::mclapply(seq_along(MethylationCallingWindows), Context_methylation,
                                  sampleSheet = Qinput,
                                  sample = MySample,
                                  genome = BSgenome.Dmelanogaster.UCSC.dm6,
                                  coverage = 20,
                                  ConvRate.thr = NULL,
                                  mc.cores = 20)

ContextMeth_all <- do.call(c, ContextMeth)

# saveRDS(ContextMeth_all, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_context_methylation.rds")
# this object is too large to upload
# subset only cytosines around promoters
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSS_extended <- resize(TSSsc_DM, width = 1000, fix = "center")

ContextMeth_all %>% filter_by_overlaps(TSS_extended) -> ContextMeth_subset
saveRDS(ContextMeth_subset, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_context_methylation.rds")
