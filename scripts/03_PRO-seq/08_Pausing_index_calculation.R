########## Pausing Index calculation for mouse and Drosophila ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.01.2025

library(QuasR)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Dmelanogaster.UCSC.dm6)

##### function to calculate RPM
RPM <- function(read_counts, all_mapped){
  norm_counts <- (1e6 / all_mapped) * read_counts
  return(norm_counts)
}

##### Mouse #####
### mouse TSSs
TSSsc_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc_MM$gene_id <- unlist(TSSsc_MM$gene_id)

genes_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected_whole_gene.rds")

# keep only genes that are longer than 600 bases
genes_MM <- genes_MM[width(genes_MM) >= 600, ]
length(genes_MM) # 24869 genes

##### Calculate pausing index (PI) for mouse
# Load PRO-seq
# this file points to the BAM files of PRO-seq in TKO mESC cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_MM_TRP_PRO-seq.txt"
genome <- "BSgenome.Mmusculus.UCSC.mm10"
paired <- "no"

### PRO-seq for TKO mESCs
PRO_seq_MM <- qAlign(sampleFile = Qinput_sample, 
                     genome = "BSgenome.Mmusculus.UCSC.mm10", 
                     paired = "no",
                     checkOnly = T)

cluObj <- makeCluster(4)

align_stats <- as.data.frame(alignmentStats(PRO_seq_MM))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

### Count reads from promoter regions
promReg <- promoters(genes_MM, upstream = 150, downstream = 150)
count_promoters <- as.data.frame(qCount(PRO_seq_MM, promReg, selectReadPosition = "start", orientation = "opposite", clObj = cluObj))

# keep only samples from 0 min
count_promoters %>%
  select(TKO_0min_R1, TKO_0min_R2) -> count_promoters

count_promoters$TKO_0min_R1 <- unlist(lapply(count_promoters$TKO_0min_R1, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_0min_R1",]$mapped))
count_promoters$TKO_0min_R2 <- unlist(lapply(count_promoters$TKO_0min_R2, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_0min_R2",]$mapped))
count_promoters$TKO_avg <- (count_promoters$TKO_0min_R1 + count_promoters$TKO_0min_R2)/2


### Count reads from gene bodies
bodyReg <- resize(GenomicRanges::shift(TSSsc_MM, 300), 300, fix = "start")
bodyReg <- bodyReg[names(bodyReg) %in% genes_MM$gene_id]
count_gene_bodies <- as.data.frame(qCount(PRO_seq_MM, bodyReg, selectReadPosition = "start", orientation = "opposite", clObj = cluObj))

# keep only samples from 0 min
count_gene_bodies %>%
  select(TKO_0min_R1, TKO_0min_R2) -> count_gene_bodies

count_gene_bodies$TKO_0min_R1 <- unlist(lapply(count_gene_bodies$TKO_0min_R1, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_0min_R1",]$mapped))
count_gene_bodies$TKO_0min_R2 <- unlist(lapply(count_gene_bodies$TKO_0min_R2, RPM, all_mapped = align_stats[rownames(align_stats) == "TKO_0min_R2",]$mapped))
count_gene_bodies$TKO_avg <- (count_gene_bodies$TKO_0min_R1 + count_gene_bodies$TKO_0min_R2)/2


### Calculate Pause Index
PI_MM <- data.frame(genes = genes_MM$gene_id, 
                    TSS_count = count_promoters$TKO_avg,
                    GB_count = count_gene_bodies$TKO_avg,
                    PI = count_promoters$TKO_avg / count_gene_bodies$TKO_avg)

saveRDS(PI_MM, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_qPRO-seq_PI_GB_TSS300-TSS600.rds")



##### Drosophila #####
### Drosophila TSSs
TSSsc_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc_DM$gene_id <- unlist(TSSsc_DM$gene_id)

genes_DM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected_whole_gene.rds")

# keep only genes that are longer than 600 bases
genes_DM <- genes_DM[width(genes_DM) >= 600, ]
length(genes_DM) # 16198 genes

##### Calculate pausing index (PI) for Drosophila
# Load PRO-seq
# this file points to the BAM files of PRO-seq in Drosophila S2 cells following a time-course TRP treatment from this study (E-MTAB-14462)
Qinput_sample <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_DM_TRP_PRO-seq.txt"
genome <- "BSgenome.Dmelanogaster.UCSC.dm6"
paired <- "no"

### PRO-seq for Drosophila S2
PRO_seq_DM <- qAlign(sampleFile = Qinput_sample, 
                     genome = "BSgenome.Dmelanogaster.UCSC.dm6", 
                     paired = "no",
                     checkOnly = T)

cluObj <- makeCluster(4)

align_stats <- as.data.frame(alignmentStats(PRO_seq_DM))
rownames(align_stats) <- unlist(lapply(rownames(align_stats), gsub, pattern = ":genome", replacement = ""))

### Count reads from promoter regions
promReg <- promoters(genes_DM, upstream = 150, downstream = 150)
count_promoters <- as.data.frame(qCount(PRO_seq_DM, promReg, selectReadPosition = "start", orientation = "opposite", clObj = cluObj))

# keep only samples from 0 min
count_promoters %>%
  select(S2_0min_R1, S2_0min_R2) -> count_promoters

count_promoters$S2_0min_R1 <- unlist(lapply(count_promoters$S2_0min_R1, RPM, all_mapped = align_stats[rownames(align_stats) == "S2_0min_R1",]$mapped))
count_promoters$S2_0min_R2 <- unlist(lapply(count_promoters$S2_0min_R2, RPM, all_mapped = align_stats[rownames(align_stats) == "S2_0min_R2",]$mapped))
count_promoters$S2_avg <- (count_promoters$S2_0min_R1 + count_promoters$S2_0min_R2)/2


### Count reads from gene bodies
bodyReg <- resize(GenomicRanges::shift(TSSsc_DM, 300), 300, fix = "start")
bodyReg <- bodyReg[names(bodyReg) %in% genes_DM$gene_id]
count_gene_bodies <- as.data.frame(qCount(PRO_seq_DM, bodyReg, selectReadPosition = "start", orientation = "opposite", clObj = cluObj))

# keep only samples from 0 min
count_gene_bodies %>%
  select(S2_0min_R1, S2_0min_R2) -> count_gene_bodies

count_gene_bodies$S2_0min_R1 <- unlist(lapply(count_gene_bodies$S2_0min_R1, RPM, all_mapped = align_stats[rownames(align_stats) == "S2_0min_R1",]$mapped))
count_gene_bodies$S2_0min_R2 <- unlist(lapply(count_gene_bodies$S2_0min_R2, RPM, all_mapped = align_stats[rownames(align_stats) == "S2_0min_R2",]$mapped))
count_gene_bodies$S2_avg <- (count_gene_bodies$S2_0min_R1 + count_gene_bodies$S2_0min_R2)/2


### Calculate Pause Index
PI_DM <- data.frame(genes = genes_DM$gene_id, 
                    TSS_count = count_promoters$S2_avg,
                    GB_count = count_gene_bodies$S2_avg,
                    PI = count_promoters$S2_avg / count_gene_bodies$S2_avg)

saveRDS(PI_DM, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_qPRO-seq_PI_GB_TSS300-TSS600.rds")
