########## Data preparation for mouse (TSSs) ##########
# Author: Kasit Chatsirisupachai & Christine Moene
# LastUpdate: 22.09.2024

library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)

##### define genes (TSSs)
REFSEQ <- readRDS("/g/krebs/krebs/analysis/RNA-seq/rds/REFSEQ_reference_transcripts.rds")

TSSs <- promoters(REFSEQ, 1, 1)
TSSs <- TSSs[seqnames(TSSs) %in% paste("chr", seq(19), sep = "")] #only keep chr 1-19

TSSs <- unique(TSSs)
names(TSSs) <- TSSs$gene_id

# re-annotate TSS with CAGE
CAGE_peaks <- import("/g/krebs/krebs/DB/CAGE/mm10_liftover+new_CAGE_peaks_phase1and2.bed.gz")
CAGE_main <- GRanges(seqnames(CAGE_peaks), CAGE_peaks$thick, strand = strand(CAGE_peaks), score = CAGE_peaks$score)

o2 <- as.matrix(findOverlaps(resize(TSSs, 100, fix = "center"), CAGE_main))

TSS_START_CAGE <- mclapply(seq_along(TSSs), function(i){
  if(i %in% o2[, 1]){
    CAGE <- CAGE_main[o2[(o2[, 1] == i), 2]] #takes CAGE_main subject hit(s) belonging to queryHit==i
    start(CAGE)[order(CAGE$score, decreasing = T)][1] #if multiple CAGE TSSs, take the one with highest CAGE_score
  } else (start(TSSs[i]))
}, mc.cores=10) #makes new list of TSS, with either the original TSS or replaced by TSS from CAGE

TSSsc <- GRanges(seqnames(TSSs), IRanges(unlist(TSS_START_CAGE), unlist(TSS_START_CAGE)), strand(TSSs))
elementMetadata(TSSsc) <- elementMetadata(TSSs)
names(TSSsc) <- names(TSSs)
seqlengths(TSSsc) <- seqlengths(TSSs)
seqinfo(TSSsc) <- seqinfo(TSSs)
saveRDS(TSSsc, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

##### TATA promoters
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc
st <- resize(TSSsc, 37, fix = "end")
st <- resize(st, 18, fix = "start")
st.seq <- getSeq(Mmusculus, st)
tbpmot1b <- vcountPattern(DNAString("TATAWAWR"), st.seq, max.mismatch = 1, fixed = "subject") > 0

TATA_promoters <- data.frame(gene_id = names(st.seq), TATA = tbpmot1b)
saveRDS(TATA_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")

