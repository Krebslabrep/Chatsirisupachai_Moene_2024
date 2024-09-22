########## Data preparation for Drosophila (TSSs) ##########
# Author: Kasit Chatsirisupachai (modified from Christine Moene's code)
# LastUpdate: 23.09.2024

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(parallel)
library(RMariaDB)
library(GenomicFeatures)
library(Biostrings)
library(AnnotationDbi)
library(rtracklayer)
library(liftOver)

##### Load AllC (position of all cytosines )
AllC <- readRDS("/g/krebs/moene/analysis/SMF/AllC_dm6_2.rds")

##### define genes (TSSs)
rs_tss <- makeTxDbFromUCSC(genome = "dm6", tablename = "ncbiRefSeq")
REFSEQ <- transcripts(rs_tss)
saveRDS(REFSEQ, "/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/data/DM/DM_REFSEQ_reference_transcripts.rds")

genelengths_DM <- width(transcripts(rs_tss))
names(genelengths_DM) <- transcripts(rs_tss)$tx_name

TSSs <- resize(transcripts(rs_tss), 1, fix = "start")

### re-annotate TSS with CAGE
# because it was originally aligned to dm3, we used liftOver
dm3_dm6_chain <- import("/g/krebs/moene/DB/DM/dm3ToDm6.over.chain")

## forward data
CAGEfw <- import("/g/krebs/moene/DB/DM/modEncode_2549/signal_data_files/2549_celniker_cage_fwd.wig")
CAGEfw_dm6_list <- liftOver(CAGEfw, dm3_dm6_chain)
summary(lengths(CAGEfw_dm6_list)) # all are 1 to 0 or 1 to 1, so we can collapse the list to a GRanges object
CAGEfw_dm6 <- unlist(CAGEfw_dm6_list)

## reverse (backwards) data
CAGEbkwd <- import("/g/krebs/moene/DB/DM/modEncode_2549/signal_data_files/2549_celniker_cage_bkwd.wig")
CAGEbkwd_dm6_list <- liftOver(CAGEbkwd, dm3_dm6_chain)
summary(lengths(CAGEbkwd_dm6_list)) #all are 1 to 0 or 1 to 1, so we can collapse the list to a GRanges object
CAGEbkwd_dm6 <- unlist(CAGEbkwd_dm6_list)

# compare CAGE data to TSSs of genes
TSSs <- resize(transcripts(rs_tss), 1, fix = "start")
ov1 <- as.matrix(findOverlaps(resize(TSSs, 100, fix = "center"), CAGEfw_dm6))
ov2 <- as.matrix(findOverlaps(resize(TSSs, 100, fix = "center"), CAGEbkwd_dm6))

TSS_START_CAGE <- mclapply(seq_along(TSSs),function(i){
  if(as.character(strand(TSSs[i])) == "+"){ # + strand genes compare to CAGEfw
    if(i %in% ov1[, 1]){
      CAGE <- CAGEfw_dm6[ov1[(ov1[, 1] == i), 2]] #takes CAGE_main subject hit(s) belonging to queryHit==i
      start(CAGE)[order(CAGE$score, decreasing = T)][1] #if multiple CAGE TSSs, take the one with highest CAGE_score
    } else(start(TSSs[i]))
  } else{                                    # - strand genes compare to CAGEbkwd
    if(i %in% ov2[,1]){
      CAGE <- CAGEbkwd_dm6[ov2[(ov2[, 1] == i), 2]] #takes CAGE_main subject hit(s) belonging to queryHit==i
      start(CAGE)[order(CAGE$score,decreasing = F)][1] #if multiple CAGE TSSs, take the one with highest (= most negative) CAGE_score
    } else(start(TSSs[i]))
  }
}, mc.cores=10) #makes new list of TSS, with either the original TSS or replaced by TSS from CAGE

TSSsc_dm6 <- GRanges(seqnames(TSSs), IRanges(unlist(TSS_START_CAGE), unlist(TSS_START_CAGE)), strand(TSSs))
elementMetadata(TSSsc_dm6) <- elementMetadata(TSSs)
names(TSSsc_dm6) <- names(TSSs)
seqlengths(TSSsc_dm6) <- seqlengths(TSSs)
seqinfo(TSSsc_dm6) <- seqinfo(TSSs)

length(TSSsc_dm6) # 34463
length(unique(TSSsc_dm6)) # 22323

TSSsc_dm6.2 <- unique(TSSsc_dm6)
norm_chr <- seqnames(Dmelanogaster)[1:8]
TSSsc_dm6.3 <- TSSsc_dm6.2[seqnames(TSSsc_dm6.2) %in% norm_chr]
names(TSSsc_dm6.3) <- TSSsc_dm6.3$tx_name
TSSsc_dm6.3$gene_id <- TSSsc_dm6.3$tx_name

TSSsc_dm6.3 <- TSSsc_dm6.3[seqnames(TSSsc_dm6.3) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4")] # only keep autosomes
length(TSSsc_dm6.3)  # 18616

saveRDS(TSSsc_dm6.3, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

##### TATA promoters
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
st <- resize(TSSsc, 37, fix = "end")
st <- resize(st, 18, fix = "start")
st.seq <- getSeq(Dmelanogaster, st)
tbpmot1b <- vcountPattern(DNAString("TATAWAWR"), st.seq, max.mismatch = 1, fixed = "subject") > 0

TATA_promoters <- data.frame(gene_id = names(st.seq), TATA = tbpmot1b)
saveRDS(TATA_promoters, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_TATA_promoters.rds")
