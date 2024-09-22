########## Composite plot for mouse SMF whole-genome data in TKO ES and other cells ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 22.09.2024

library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(plyranges)

### Load data
TSSsc_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TATA_promoters_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_TATA_promoters.rds")
Top_promoters_MM_PolII_ChIP <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_quantile_list.rds")
Top_promoters_MM_other_types <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_promoter_ChIP-seq_top_5_list_other_cell_types.rds")
av_met_MM <- readRDS("/g/krebs/moene/for_manuscript/analysis/R_objects/MM/average_meth.rds")    # average methylation from Sonmezer et al., 2021
AllC_MM <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/data/MM/AllC.rds")          # position of all cytosines
baits.i <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/data/MM/baitsi.rds")        # position of baits

av_met_MM <- av_met_MM[baits.i, ]

### merge methylation data from C2C12 and MEL
# load the C2C12 and MEL data from Kreibich et al., 2023
C2C12_met_R1 <- readRDS("/g/krebs/krebs/analysis/SMF/MM/methCall/Context_methylation_call_SMF_MM_SMF_MM_C2C12_NO_R1.txt.rds")
C2C12_met_R1 <- C2C12_met_R1[baits.i, ]
C2C12_met_R2 <- readRDS("/g/krebs/krebs/analysis/SMF/MM/methCall/Context_methylation_call_SMF_MM_SMF_MM_C2C12_NO_R2.txt.rds")
C2C12_met_R2 <- C2C12_met_R2[baits.i, ]

MEL_met_R1 <- readRDS("/g/krebs/krebs/analysis/SMF/MM/methCall/Context_methylation_call_SMF_MM_SMF_MM_MEL_NO_R1.txt.rds")
MEL_met_R1 <- MEL_met_R1[baits.i, ]
MEL_met_R2 <- readRDS("/g/krebs/krebs/analysis/SMF/MM/methCall/Context_methylation_call_SMF_MM_SMF_MM_MEL_NO_R2.txt.rds")
MEL_met_R2 <- MEL_met_R2[baits.i, ]

# combine everything in one comparison matrix
av_met_MM <- cbind(av_met_MM, rowMeans(cbind(C2C12_met_R1$V1,C2C12_met_R2$V1)), rowMeans(cbind(MEL_met_R1$V1,MEL_met_R2$V1)))
colnames(av_met_MM) <- c("ES", "NP", "TKO", "C2C12", "MEL")
av_met_MM[av_met_MM == "NaN"] <- NA

### Add methylation data to AllC_MM
# select only for bait capture regions
AllC_MM <- AllC_MM[baits.i]
AllC_MM$V1 <- NULL
AllC_MM$ES_metNULLAllC_MM$ES_meth <- av_met_MM[,"ES"]
AllC_MM$NP_meth <- av_met_MM[,"NP"]
AllC_MM$TKO_meth <- av_met_MM[,"TKO"]
AllC_MM$C2C12_meth <- av_met_MM[,"C2C12"]
AllC_MM$MEL_meth <- av_met_MM[,"MEL"]

# the AllC_MM object is too large to upload
# subset this object by overlap with the promoter regions
TSS_extended <- resize(TSSsc_MM, fix = "center", width = 1000)
AllC_MM %>% filter_by_overlaps(TSS_extended) -> AllC_MM_subset
saveRDS(AllC_MM_subset, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_avg_meth_all_cell_types.rds")


########## Plots ##########
AllC_MM_subset <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_avg_meth_all_cell_types.rds")

# select only GC index (for other cell types than TKO)
AllC_MM_GC <- AllC_MM_subset[AllC_MM_subset$type == "DGCHN",]

### Function to subset cytosines from promoters of interest (e.g. top 5%) and plot profile of SMF signal
plot_composite_SMF <- function(TSSs, TSSsc_MM, AllC_MM, meth_call, width = 500, plot_title, save_file_name){
  
  # subset cytosines
  plotC <- as.matrix(findOverlaps(resize(TSSsc_MM, fix = "center", width = width)[TSSs], AllC_MM, ignore.strand = TRUE))
  
  # calculate distance between Cs and corresponding TSSsc
  Neg.i <- strand(TSSsc_MM[TSSs][plotC[,1]]) == "-"
  orientation <- rep(1, length(TSSsc_MM[TSSs][plotC[,1]]))
  orientation[as.vector(Neg.i)] <- -1
  
  dist_from_TSS <- (start(AllC_MM[plotC[,2]]) - start(TSSsc_MM[TSSs][plotC[,1]])) * orientation
  
  SMF_df <- data.frame(distance = dist_from_TSS, SMF = (1 - AllC_MM@elementMetadata[[meth_call]][plotC[,2]]) * 100)
  SMF_df <- SMF_df[complete.cases(SMF_df),]
  
  SMF_df_avg <- aggregate(.~distance, data=SMF_df, mean)
  
  # create a trend line
  avg_profile <- caTools::runmean(SMF_df_avg$SMF, 20, endrule = "constant")
  avg_profile <- data.frame(distance = SMF_df_avg$distance, SMF = avg_profile)
  
  # plot
  # pdf(paste0("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/analysis/SMF/MM_analysis/composite_plots/", save_file_name),
  #    height = 4.5, width = 4.5, useDingbats = FALSE)
  p <- ggplot(SMF_df, aes(x = distance, y = SMF)) +
    geom_point(size = 0.2, alpha = 0.2) + 
    geom_line(data = avg_profile, aes(x = distance, y = SMF), linewidth = 1, color = "red") +
    #stat_summary_bin(aes(y = SMF), fun = mean, colour="red", geom = "smooth", bins = 100) +
    xlab("Position relative to TSS") +
    ylab("SMF (%)") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.title = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(colour = "black"))
  print(p)
  dev.off()
  
}


# Top 5% mouse TSSs (TKO) *** Figure 1B ***
plot_composite_SMF(TSSs = Top_promoters_MM_PolII_ChIP$top5,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_subset,
                   meth_call = "TKO_meth",
                   width = 500,
                   plot_title = "Top 5% mouse TSSs",
                   save_file_name = "SMF_MM_TKO_top5_percent_TSSs.pdf")

# Top 5% TATA-containing mouse TSSs (TKO) *** Figure 1C ***
plot_composite_SMF(TSSs = Top_promoters_MM_PolII_ChIP$top5[Top_promoters_MM_PolII_ChIP$top5 %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id],
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_subset,
                   meth_call = "TKO_meth",
                   width = 500,
                   plot_title = "Top 5% mouse TSSs (TATA)",
                   save_file_name = "SMF_MM_TKO_top5_percent_TSSs_TATA.pdf")

# Bottom 10% mouse TSSs (TKO) *** Figure 1D ***
plot_composite_SMF(TSSs = Top_promoters_MM_PolII_ChIP$bottom10,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_subset,
                   meth_call = "TKO_meth",
                   width = 500,
                   plot_title = "Bottom 10% mouse TSSs",
                   save_file_name = "SMF_MM_TKO_bottom10_percent_TSSs.pdf")

### Other cell types (*** Figure S1 ***)
# Top 5% mouse TSSs (ES)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_mESC,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "ES_meth",
                   width = 500,
                   plot_title = "Embryonic stem cells",
                   save_file_name = "SMF_MM_ES_top5_percent_TSSs.pdf")

# Top 5% TATA-containing mouse TSSs (ES)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_mESC[Top_promoters_MM_other_types$top5_mESC %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id],
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "ES_meth",
                   width = 500,
                   plot_title = "Embryonic stem cells (TATA)",
                   save_file_name = "SMF_MM_ES_top5_percent_TSSs_TATA.pdf")

# Top 5% mouse TSSs (NP)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_NP,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "NP_meth",
                   width = 500,
                   plot_title = "Neural progenitor cells",
                   save_file_name = "SMF_MM_NP_top5_percent_TSSs.pdf")

# Top 5% TATA-containing mouse TSSs (NP)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_NP[Top_promoters_MM_other_types$top5_NP %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id],
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "NP_meth",
                   width = 500,
                   plot_title = "Neural progenitor cells (TATA)",
                   save_file_name = "SMF_MM_NP_top5_percent_TSSs_TATA.pdf")

# Top 5% mouse TSSs (C2C12)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_C2C12,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "C2C12_meth",
                   width = 500,
                   plot_title = "Myoblasts",
                   save_file_name = "SMF_MM_C2C12_top5_percent_TSSs.pdf")

# Top 5% TATA-containing mouse TSSs (C2C12)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_C2C12[Top_promoters_MM_other_types$top5_C2C12 %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id],
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "C2C12_meth",
                   width = 500,
                   plot_title = "Myoblasts (TATA)",
                   save_file_name = "SMF_MM_C2C12_top5_percent_TSSs_TATA.pdf")

# Top 5% mouse TSSs (MEL)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_MEL,
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "MEL_meth",
                   width = 500,
                   plot_title = "Murine erythroleukemia",
                   save_file_name = "SMF_MM_MEL_top5_percent_TSSs.pdf")

# Top 5% TATA-containing mouse TSSs (MEL)
plot_composite_SMF(TSSs = Top_promoters_MM_other_types$top5_MEL[Top_promoters_MM_other_types$top5_MEL %in% TATA_promoters_MM[TATA_promoters_MM$TATA == TRUE,]$gene_id],
                   TSSsc_MM = TSSsc_MM,
                   AllC_MM = AllC_MM_GC,
                   meth_call = "MEL_meth",
                   width = 500,
                   plot_title = "Murine erythroleukemia (TATA)",
                   save_file_name = "SMF_MM_MEL_top5_percent_TSSs_TATA.pdf")

