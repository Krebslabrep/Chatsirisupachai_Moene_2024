########## Plot single-molecule of the TRP analysis (For representative TATA-less promoter NM_023665.3) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringi)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

### Qinput
# this file points toward BAM files of amplicon SMF in TKO mESCs from this study (E-MTAB-14461)
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_amplicon_SMF_MM_TRP_experiment.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])

### amplicon regions
amplicons <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_amplicon_regions.rds")
amplicons <- amplicons[grep("^NM|^NR", amplicons$TFBSname)]

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc_of_interest <- TSSsc[unlist(TSSsc$gene_id) %in% unique(amplicons$TFBSname)]

### Load SM sorting merged freq
merged_freq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_TKO_DE_amplicon_TRP_promoter_states_freq_matrix_collapsed.rds")

TSS <- "NM_023665.3"
Region_of_interest <- resize(TSSsc_of_interest[TSS], width = 1000, fix = "center")

# call context methylation
Methylation <- CallContextMethylation(sampleSheet = Qinput, 
                                      sample = MySample, 
                                      genome = BSgenome.Mmusculus.UCSC.mm10, 
                                      RegionOfInterest = Region_of_interest, 
                                      coverage = 20, 
                                      ConvRate.thr = 0.2)

Methylation[[1]]

### Randomly subset reads
# Set a seed for reproducibility
set.seed(123)

# Set sample size
sample_size <- 3000

### Randomly sample row indices
# DMSO
num_rows <- nrow(Methylation[[2]]$DMSO_DE_)
sampled_indices <- sample(1:num_rows, sample_size)

# Subset the sparse matrix
subset_matrix <- Methylation[[2]]$DMSO_DE_[sampled_indices, ]

# Replace Methylation object with the subset
Methylation[[2]]$DMSO_DE_ <- subset_matrix

# TRP
num_rows <- nrow(Methylation[[2]]$TRP_DE_30min)
sampled_indices <- sample(1:num_rows, sample_size)

# Subset the sparse matrix
subset_matrix <- Methylation[[2]]$TRP_DE_30min[sampled_indices, ]

# Replace Methylation object with the subset
Methylation[[2]]$TRP_DE_30min <- subset_matrix

### Plot in a relative location to TSS
## rearrange genomic location
TSS_coordinate <- start(TSSsc[TSS])
TSS_strand <- as.character(strand(TSSsc[TSS]))

if(TSS_strand == "+"){
  start(Methylation[[1]]@ranges) <- (start(Methylation[[1]]@ranges) - TSS_coordinate)
  end(Methylation[[1]]@ranges) <- (end(Methylation[[1]]@ranges) - TSS_coordinate)
} else if (TSS_strand == "-"){
  start(Methylation[[1]]@ranges) <- -(start(Methylation[[1]]@ranges) - TSS_coordinate)
  end(Methylation[[1]]@ranges) <- -(end(Methylation[[1]]@ranges) - TSS_coordinate)
  
}

colnames(Methylation[[2]]$DMSO_DE_) <- start(Methylation[[1]]@ranges)
colnames(Methylation[[2]]$TRP_DE_30min) <- start(Methylation[[1]]@ranges)

## region of interest
my_ROI <- GRanges(seqnames = TSS, ranges = IRanges(start = -200, end = 200), strand = "+")

# make bins
# use TSS location at 0 (because we want to plot in a relative space of the TSS)
bins <- MakeBins(RegionsOfInterest = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"), 
                 BinType = "Promoter", StrandAware = TRUE)
bins$TF <- bins$name

## Average SMF plot
PlotAvgSMF(MethGR = Methylation[[1]],
           MethSM = Methylation[[2]],
           RegionOfInterest = my_ROI,
           SortedReads = NULL,
           ShowContext = FALSE,
           TFBSs = bins,
           SortingBins = NULL) -> Avg_pl
Avg_pl + geom_vline(xintercept = 0,
                    linetype = "dashed", 
                    color = "#737373") +
  ggtitle("Rsrp1 (chr1:134923425-134923824)") -> Avg_pl


## Sort reads
# sort reads by promoter bins
SortedReads_Promoter <- SortReadsBySinglePromoter(MethSM = Methylation[[2]], 
                                                  TSS = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"),
                                                  coverage = 20)

## Single-molecule stack plot
PlotSM(MethSM = Methylation[[2]], 
       RegionOfInterest = my_ROI, 
       SortedReads = SortedReads_Promoter, 
       sorting.strategy = "promoter") -> SM_pl


## State quantification plot
Strand <- as.character(strand(my_ROI))
states <- Promoterstates()
OrderedReads <- lapply(SortedReads_Promoter, function(sR){sR[as.character(unlist(states))]})

Reduce(rbind,
       lapply(seq_along(OrderedReads), function(i){
         Reduce(rbind,
                lapply(seq_along(OrderedReads[[i]]), function(j){
                  if(!is.null(OrderedReads[[i]][[j]])){
                    tibble(ReadID = OrderedReads[[i]][[j]], 
                           Pattern = names(OrderedReads[[i]])[j], 
                           Sample = names(OrderedReads)[i])
                  }
                }))
       })) -> OrderedReads_tbl


full_join(OrderedReads_tbl, rownames_to_column(data.frame(Pattern = unlist(states)), "State"), "Pattern") %>% 
  na.omit() %>% 
  separate(Pattern, into = c(paste0("Bin", seq(unique(nchar(unlist(states)))))), sep = "(?<=.)", extra = 'drop') %>%
  gather(Bin, Methylation, -ReadID, -Sample, -State) -> PlottingDF
PlottingDF$ReadID <- factor(PlottingDF$ReadID, levels = unlist(OrderedReads))


# add new column for colour purpose
colour_set <- colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,2,3,1,5,4)]#[c(2,3,4,9)]
colour_set <- c("#FFFFFF", colour_set)


##### Simpler state quantification plot
simple_prom_state <- function(i, PlottingDF){
  tmp <- PlottingDF[i,]
  
  if(grepl("unassigned", tmp$State)){return("unassigned")}
  else if(grepl("nucleosome", tmp$State)){return("nucleosome")}
  else if(grepl("unbound", tmp$State)){return("unbound")}
  else if(tmp$State == "PIC"){return("PIC")}
  else if(tmp$State == "PIC.polII"){return("PIC.polII")}
  else if(tmp$State == "polII"){return("polII")}
}

PlottingDF$Promoter_State <- unlist(lapply(seq_along(1:nrow(PlottingDF)), simple_prom_state, PlottingDF = PlottingDF))
PlottingDF$Promoter_State <- factor(PlottingDF$Promoter_State, levels = c("polII", "PIC.polII", "PIC", "unbound", "nucleosome", "unassigned"))

PlottingDF %>%
  ggplot(aes(x=1, y=ReadID)) + 
  geom_tile(aes(fill=Promoter_State), height=1, width=0.5) +
  facet_wrap(~Sample, scales = "free_y", dir = 'v') +
  ylab("") +
  xlab("") +
  scale_discrete_manual(aesthetics = "fill", values = rev(colour_set[2:7])) +
  theme_classic() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank()
        #strip.background = element_blank(),
        #strip.text.x = element_blank()
  ) -> StateQuant_pl


##### *** Figure S4 *** #####
layout <- "
  #A
  CB
  "
Avg_pl + SM_pl + StateQuant_pl +
  patchwork::plot_layout(ncol = 2, design = layout, widths = c(0.05, 1), heights = c(1, 0.8), guides = "collect") -> FinalPlot


#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/TRP_analysis/Rsrp1_example_simpler_quantification_plot.pdf", 
#    width = 8, height = 8, useDingbats = FALSE)
print(FinalPlot)
dev.off()


### State frequency
PlottingDF %>% 
  group_by(Sample) %>% 
  count(Promoter_State) %>% 
  mutate(frequency = prop.table(n) * 100) -> Promoter_State_Frequency



