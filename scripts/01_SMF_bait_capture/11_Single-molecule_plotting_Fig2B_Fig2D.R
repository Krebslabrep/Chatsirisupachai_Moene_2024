########## Plot single-molecule examples for mouse and Drosophila promoters ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 29.08.2024

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(stringi)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

### This script is used to plot single molecule examples of mouse (Skp1; Fig 2D) and Drosophila (Svil; Fig 2B) promoters
### We used amplicon data to increase the number of molecules for plotting.

##### Mouse #####
### Qinput
# this file points toward BAM files of amplicon SMF in TKO mESCs from this study (E-MTAB-14461)
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_amplicon_SMF_MM_TRP_experiment.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
MySample <- MySample[1:2]   # since we are only interested in DMSO condition in this script

### amplicon regions
amplicons <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_amplicon_regions.rds")
amplicons <- amplicons[grep("^NM|^NR", amplicons$TFBSname)]

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/MM_REFSEQ_reference_transcripts_CAGE_corrected.rds")
TSSsc_of_interest <- TSSsc[unlist(TSSsc$gene_id) %in% unique(amplicons$TFBSname)]

### Load SM sorting merged freq
SM_all <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/amplicon_SMF/MM_TKO_DE_amplicon_TRP_sorting_output.rds")

TSS <- "NM_011543.4" # (Skp1)
Region_of_interest <- resize(TSSsc_of_interest[TSS], width = 1000, fix = "center")

# call context methylation
Methylation <- CallContextMethylation(sampleSheet = Qinput, 
                                      sample = MySample, 
                                      genome = BSgenome.Mmusculus.UCSC.mm10, 
                                      RegionOfInterest = Region_of_interest, 
                                      coverage = 20, 
                                      ConvRate.thr = 0.2)

### Randomly subset reads
# Set a seed for reproducibility
set.seed(123)

# Set sample size
sample_size <- 3000

# Randomly sample row indices
num_rows <- nrow(Methylation[[2]]$DMSO_DE_)
sampled_indices <- sample(1:num_rows, sample_size)

# Subset the sparse matrix
subset_matrix <- Methylation[[2]]$DMSO_DE_[sampled_indices, ]

# Replace Methylation object with the subset
Methylation[[2]]$DMSO_DE_ <- subset_matrix

### Plot in a relative location to TSS
## rearrange genomic location
TSS_coordinate <- start(TSSsc[TSS])
start(Methylation[[1]]@ranges) <- start(Methylation[[1]]@ranges) - TSS_coordinate
end(Methylation[[1]]@ranges) <- end(Methylation[[1]]@ranges) - TSS_coordinate

colnames(Methylation[[2]]$DMSO_DE_) <- start(Methylation[[1]]@ranges)

## region of interest
my_ROI <- GRanges(seqnames = TSS, ranges = IRanges(start = -200, end = 200), strand = "+")

# make bins
# use TSS location at 0 (because we want to plot in a relative space of the TSS)
bins <- MakeBins(RegionsOfInterest = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"), 
                 BinType = "Promoter", StrandAware = TRUE)
bins$TF <- bins$name

## Sort reads
# sort reads by promoter bins
SortedReads_Promoter <- SortReadsBySinglePromoter(MethSM = Methylation[[2]], 
                                                  TSS = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"),
                                                  coverage = 20)

SortedReads_Promoter$DMSO_DE_ %>% unlist() %>% length() -> num_reads 

bins_df <- data.frame(bins)

## Single-molecule stack plot
PlotSM(MethSM = Methylation[[2]], 
       RegionOfInterest = my_ROI, 
       SortedReads = SortedReads_Promoter, 
       sorting.strategy = "promoter") +
  xlab("position relative to TSS (bp)") +
  xlim(c(-200, 200)) + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks = element_blank()
  ) -> SM_pl

ggplot() + geom_rect(bins_df, 
                     mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 150), 
                     inherit.aes = FALSE) + 
  xlab("position relative to TSS (bp)") +
  xlim(c(-200, 200)) + 
  ylim(c(-100, 150)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) -> bin_pl

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


# add state
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

# colour set
colour_set <- colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,2,3,1,5,4)]#[c(2,3,4,9)]

# plot
PlottingDF %>%
  ggplot(aes(x=1, y=ReadID)) + 
  geom_tile(aes(fill=Promoter_State), height=1, width=0.5) +
  facet_wrap(~Sample, scales = "free_y", dir = 'v') +
  ylab("") +
  xlab("") +
  scale_discrete_manual(aesthetics = "fill", values = rev(colour_set)) +
  theme_classic() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
  ) -> StateQuant_pl


##### *** Figure 2D *** #####
## Combining plots
layout <- "
  AB
  C#
  "
SM_pl + StateQuant_pl + bin_pl +
  patchwork::plot_layout(ncol = 2, design = layout, widths = c(1, 0.05), heights = c(1, 0.1), guides = "collect") -> FinalPlot

#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/MM_analysis/composite_plots/SKP1_example_plot_with_bins.pdf", 
#    width = 6, height = 4.5, useDingbats = FALSE)
print(FinalPlot)
dev.off()


### State frequency
PlottingDF %>% 
  group_by(Sample) %>% 
  count(Promoter_State) %>% 
  mutate(frequency = prop.table(n) * 100) -> Promoter_State_Frequency




##### Drosophila #####
### Qinput
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_amplicon_SMF_DM_TRP_experiment.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
MySample <- MySample[1:2]   # since we are only interested in DMSO condition in this script

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Load SM sorting merged freq
merged_freq <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/amplicon_SMF/DM_S2_DE_amplicon_TRP_promoter_states_freq_matrix_collapsed.rds")

TSS <- "NM_001031933.3"   # Svil
Region_of_interest <- resize(TSSsc[TSS], width = 1000, fix = "center")

# call context methylation
Methylation <- CallContextMethylation(sampleSheet = Qinput, 
                                      sample = MySample, 
                                      genome = BSgenome.Dmelanogaster.UCSC.dm6, 
                                      RegionOfInterest = Region_of_interest, 
                                      coverage = 20, 
                                      ConvRate.thr = 0.2)

### Randomly subset reads
# Set a seed for reproducibility
set.seed(123)

# Set sample size
sample_size <- 3000

# Randomly sample row indices
num_rows <- nrow(Methylation[[2]]$DMSO_DE_S2)
sampled_indices <- sample(1:num_rows, sample_size)

# Subset the sparse matrix
subset_matrix <- Methylation[[2]]$DMSO_DE_S2[sampled_indices, ]

# Replace Methylation object with the subset
Methylation[[2]]$DMSO_DE_S2 <- subset_matrix

### Plot in a relative location to TSS
## rearrange genomic location
TSS_coordinate <- start(TSSsc[TSS])
start(Methylation[[1]]@ranges) <- start(Methylation[[1]]@ranges) - TSS_coordinate
end(Methylation[[1]]@ranges) <- end(Methylation[[1]]@ranges) - TSS_coordinate

colnames(Methylation[[2]]$DMSO_DE_S2) <- start(Methylation[[1]]@ranges)

## region of interest
my_ROI <- GRanges(seqnames = TSS, ranges = IRanges(start = -200, end = 200), strand = "+")

# make bins
# use TSS location at 0 (because we want to plot in a relative space of the TSS)
bins <- MakeBins(RegionsOfInterest = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"), 
                 BinType = "Promoter", StrandAware = TRUE)
bins$TF <- bins$name

## Sort reads
# sort reads by promoter bins
SortedReads_Promoter <- SortReadsBySinglePromoter(MethSM = Methylation[[2]], 
                                                  TSS = GRanges(seqnames = TSS, ranges = IRanges(start = 0), strand = "+"),
                                                  coverage = 20)

SortedReads_Promoter$DMSO_DE_S2 %>% unlist() %>% length() -> num_reads 

bins_df <- data.frame(bins)

## Single-molecule stack plot
PlotSM(MethSM = Methylation[[2]], 
       RegionOfInterest = my_ROI, 
       SortedReads = SortedReads_Promoter, 
       sorting.strategy = "promoter") +
  xlab("position relative to TSS (bp)") +
  xlim(c(-200, 200)) + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks = element_blank()
  ) -> SM_pl

ggplot() + geom_rect(bins_df, 
                     mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 150), 
                     inherit.aes = FALSE) + 
  xlab("position relative to TSS (bp)") +
  xlim(c(-200, 200)) + 
  ylim(c(-100, 150)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) -> bin_pl

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


# add state
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

# colour set
colour_set <- colorRampPalette(brewer.pal(9,"Set1"))(9)[c(9,2,3,1,5,4)]#[c(2,3,4,9)]

# plot
PlottingDF %>%
  ggplot(aes(x=1, y=ReadID)) + 
  geom_tile(aes(fill=Promoter_State), height=1, width=0.5) +
  facet_wrap(~Sample, scales = "free_y", dir = 'v') +
  ylab("") +
  xlab("") +
  scale_discrete_manual(aesthetics = "fill", values = rev(colour_set)) +
  theme_classic() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
  ) -> StateQuant_pl


##### *** Figure 2B *** #####
## Combining plots
layout <- "
  AB
  C#
  "
SM_pl + StateQuant_pl + bin_pl +
  patchwork::plot_layout(ncol = 2, design = layout, widths = c(1, 0.05), heights = c(1, 0.1), guides = "collect") -> FinalPlot


#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/SMF/DM_analysis/composite_plots/Svil_example_plot_with_bins.pdf", 
#    width = 6, height = 4.5, useDingbats = FALSE)
print(FinalPlot)
dev.off()


### State frequency
PlottingDF %>% 
  group_by(Sample) %>% 
  count(Promoter_State) %>% 
  mutate(frequency = prop.table(n) * 100) -> Promoter_State_Frequency
