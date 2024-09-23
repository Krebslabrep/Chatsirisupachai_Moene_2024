########## Get state frequencies out of SM sorting output (Drosophila) ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.09.2024

### Qinput
# this file points toward BAM files of whole-genome SMF in S2 cells from Krebs et al., 2017
Qinput <- "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Qinput_files/Qinput_whole-genome_SMF_DM.txt"
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])

### TSSs
TSSsc <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/DM_REFSEQ_reference_transcripts_CAGE_corrected.rds")

### Load SM sorting output
SM_all <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_sorting_output.rds")

### get SM sorting output as a data frame
SM_df <- data.frame(SM_all[[3]])
SM_df <- SM_df[complete.cases(SM_df),]

### 16 possible states
allPos <- expand.grid(c(0,1), c(0,1), c(0,1), c(0,1))
states <- names(table(apply(allPos, 1, function(x){(paste(as.character((x)), collapse = ""))})))


##### Main function to create a state count matrix of all promoters for each sample
state_matrix_generator <- function(SM_df, TSSsc, states, count_or_freq){
  
  ### Get samples
  samples <- unique(SM_df$Sample)
  samples <- samples[!is.na(samples)]
  
  ### working on each sample
  stateM_list <- lapply(samples, function(sample){
    
    print(paste0("Creating state matrix for: ", sample))
    
    # subset df to select only a sample we are working with
    SM_df_sample <- SM_df[SM_df$Sample == sample,]
    
    # create matrix to store states
    stateM <- matrix(nrow = length(TSSsc), ncol = length(states))
    rownames(stateM) <- names(TSSsc)
    colnames(stateM) <- states
    
    # get TSSs
    sample_TSSs <- unique(SM_df_sample$TSS)
    
    # select counts or frequencies
    if(count_or_freq == "count"){
      # loop through all TSSs and add counts to the matrix
      for(TSS in sample_TSSs){
        tmp_df <- SM_df_sample[SM_df_sample$TSS == TSS, ]
        for(state in tmp_df$State){
          stateM[TSS, state] <- tmp_df[tmp_df$State == state, ]$Counts
        }
        stateM[is.na(stateM)] <- 0
      }
      return(stateM)
    } else if(count_or_freq == "freq"){
      # loop through all TSSs and add counts to the matrix
      for(TSS in sample_TSSs){
        tmp_df <- SM_df_sample[SM_df_sample$TSS == TSS, ]
        for(state in tmp_df$State){
          stateM[TSS, state] <- tmp_df[tmp_df$State == state, ]$Freqs
        }
        stateM[is.na(stateM)] <- 0
      }
      return(stateM)
    } else {
      stop("Please provide either 'count' or 'freq' ")
    }
  })
  
  names(stateM_list) <- samples
  
  return(stateM_list)
}

count_matrix_all_samples <- state_matrix_generator(SM_df, TSSsc, states, count_or_freq = "count")
freq_matrix_all_samples <- state_matrix_generator(SM_df, TSSsc, states, count_or_freq = "freq")


##### Collapse states into promoter states
### list promoter states
promoter_states <- Promoterstates()
promoter_statesF <- as.factor(unlist(lapply(seq_along(promoter_states), function(i){
  rep(names(promoter_states[i]), length(promoter_states[[i]]))
}))[order(unlist(promoter_states))])
promoter_statesF <- factor(promoter_statesF, levels = names(promoter_states))

### Function to collapse promoter states for each sample
collapse_states <- function(state_matrix, promoter_statesF){
  return(t(apply(state_matrix, 1, function(x){tapply(x, promoter_statesF, sum)})))
}

### Promoter state count matrix
grouped_count_mats <- lapply(count_matrix_all_samples, collapse_states, promoter_statesF = promoter_statesF)

### Promoter state frequency matrix
grouped_freq_mats <- lapply(freq_matrix_all_samples, collapse_states, promoter_statesF = promoter_statesF)

# Identify promoters that have been detected and sorted by SMF in each replicate
whole_genome_R1 <- grouped_freq_mats$SMF_DM_DE_S2_R1[round(rowSums(grouped_freq_mats$SMF_DM_DE_S2_R1)) == 100,]
dim(whole_genome_R1)     # 6021

whole_genome_R2 <- grouped_freq_mats$SMF_DM_DE_S2_R2[round(rowSums(grouped_freq_mats$SMF_DM_DE_S2_R2)) == 100,]
dim(whole_genome_R2)     # 6003

# select only promoters that are shared by both replicates
genes <- rownames(whole_genome_R1)[rownames(whole_genome_R1) %in% rownames(whole_genome_R2)]    # 5912 common genes between replicates
whole_genome_R1 <- whole_genome_R1[rownames(whole_genome_R1) %in% genes,]
whole_genome_R2 <- whole_genome_R2[rownames(whole_genome_R2) %in% genes,]

whole_genome_R1 <- whole_genome_R1[order(row.names(whole_genome_R1)), ]
whole_genome_R2 <- whole_genome_R2[order(row.names(whole_genome_R2)), ]


grouped_freq_mats <- list(SMF_DM_DE_S2_R1 = whole_genome_R1, SMF_DM_DE_S2_R2 = whole_genome_R2)

saveRDS(grouped_freq_mats, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix.rds")



# make average
whole_genome_avg <- (whole_genome_R1 + whole_genome_R2)/2

saveRDS(as.data.frame(whole_genome_avg), "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/genome_wide_SMF/DM_S2_DE_promoter_states_freq_matrix_avg.rds")

