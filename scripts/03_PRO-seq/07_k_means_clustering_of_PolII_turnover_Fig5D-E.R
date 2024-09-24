########## k-means clustering of the promoters by Pol II turnover rates ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 24.09.2024

library(tidyverse)
library(ggplot2)
library(pheatmap)

### Load data
### Top 5
mouse_top5_FC <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/mouse_TKO_TRP_promoter_counts_SpikeIn_norm_avg_top5_FC.rds")
droso_top5_FC <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm_avg_top5_FC.rds")

droso_top5_FC %>% drop_na() -> droso_top5_FC_1
mouse_top5_FC %>% drop_na() -> mouse_top5_FC_1

# mouse
MM_obj <- pheatmap(mouse_top5_FC_1, 
                   kmeans_k = 4, 
                   cluster_cols = F, 
                   display_numbers = T)
MM_obj$kmeans
MM_obj$kmeans$cluster %>% as.data.frame() %>% rename(cluster = names(.)[1]) -> MM_cluster

saveRDS(MM_obj, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/MM_K_means_clustering.rds")
saveRDS(MM_cluster, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/MM_K_means_clustering_TSS_cluster_mapping.rds")

# droso
DM_obj <- pheatmap(droso_top5_FC_1, 
                   kmeans_k = 4, 
                   cluster_cols = F, 
                   display_numbers = T)
DM_obj$kmeans
DM_obj$kmeans$cluster %>% as.data.frame() %>% rename(cluster = names(.)[1]) -> DM_cluster

saveRDS(DM_obj, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/DM_K_means_clustering.rds")
saveRDS(DM_cluster, "/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/DM_K_means_clustering_TSS_cluster_mapping.rds")


##### Heatmap with clusters #####
mouse_obj <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/MM/PRO_seq/MM_K_means_clustering.rds")
droso_obj <- readRDS("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/DM/PRO_seq/DM_K_means_clustering.rds")

### Mouse
mouse_obj$kmeans$cluster[mouse_obj$kmeans$cluster == 1] %>% names() -> cluster_1
mouse_obj$kmeans$cluster[mouse_obj$kmeans$cluster == 2] %>% names() -> cluster_2
mouse_obj$kmeans$cluster[mouse_obj$kmeans$cluster == 3] %>% names() -> cluster_3
mouse_obj$kmeans$cluster[mouse_obj$kmeans$cluster == 4] %>% names() -> cluster_4

# sort data by cluster
TSS_order <- c(cluster_3, cluster_4, cluster_1, cluster_2)
mouse_top5_FC_2 <- mouse_top5_FC_1[TSS_order,]

# use 1 as a maximum value
mouse_top5_FC_2[mouse_top5_FC_2 > 1] <- 1

# change colnames
colnames(mouse_top5_FC_2) <- c("0 min", "2.5 min", "5 min", "10 min", "20 min")

# define the annotation
annotation_row = data.frame(
  Cluster = factor(rep(c("Cluster3", "Cluster4", "Cluster1", "Cluster2"), 
                       c(length(cluster_3), length(cluster_4), length(cluster_1), length(cluster_2))))
)

row.names(annotation_row) <- row.names(mouse_top5_FC_2)

cluster_colours <- list(Cluster = c(Cluster1 = "#fecc5c", Cluster2 = "#ffffb2", Cluster3 = "#e31a1c", Cluster4 = "#fd8d3c"))

##### *** Figure 5E *** #####
p <- pheatmap(mouse_top5_FC_2, 
              cluster_rows = FALSE, 
              cluster_cols = FALSE, 
              show_rownames = FALSE, 
              color = colorRampPalette(c("white", "#a50f15"))(50),
              annotation_row = annotation_row,
              annotation_colors = cluster_colours,
              angle_col = 90,
              annotation_legend = FALSE)
print(p)
dev.off()


### Drosophila
droso_obj$kmeans$cluster[droso_obj$kmeans$cluster == 1] %>% names() -> cluster_1
droso_obj$kmeans$cluster[droso_obj$kmeans$cluster == 2] %>% names() -> cluster_2
droso_obj$kmeans$cluster[droso_obj$kmeans$cluster == 3] %>% names() -> cluster_3
droso_obj$kmeans$cluster[droso_obj$kmeans$cluster == 4] %>% names() -> cluster_4

# sort data by cluster
TSS_order <- c(cluster_1, cluster_2, cluster_3, cluster_4)
droso_top5_FC_2 <- droso_top5_FC_1[TSS_order,]

# use 1 as a maximum value
droso_top5_FC_2[droso_top5_FC_2 > 1] <- 1

# define the annotation
annotation_row = data.frame(
  Cluster = factor(rep(c("Cluster1", "Cluster2", "Cluster3", "Cluster4"), 
                       c(length(cluster_1), length(cluster_2), length(cluster_3), length(cluster_4))))
)

row.names(annotation_row) <- row.names(droso_top5_FC_2)

cluster_colours <- list(Cluster = c(Cluster1 = "#225ea8", Cluster2 = "#41b6c4", Cluster3 = "#a1dab4", Cluster4 = "#ffffcc"))

##### *** Figure 5D *** #####
p <- pheatmap(droso_top5_FC_2, 
              cluster_rows = FALSE, 
              cluster_cols = FALSE, 
              show_rownames = FALSE, 
              color = colorRampPalette(c("white", "#225ea8"))(50),
              annotation_row = annotation_row,
              annotation_colors = cluster_colours)
print(p)
dev.off()


