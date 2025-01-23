########## Visualising Pol II occupancy, turnover time, and initiation time ##########
# Author: Kasit Chatsirisupachai
# LastUpdate: 23.01.2025

library(tidyverse)
library(ggplot2)
library(ggExtra)
library(readr)

### Read data
MM <- read_csv("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/output/Kasit2024_table_parameters_mouse.csv")
DM <- read_csv("/g/krebs/chatsiri/mouse_droso_PolII/Chatsirisupachai_Moene_2024/data/Modelling/output/Kasit2024_table_parameters_drosophila.csv")

MM %>%
  mutate(species = "Mouse") -> MM

DM %>% 
  mutate(species = "Drosophila") -> DM

df <- rbind(MM, DM)
df

df %>%
  filter(is.finite(`initiation_time_[min]`)) %>%
  mutate(log2_initiation = log2(`initiation_time_[min]`)) %>%
  arrange(species) -> df

# median values
df %>% 
  group_by(species) %>% 
  summarise(median(log2_initiation), median(q)) -> median_values

### plot
p <- ggplot(df, aes(x = q, y = log2_initiation, color = species)) +
  geom_point(size = 2.2, alpha = 0.25, stroke = NA) +
  labs(x = "Pol II occupancy (q)", y = "log2(initiation time (min))") +
  xlim(0, 1) +
  scale_color_manual(values = c("Drosophila" = "darkblue", "Mouse" = "darkred")) +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p <- ggMarginal(p, 
                type = "density", 
                margins = "both",
                size = 5,
                groupColour = TRUE,
                groupFill = TRUE)


#pdf("/g/krebs/chatsiri/mouse_droso_PolII/re_analysis/analysis/Modelling/Visualise_PolII_occupancy_and_Initiation_time.pdf",
#    width = 4.5, height = 4.5, useDingbats = FALSE)
print(p)
dev.off()


