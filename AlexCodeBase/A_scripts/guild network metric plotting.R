# plotting network stats through time
# Toarcian
# Dunhill et al 2022

# libraries ----
library(tidyverse)
library(patchwork)

# data and data subset to guild level  ----
stats <- read_csv("AlexCodeBase/A_data/metrics_time.csv")
guilds <- stats[grep(TRUE,stats[,"resolution"] == "guild"),]

# data management to create axis labels ----
guilds <- guilds %>% 
  mutate(time2 = case_when(
    time == 1 ~ "pre-extinction",
    time == 2 ~ "post-extinction",
    time == 3 ~ "early recovery",
    time == 4 ~ "late recovery"))

# subsets of Structure metrics and Motif metrics ----
guildsStructureDat <- guilds %>% 
  select(time, time2, taxa, connectance, max_tl, generality, vulnerability) %>% 
  pivot_longer(-c(time, time2), names_to = "Metric", values_to = "Value")

guildsMotifDat <- guilds %>% 
  select(time, time2, s1, s2, s4, s5) %>% 
  pivot_longer(-c(time, time2), names_to = "Metric", values_to = "Value")

# plots ----
structPlot <- ggplot(guildsStructureDat, aes(x = time2, y =Value, group = Metric))+
  geom_line()+
  facet_wrap(~ Metric, scales = "free_y", ncol = 1)+
  labs(x = NULL, y = NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

structPlot

motifPlot <- ggplot(guildsMotifDat, aes(x = time2, y =Value, group = Metric))+
  geom_line()+
  facet_wrap(~ Metric, scales = "free_y", ncol = 1)+
  labs(x = NULL, y = NULL)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

motifPlot


# patchwork layout ----
#structPlot+motifPlot
