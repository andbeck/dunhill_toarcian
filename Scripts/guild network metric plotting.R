# plotting network stats through time
# Toarcian
# Dunhill et al 2022

# libraries ----
library(tidyverse)
library(patchwork)

# data and data subset to guild level  ----
## this is andrew based data using GuildWebsBuildPlotAnalyse and calc_select_stats
guilds <- read_csv("Data/metrics_time2.csv")

## these are alex's original data
## did not use calc_select_stats
## this uses calc_network_stats 

guilds2 <- read_csv("Data/metrics_time.csv") %>% filter(resolution == "guild") %>% 
  select(-resolution)

colnames(guilds)
colnames(guilds2)

guilds %>% select(sd_normalized_in_degree, sd_normalized_out_degree)
guilds2 %>% select(generality, vulnerability)


# data management to create axis labels ----
guilds <- guilds %>% 
  mutate(time2 = case_when(
    time == 1 ~ "pre-extinction",
    time == 2 ~ "post-extinction",
    time == 3 ~ "early recovery",
    time == 4 ~ "late recovery")) %>% 
  mutate(time2 = factor(time2, levels = c("pre-extinction",
                                          "post-extinction",
                                          "early recovery",
                                          "late recovery")))


# subsets of Structure metrics and Motif metrics ----
guildsStructureDat <- guilds %>% 
  select(time, time2, taxa, connectance, max_tl_std, sd_normalized_in_degree, sd_normalized_out_degree) %>% 
  pivot_longer(-c(time, time2), names_to = "Metric", values_to = "Value") %>% 
  mutate(Metric = case_when(
    Metric == "taxa" ~ "a-Guild Richness",
    Metric == "connectance" ~ "b-Connectance",
    Metric == "max_tl_std" ~ "c-Max Trophic Chain Length",
    Metric == "sd_normalized_in_degree" ~ "d-Generality",
    Metric == 'sd_normalized_out_degree' ~ "e-Vulnerability"
  ))

guildsMotifDat <- guilds %>% 
  select(time, time2, norm_mot_lin, norm_mot_omn, norm_mot_ap_comp, norm_mot_dir_comp) %>% 
  pivot_longer(-c(time, time2), names_to = "Metric", values_to = "Value") %>% 
  mutate(Metric = case_when(
    Metric == "norm_mot_lin" ~ "Linear Food Chain",
    Metric == "norm_mot_omn" ~ "Omnivory",
    Metric == "norm_mot_ap_comp" ~ "Apparent Competition",
    Metric == "norm_mot_dir_comp" ~ "Competition"
  ))

guildsMotifDat2 <- guilds2 %>% 
  select(time, s1, s2, s4, s5) %>% 
  pivot_longer(-c(time), names_to = "Metric", values_to = "Value") %>% 
  mutate(Metric = case_when(
    Metric == "s1" ~ "Linear Food Chain",
    Metric == "s2" ~ "Omnivory",
    Metric == "s4" ~ "Apparent Competition",
    Metric == "s5" ~ "Competition"
  )) %>% 
  mutate(time2 = case_when(
    time == 1 ~ "pre-extinction",
    time == 2 ~ "post-extinction",
    time == 3 ~ "early recovery",
    time == 4 ~ "late recovery")) %>% 
  mutate(time2 = factor(time2, levels = c("pre-extinction",
                                          "post-extinction",
                                          "early recovery",
                                          "late recovery")))


# plots ----
structPlot <- ggplot(guildsStructureDat, aes(x = time2, y =Value, group = Metric))+
  geom_line()+
  facet_wrap(~ Metric, scales = "free_y", ncol = 1)+
  labs(x = NULL, y = NULL, title = "(c)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# structPlot

# motifPlot <- ggplot(guildsMotifDat, aes(x = time2, y =Value, group = Metric))+
#   geom_line()+
#   facet_wrap(~ Metric, scales = "free_y", ncol = 1)+
#   labs(x = NULL, y = NULL)+
#   theme_bw(base_size = 15)+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
# 
# motifPlot

motifPlot2 <- ggplot(guildsMotifDat2, aes(x = time2, y =Value, group = Metric))+
  geom_line()+
  facet_wrap(~ Metric, scales = "free_y", ncol = 1)+
  labs(x = NULL, y = NULL, title = "(d)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# motifPlot2


# patchwork layout ----
# structPlot+motifPlot
