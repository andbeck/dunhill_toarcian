### This code takes the outputs from GenerateExtinctionSequences 
# and builds graphs for manuscript
# YOU MUST RUN GenerateExtinctionSequences in order 
# to create the wrkWebs_allSeqs.RData file loaded here.

# we use functions from Shaw et al and cheddar and Anubhav et al to estimate 
## structural network properties
## motifs
## True Skill Statistic Metrics.

# CRAN R libraries ----
## libraries for network things ----
library(cheddar)
library(igraph)
library(NetIndices)
## libraries for graph things ----
library(tidyverse)
library(ggrepel)
library(patchwork)

# raw Guild data and other metric functions ----
# extra data on guilds and functions for calculating
# network metrics and motifs (see ReadMe for more information)
load("Data/ToarcianWebs_Guild_May2021.RData")
source("Scripts/NewMethod_Functions_update4Publication.R")
source("Scripts/pfim_scripts.R")

# LOAD wrkWebs here to start anaysis of structure/motifs and TSS ----
# These data were created by GenerateExtinctionSequences
# This must have been run first.
# These data are provided... you do not need to run them.

load("Data/wrkWebs_allSeqs.RData")

# check to see that data exists
ifelse(exists("wrkWebs_allSeqs", 1), 
       "PROCEED - Your Data exist", 
       "STOP, GenerateExtinctionSequences.R needed")

# setup masterSpecies list for TSS analysis---- 

# master species list contains all species in the pre and post webs
# the extWeb is then the set of 21 trophic links from the sims
# the testWeb is then the set of 21 trophic links from the raw G2 data
# we compare presence/absence of links in sims to raw, with core set of absences
# defined by masterSpeces 0's in raw data.
G1 <- read_csv("Data/G1_Guilds.csv")
G2 <- read_csv("Data/G2_Guilds.csv")

masterSpecies <- bind_rows(G1, G2) %>% 
  unique() %>% 
  rename("node" = "Guild") %>% 
  select(node, clade) %>% unique()

# Create collection data frame for metrics by sequence ----
masterStats <- data.frame(matrix(ncol = 4, nrow = 0))
names(masterStats) <- c("metric", "stat", "estimate", "trait")

# Generate and collect all metrics and motif data into single data frame ----
# applying the jack motif, jack nets and TSS
# see pfim scripts

for(i in 1:length(wrkWebs_allSeqs)){
  code = names(wrkWebs_allSeqs[i])
  
  # apply Jacks Motif Stats to the webs by each trait sequence
  motifStats <- sapply(wrkWebs_allSeqs[[i]], function(x) jackMotifs(x)) %>% t() %>% 
    data.frame() %>% 
    pivot_longer(everything(), names_to = "metric", values_to = "estimate") %>% 
    group_by(metric) %>% 
    summarise(
      meanVal = mean(estimate),
      sdVal = sd(estimate)
    )
  
  # apply Jack Net Stats selection to the webs by each trait sequence
  netStats_Jack <- sapply(wrkWebs_allSeqs[[i]], function(x) jackNetworks(x)) %>% t() %>%
    data.frame() %>% 
    pivot_longer(everything(), names_to = "metric", values_to = "estimate") %>% 
    group_by(metric) %>% 
    summarise(
      meanVal = mean(estimate),
      sdVal = sd(estimate)
    )
  
  # apply TSS stat from Anubhav to the webs by trait sequence
  tssStats <- sapply(wrkWebs_allSeqs[[i]], function(x) TSS_func_struc(x)) %>% 
    as_tibble(.) %>% 
    rename("TSS" = value) %>% 
    summarise(
      meanVal = mean(TSS),
      sdVal = sd(TSS)
    ) %>% 
    mutate(metric = "TSS") %>% 
    select(metric, meanVal, sdVal)
  
  outFull <- bind_rows(motifStats, netStats_Jack, tssStats) %>% 
    mutate(trait = paste(code))
  masterStats <- rbind(masterStats, outFull)
}

# define trait order for plotting on axes ----
traitOrder <- c("rand",
                "size_b2s","size_s2b",
                "tier_i2p","tier_p2i",
                "mot_fn","mot_nf",
                "calc_h2l","calc_l2h",
                "gen_l2h","gen_h2l",
                "vuln_l2h","vuln_h2l")

# define metric order for plotting ----
metricOrder <- c("connectance","mean_normalized_degree",
                 "sd_normalized_in_degree","sd_normalized_out_degree",
                 "mean_tl_std","max_tl_std","soi_std",
                 "size","diameter","mean_norm_btw",
                 "norm_mot_lin","norm_mot_omn","norm_mot_ap_comp","norm_mot_dir_comp",
                 "TSS")

# Plots of Networks (not used)  ----

# # 4 samples of the 21 spp final networks for each sequence
# 
# # reproducibility for subsampling networks visualisations
# set.seed(123)
# 
# # create network vis's (not used in MS)
# w <-  4
# par(mfrow = c(length(wrkWebs_allSeqs),w+1), mar = c(1,1,1,1))
# for(i in c(1,4,5,6,7,2,3,12,13,10,11,8,9)){
#   webs <- wrkWebs_allSeqs[[i]]
#   plot(1,1, axes = FALSE, type = "n")
#   text(1,1,label = names(wrkWebs_allSeqs)[i])
#   for (j in sample(1:length(webs), w , replace = FALSE)){
#     PlotWebByLevel(webs[[j]], main = "", frame.plot = TRUE)
#   }
# }

# Plots in Manuscript ----

## master data ----
mm <- masterStats %>% 
  mutate(trait = factor(trait, levels = traitOrder)) %>% 
  mutate(metric = factor(metric, levels = metricOrder)) %>% 
  filter(metric != "size")

## master reference data ----
postMetrics <- jackNetworks(postCom_Guild) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "metric") %>%
  rename("estimate" = ".") %>% 
  bind_rows(data.frame(metric = c("TSS"), estimate = c(NA))) %>% 
  mutate(metric = factor(metric, levels = metricOrder)) %>% 
  filter(metric != "size")

## subsets for plotting ----

# network metrics
nets <- mm %>% filter(metric != "TSS") %>% 
  filter(!grepl("mot", metric))
# motif metrics
mots <- mm %>% filter(grepl("mot", metric))
# True Skill metric
TSSs <- mm %>% filter(metric == "TSS") %>% 
  mutate(trait2 = case_when(
    trait == "rand" ~ "Random",
    trait == "size_b2s" ~ "Size (Large-Small)",
    trait == "size_s2b" ~ "Size (Small-Big)",
    trait == "tier_i2p" ~ "Tiering (infaunal-pelagic)",
    trait == "tier_p2i" ~ "Tiering (pelagic-infaunal)",
    trait == "mot_fn" ~ "Motility (fast-none)",
    trait == "mot_nf" ~ "Motility (none-fast)",
    trait == "calc_l2h" ~ "Calcified (low-high)",
    trait == "calc_h2l" ~ "Calcified (high-low)",
    trait == "gen_l2h" ~ "Generalism (low-high)",
    trait == "gen_h2l" ~ "Generalism (high-low)",
    trait == "vuln_l2h" ~ "Vulnerability (low-high)",
    trait == "vuln_h2l" ~ "Vulnerability (high-low)"
  ))

# post (empirical) motif and network metrics
post_nets <- postMetrics %>%
  filter(metric %in% c("connectance", "max_tl_std",
                       "sd_normalized_in_degree",
                       "sd_normalized_out_degree")) %>% 
  droplevels()

post_mots <- postMetrics %>% filter(grepl("mot", metric))

# use these
netsUse <- nets %>% 
  filter(metric %in% c("connectance", "max_tl_std",
                       "sd_normalized_in_degree",
                       "sd_normalized_out_degree")) %>% 
  droplevels()

## labels for network and motif metrics ----

net_labels <-  c(connectance = "Connectance",
                 sd_normalized_in_degree = "SD In Degree",
                 sd_normalized_out_degree = "SD Out Degree",
                 max_tl_std = "Max Trophic Level")

mot_labels <-  c(norm_mot_lin = "Linear Food Chain",
                 norm_mot_omn = "Omnivory",
                 norm_mot_ap_comp = "Apparent Competition",
                 norm_mot_dir_comp = "Direct Competition")


## Graph Format 1 (not used): Barplots of raw metrics ----

# net_plots <- ggplot(netsUse, aes(x = trait, y = meanVal, 
#                                  fill = trait))+
#   geom_col()+
#   ylab("Metric Score")+
#   geom_errorbar(aes(ymin = meanVal - sdVal,
#                     ymax = meanVal + sdVal),
#                 width = 0, size = 0.25)+
#   facet_wrap(~metric, scales = "free_y", ncol = 2, 
#              labeller = as_labeller(net_labels))+
#   geom_hline(data = post_nets, aes(yintercept = estimate))+
#   theme(axis.text.x = element_text(angle = 90),
#         strip.text.x = element_text(size = 5))+
#   guides(fill = "none")
# 
# mot_plots <- ggplot(mots,aes(x = trait, y = meanVal, fill = trait))+
#   geom_col()+
#   ylab("Metric Score")+
#   geom_errorbar(aes(ymin = meanVal - sdVal,
#                     ymax = meanVal + sdVal),
#                 width = 0, size = 0.25)+
#   facet_wrap(~metric, scales = "free_y", ncol = 2, 
#              labeller = as_labeller(mot_labels))+
#   geom_hline(data = post_mots,
#              aes(yintercept = estimate))+
#   theme(axis.text.x = element_text(angle = 90),
#         strip.text.x = element_text(size = 5))+
#   guides(fill = "none")
# 
# TSS Only
TSS_graph <- ggplot(TSSs,
                    aes(x = fct_reorder(trait2, meanVal), y = meanVal, fill = trait2))+
  geom_col()+
  ylab("TSS Score")+xlab("Sequence")+
  geom_errorbar(aes(ymin = meanVal - sdVal,
                    ymax = meanVal + sdVal),
                width = 0, linewidth = 0.25)+
  labs(x = NULL)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(size = 5))+
  guides(fill = "none")+
  coord_flip()

TSS_graph

# ### barplot version 1 ----
# TSS_graph/(net_plots+mot_plots)

## Data for Graph Format 2 and 3: difference to reference ----

## generate differences from empirical ----

# metrics
mets <- rep(metricOrder[metricOrder!="size"], each = 13)

# reference estimate
refEst <- rep(postMetrics$estimate, each = 13) %>% 
  tibble(.,mets) %>% 
  rename("refEst" = ".")

# difference calculation
diffDat <- bind_cols(arrange(mm, metric), refEst) %>% 
  mutate(diffEst = meanVal-refEst) %>% 
  filter(metric !="TSS")

# 20% of max difference to draw boundary lines
diffBounds <- diffDat %>% 
  group_by(metric) %>% 
  summarise(
    interval = 0.2*(max(abs(diffEst)))
  )

# relabelling and filtering difference data for net stats
diff_nets <- diffDat %>% filter(metric != "TSS") %>% 
  filter(!grepl("mot", metric)) %>% 
  mutate(trait2 = case_when(
    trait == "rand" ~ "Random",
    trait == "size_b2s" ~ "Size (Large-Small)",
    trait == "size_s2b" ~ "Size (Small-Big)",
    trait == "tier_i2p" ~ "Tiering (infaunal-pelagic)",
    trait == "tier_p2i" ~ "Tiering (pelagic-infaunal)",
    trait == "mot_fn" ~ "Motility (fast-none)",
    trait == "mot_nf" ~ "Motility (none-fast)",
    trait == "calc_l2h" ~ "Calcified (low-high)",
    trait == "calc_h2l" ~ "Calcified (high-low)",
    trait == "gen_l2h" ~ "Generalism (low-high)",
    trait == "gen_h2l" ~ "Generalism (high-low)",
    trait == "vuln_l2h" ~ "Vulnerability (low-high)",
    trait == "vuln_h2l" ~ "Vulnerability (high-low)"
  ))

# relabelling and filtering difference data for motif stats
diff_mots <- diffDat %>% filter(grepl("mot", metric)) %>% 
  mutate(trait2 = case_when(
    trait == "rand" ~ "Random",
    trait == "size_b2s" ~ "Size (Large-Small)",
    trait == "size_s2b" ~ "Size (Small-Big)",
    trait == "tier_i2p" ~ "Tiering (infaunal-pelagic)",
    trait == "tier_p2i" ~ "Tiering (pelagic-infaunal)",
    trait == "mot_fn" ~ "Motility (fast-none)",
    trait == "mot_nf" ~ "Motility (none-fast)",
    trait == "calc_l2h" ~ "Calcified (low-high)",
    trait == "calc_h2l" ~ "Calcified (high-low)",
    trait == "gen_l2h" ~ "Generalism (low-high)",
    trait == "gen_h2l" ~ "Generalism (high-low)",
    trait == "vuln_l2h" ~ "Vulnerability (low-high)",
    trait == "vuln_h2l" ~ "Vulnerability (high-low)"
  ))

# boundaries for visualisation
bounds_nets <- diffBounds %>% filter(metric != "TSS") %>% 
  filter(!grepl("mot", metric))
bounds_mots <- diffBounds %>% filter(grepl("mot", metric))


# ### Plots Format 2 (not used) - difference barplots with boundaries ----
# 
# # plot difference with 20% max difference area'd
# diffPlot_nets <- ggplot(diff_nets, aes(x = trait, y = diffEst, 
#                                        fill = trait))+
#   geom_col()+
#   geom_hline(data = bounds_nets, aes(yintercept = 0+interval), col = "blue")+
#   geom_hline(data = bounds_nets, aes(yintercept = 0-interval), col = "blue")+
#   facet_wrap(~metric, scales = "free_y", ncol = 3)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))+
#   guides(fill = "none")
# 
# diffPlot_mots <- ggplot(diff_mots, aes(x = trait, y = diffEst, 
#                                        fill = trait))+
#   geom_col()+
#   geom_hline(data = bounds_mots, aes(yintercept = 0+interval), col = "blue")+
#   geom_hline(data = bounds_mots, aes(yintercept = 0-interval), col = "blue")+
#   facet_wrap(~metric, scales = "free_y", ncol = 2)+
#   theme(axis.text.x = element_text(angle = 90))+
#   guides(fill = "none")
# 
# TSS_graph+(diffPlot_nets/diffPlot_mots)


## Plots version 3 (used) - Williams and Martinez Style ----
# plot differences from reference among sequences as dots 
# for each metric and identify closest trait order for each metric

### filter and select to find network metric and sequence for closest ----
nets_close <- diff_nets %>% 
  group_by(metric) %>% 
  filter(abs(diffEst - 0) == min(abs(diffEst - 0))) %>% 
  select(metric, trait, trait2)

### filter and select to find motif and sequence for closest ----
mots_close <- diff_mots %>% group_by(metric) %>% 
  filter(abs(diffEst - 0) == min(abs(diffEst - 0)))  %>% 
  select(metric, trait, trait2)

### create closest combined data frame ----
# manage labelling
all_close <- bind_rows(nets_close, mots_close)
all_diffs <- bind_rows(diff_nets, diff_mots) %>% 
  mutate(mets = factor(mets, levels = 
                         c("connectance", "max_tl_std","mean_tl_std","soi_std",
                           "diameter","mean_norm_btw",
                           "mean_normalized_degree", "sd_normalized_in_degree",
                           "sd_normalized_out_degree",
                           "norm_mot_ap_comp", "norm_mot_dir_comp",
                           "norm_mot_lin", "norm_mot_omn")))

xlabs <- c("Connectance", "Max Trophic Level","Mean Trophic Level", "Omnivory Index",
           "Diameter","Betweeness", "Degree", "In Degree","Out Degree", 
           "Apparent Competition Motif","Direct Competition Motif",
           "Linear Food Chain Motif","Omnivory Motif")

### Version 3 Plot (used) ----
ggplot(all_diffs, aes(x = mets, y = diffEst, 
                      group = trait2, colour = trait2))+
  geom_jitter(height = 0, width = 0.1, size =2, alpha = 0.75)+
  scale_x_discrete(labels = xlabs)+
  labs(x = "", y = "Difference from Empirical", 
       colour = "Sequence")+
  geom_text(data = all_close, aes(x = metric, y  = 0.5, 
                                  label = trait2),
            inherit.aes = FALSE, size = 2, angle = 45, 
            fontface = "italic")+
  geom_hline(yintercept = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        legend.position = 'top',
        legend.title = element_blank())


## Create Table of closest sequence for each metric ----

# max TSS score
TSS_res <- mm %>% filter(metric == "TSS") %>% 
  filter(meanVal == max(meanVal)) %>% 
  select(metric, trait)

# closest metrics
nets_res <- diff_nets %>% 
  group_by(metric) %>% 
  filter(abs(diffEst - 0) == min(abs(diffEst - 0))) %>% 
  select(metric, trait)

# closests motifs
mots_res <- diff_mots %>% group_by(metric) %>% 
  filter(abs(diffEst - 0) == min(abs(diffEst - 0)))  %>% 
  select(metric, trait)

#combine them
out_res <- bind_rows( nets_res, mots_res)  %>% 
  filter(metric %in% c("connectance", "max_tl",
                       "sd_normalized_in_degree",
                       "sd_normalized_out_degree",
                       "norm_mot_lin",
                       "norm_mot_omn",
                       "norm_mot_ap_comp",
                       "norm_mot_dir_comp"
  )) %>% 
  mutate(trait2 = case_when(
    trait == "rand" ~ "Random",
    trait == "size_b2s" ~ "Size (Large-Small)",
    trait == "size_s2b" ~ "Size (Small-Big)",
    trait == "tier_i2p" ~ "Tiering (infaunal-pelagic)",
    trait == "tier_p2i" ~ "Tiering (pelagic-infaunal)",
    trait == "mot_fn" ~ "Motility (fast-none)",
    trait == "mot_nf" ~ "Motility (none-fast)",
    trait == "calc_l2h" ~ "Calcified (low-high)",
    trait == "calc_h2l" ~ "Calcified (high-low)",
    trait == "gen_l2h" ~ "Generalism (low-high)",
    trait == "gen_h2l" ~ "Generalism (high-low)",
    trait == "vuln_l2h" ~ "Vulnerability (low-high)",
    trait == "vuln_h2l" ~ "Vulnerability (high-low)"
  ))

out_res

# generate table and write out which trait sequence is closest for each metric ----
# and how many tmes particular trait sequences are closest

out_res_report <- with(out_res, table(trait2)) %>% data.frame() %>% 
  filter(Freq>0) %>%
  arrange(desc(Freq))

out_res
out_res_report

# # uncomment to write
# write_csv(out_res, file = "metric-sequence_association.csv")
# write_csv(out_res_report, file = "closest_sequence_count.csv")