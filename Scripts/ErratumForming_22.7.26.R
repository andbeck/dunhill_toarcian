# Erratum - Error - Generate Extinction Sequences

library(tidyverse)
library(cheddar)

# ORIGINAL ------------------------------------------------

# data and functions ----
load("Data/ToarcianWebs_Guild_May2021.RData")
source("Scripts/NewMethod_Functions_update4publication.R")

# replicates of sequences of extinctions
# for each trait, we generate 50 unique
# stratified random sequences of plausible extinction orders
reps <- 1

# SEQUENCES ----

# define levels of motility and tiering
pre_meta_Guild_use <- pre_meta_Guild |>
  mutate(motility_fn = factor(motility, levels = c("fast", "facultative", "slow", "nonmotile", NA))) |>
  mutate(motility_nf = factor(motility, levels = c(NA, "nonmotile", "slow", "facultative", "fast"))) |>
  # i to p and p to i
  mutate(tiering_i2p = factor(tiering_simple, levels = c("infaunal", "epifaunal", "pelagic", "primary"))) |>
  mutate(tiering_p2i = factor(tiering_simple, levels = c("primary", "pelagic", "epifaunal","infaunal")))

# # testing arrangements
names(pre_meta_Guild_use)
pre_meta_Guild_use |> arrange(tiering_i2p) |> pull(tiering_i2p)
pre_meta_Guild_use |> arrange(tiering_p2i) |> pull(tiering_p2i)
pre_meta_Guild_use |> arrange(motility_fn) |> pull(motility_fn)
pre_meta_Guild_use |> arrange(motility_nf) |> pull(motility_nf)
pre_meta_Guild_use |> arrange(size) |> pull(size)

# Make Sets to Compare -----------
# create random sequence of extinctions for each replicate ----
set.seed(128) # necessary for replication with random sampling

## random order ----
randOrd <- 1:reps |> map(~sample(NonBasalNodes(preCom_Guild))) |> 
  set_names("randOrd")

## non random orders ----

### motility - fast-> non (1) ; non->fast (2) ----
set.seed(128)
motOrd_fast_non <- 1:reps |> map(~pre_meta_Guild_use |>
                                    arrange(motility_fn) |>
                                    group_by(motility_fn) |>
                                    sample_frac() |>
                                    filter(node != "BASAL NODE") |>
                                    pull(node)) |> 
  set_names("motOrd_fast_non")

set.seed(128)
motOrd_non_fast <- 1:reps |> map(~pre_meta_Guild_use |>
                                    arrange(motility_nf) |>
                                    group_by(motility_nf) |>
                                    sample_frac() |>
                                    filter(node != "BASAL NODE") |>
                                    pull(node)) |> 
  set_names("motOrd_non_fast")


### infaunal, epifaunal, pelagic, primary ----
set.seed(128)
tierOrd_i2p <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(tiering_i2p) |>
                                group_by(tiering_i2p) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("tierOrd_i2p")

set.seed(128)
tierOrd_p2i <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(tiering_p2i) |>
                                group_by(tiering_p2i) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("tierOrd_p2i")

### Size ----
set.seed(128)
sizeOrd_b2s <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(size) |>
                                group_by(size) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("sizeOrd_b2s")

set.seed(128)
sizeOrd_s2b <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(desc(size)) |>
                                group_by(size) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("sizeOrd_s2b")


### Vulnerability ----
set.seed(128)
vulnOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(vulnerability) |>
                                group_by(vulnerability) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("vulnOrd_l2h")

set.seed(128)
vulnOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(desc(vulnerability)) |>
                                group_by(vulnerability) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("vulnOrd_h2l")
  
### Generality ----
set.seed(128)
genOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                               arrange(generality) |>
                               group_by(generality) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("genOrd_l2h")

set.seed(128)
genOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                               arrange(desc(generality)) |>
                               group_by(generality) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("genOrd_h2l")

### Calcification ----
set.seed(128)
calcOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(calcification) |>
                                group_by(calcification) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("calcOrdh2l")

set.seed(128)
calcOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(desc(calcification)) |>
                                group_by(calcification) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("calcOrdl2h")
  

origCheck <- data.frame(motOrd_fast_non, motOrd_non_fast,
           tierOrd_i2p, tierOrd_p2i,
           sizeOrd_s2b, sizeOrd_b2s,
           genOrd_l2h, genOrd_h2l,
           vulnOrd_h2l, vulnOrd_l2h,
           calcOrd_l2h, calcOrd_h2l,
           randOrd)

# SUGGESTED ------------------------------------------------
## For categorical variables (size, calcification):
##Ordering was explicitly defined prior to sequence generation, similar to the approach used for tiering and motility.

ss_pre_meta_Guild_use <- pre_meta_Guild |>
  mutate(ss_motility_fn = factor(motility, levels = c("fast", "facultative", "slow", "nonmotile", NA))) |>
  mutate(ss_motility_nf = factor(motility, levels = c(NA, "nonmotile", "slow", "facultative", "fast"))) |>
  # i to p and p to i
  mutate(ss_tiering_i2p = factor(tiering_simple, levels = c("infaunal", "epifaunal", "pelagic", "primary"))) |>
  mutate(ss_tiering_p2i = factor(tiering_simple, levels = c("primary", "pelagic", "epifaunal","infaunal"))) |> 
  # size
  mutate(ss_size_b2s = factor(size, levels = c("gigantic", "large", "medium", "small", "tiny", NA)))  |> 
  mutate(ss_size_s2b = factor(size, levels = c(NA, "tiny", "small", "medium", "large", "gigantic"))) |>
  # calcification
  mutate(ss_calc_h2l = factor(calcification, levels = c("heavy", "moderate", "light", NA))) |>
  mutate(ss_calc_l2h = factor(calcification, levels = c(NA, "light", "moderate", "heavy")))

# reset seed
set.seed(128)

ss_randOrd <- 1:reps |> map(~sample(NonBasalNodes(preCom_Guild))) |> 
  set_names("randOrd")

# PASTED AND CHANGED
set.seed(128)
ss_motOrd_fast_non <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                   arrange(ss_motility_fn) |>
                                   group_by(ss_motility_fn) |>
                                   sample_frac() |>
                                   filter(node != "BASAL NODE") |>
                                   pull(node)) |> 
  set_names("ss_motOrd_fast_non")

set.seed(128)
ss_motOrd_non_fast <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                   arrange(ss_motility_nf) |>
                                   group_by(ss_motility_nf) |>
                                   sample_frac() |>
                                   filter(node != "BASAL NODE") |>
                                   pull(node)) |> 
  set_names("ss_motOrd_non_fast")


### infaunal, epifaunal, pelagic, primary ----
set.seed(128)
ss_tierOrd_i2p <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(ss_tiering_i2p) |>
                               group_by(ss_tiering_i2p) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_tierOrd_i2p")

set.seed(128)
ss_tierOrd_p2i <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(ss_tiering_p2i) |>
                               group_by(ss_tiering_p2i) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_tierOrd_p2i")

### Size ----
set.seed(128)
ss_sizeOrd_b2s <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(ss_size_s2b) |>
                               group_by(ss_size_s2b) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_sizeOrd_b2s")

set.seed(128)
ss_sizeOrd_s2b <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(desc(ss_size_b2s)) |>
                               group_by(ss_size_b2s) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_sizeOrd_s2b")

### Vulnerability ----
### use - instead of desc() (desc() was not working properly)
set.seed(128)
ss_vulnOrd_l2h <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                  arrange(vulnerability) |>
                                  group_by(vulnerability) |>
                                  sample_frac() |>
                                  filter(node != "BASAL NODE") |>
                                  pull(node)) |> 
  set_names("ss_vulnOrd_l2h")

set.seed(128)
ss_vulnOrd_h2l <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                  arrange(-vulnerability) |>
                                  group_by(-vulnerability) |>
                                  sample_frac() |>
                                  filter(node != "BASAL NODE") |>
                                  pull(node)) |> 
  set_names("ss_vulnOrd_h2l")

### Generality ----
### use - instead of desc() (desc() was not working properly)
set.seed(128)
ss_genOrd_l2h <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                 arrange(generality) |>
                                 group_by(generality) |>
                                 sample_frac() |>
                                 filter(node != "BASAL NODE") |>
                                 pull(node)) |> 
  set_names("ss_genOrd_l2h")

set.seed(128)
ss_genOrd_h2l <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                                 arrange(-generality) |>
                                 group_by(-generality) |>
                                 sample_frac() |>
                                 filter(node != "BASAL NODE") |>
                                 pull(node)) |> 
  set_names("ss_genOrd_h2l")

### Calcification ----
set.seed(128)
ss_calcOrd_h2l <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(ss_calc_h2l) |>
                               group_by(ss_calc_h2l) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_calcOrdh2l")

set.seed(128)
ss_calcOrd_l2h <- 1:reps |> map(~ss_pre_meta_Guild_use |>
                               arrange(desc(ss_calc_l2h)) |>
                               group_by(ss_calc_l2h) |>
                               sample_frac() |>
                               filter(node != "BASAL NODE") |>
                               pull(node)) |> 
  set_names("ss_calcOrdl2h")




ssCheck <- data.frame(ss_motOrd_fast_non, ss_motOrd_non_fast,
           ss_tierOrd_i2p, ss_tierOrd_p2i,
           ss_sizeOrd_s2b, ss_sizeOrd_b2s,
           ss_genOrd_l2h, ss_genOrd_h2l,
           ss_vulnOrd_h2l, ss_vulnOrd_l2h,
           ss_calcOrd_l2h, ss_calcOrd_h2l,
           ss_randOrd)


names(origCheck);names(ssCheck)

## YES Match - no change
bind_cols(pull(origCheck, motOrd_fast_non),
          pull(ssCheck, ss_motOrd_fast_non))
bind_cols(pull(origCheck, tierOrd_p2i),
          pull(ssCheck, ss_tierOrd_p2i))

## YES Match: changed in ss to match mot and tier type ordering
bind_cols(pull(origCheck, sizeOrd_s2b),
          pull(ssCheck, ss_sizeOrd_s2b))

## NO MATCH?
bind_cols(pull(origCheck, calcOrdl2h),
          pull(ssCheck, ss_calcOrdl2h))

## No Match: changed desc() to "-"
### pos should be the same
bind_cols(pull(origCheck, genOrd_h2l),
          pull(ssCheck, ss_genOrd_h2l))

bind_cols(pull(origCheck,genOrd_l2h),
          pull(ssCheck, ss_genOrd_l2h))

### use - instead of desc()
bind_cols(pull(origCheck,genOrd_h2l),
          pull(ssCheck, ss_genOrd_h2l))

