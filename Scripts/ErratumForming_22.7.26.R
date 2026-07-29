# Erratum - Error - Generate Extinction Sequences

library(tidyverse)
library(cheddar)

# data and functions ----
load("Data/ToarcianWebs_Guild_May2021.RData")
source("Scripts/NewMethod_Functions_update4publication.R")

# FIXED ------------------------------------------------
## Using group_by(x) for standard ordering and group_by(desc(x)) 
## for reverse ordering provides a clean, predictable syntax across 
## all variable types:
## Numeric gets recoded from (1,2,3) to (-3, -2, -1)
## Factors get recoded from Level 1 -> Last; Last level -> Level 1

## arrange is not necessary - it is overriden by group_by and sample_frac
## Moving filter(node != "BASAL NODE") before group_by() ensures you don't 
## spend random draws on nodes you end up dropping!)

## Continuous variables ARE problematic: SHOULD have categorised them all
## desc() works but
## while they are discrete classes by taxa, there 
## are many isolated species (e.g. only one belemite with one size)
## SO: consider using categories (n=4) via cut_number in dplyr
## or custom cut.  So that internal randomisation is better.

## set.seed(128) before every ordering

## IMPORTANT: make sure all factors are coded with levels correctly ----
names(pre_meta_Guild)

### get current levels and adjust ----
pre_meta_Guild |> 
  select(where(is.factor)) |> 
  map(levels)

### adjusting ----

### size is OK: defaults from LARGE to SMALL
### calcification is OK: defaults from HEAVY to LIGHT

pre_meta_Guild_use <- pre_meta_Guild |>
  # OMIT BASAL NODE HERE
  filter(node != "BASAL NODE") |> 
  # motility needs ordering
  mutate(motility = factor(motility, levels = c("fast", "facultative", "slow", "nonmotile", NA))) |>
  # tiering needs ordering
  mutate(tiering = factor(tiering_simple, levels = c("infaunal", "epifaunal", "pelagic", "primary"))) |> 
  # cut_number works well with vuln data
  mutate(vulnerability_cat = cut_number(vulnerability, 
                                        n = 3,
                                        labels = c("low", "medium", "high"))) |> 
  # need to use cut() with generality because 68% is single value
  mutate(generaity_cat = cut(generality,
                       breaks = c(-Inf, 0.5, 3.0, Inf),
                       labels = c("low", "medium", "high")))



## Generate Sequences
reps <- 1

# create random sequence of extinctions for each replicate ----
## random order ----
set.seed(128)
randOrd <- 1:reps |> map(~sample(NonBasalNodes(preCom_Guild))) |> 
  set_names("randOrd")

## non random orders: FACTORS- Motility, Tiering, Calcification ----

### MOTILITY - fast-> non (1) ; non->fast (2) ----
set.seed(128)
motOrd_fast_non <- 1:reps |> map(~pre_meta_Guild_use |>
                                   group_by(motility) |>
                                   sample_frac() |>
                                   pull(node)) |> 
  set_names("motOrd_fast_non")

set.seed(128)
motOrd_non_fast <- 1:reps |> map(~pre_meta_Guild_use |>
                                   group_by(desc(motility)) |>
                                   sample_frac() |>
                                   pull(node)) |> 
  set_names("motOrd_non_fast")

# # TEST: species sampled are at opposite ends
# bind_cols(motOrd_fast_non, motOrd_non_fast) |> 
#   head()
# bind_cols(motOrd_fast_non, motOrd_non_fast) |> 
#   tail()

### Tiering - infaunal to primary (1) ; primary to infaunal (2) ----
set.seed(128)
tierOrd_i2p <- 1:reps |> map(~pre_meta_Guild_use |>
                                   group_by(tiering) |>
                                   sample_frac() |>
                                   pull(node)) |> 
  set_names("tierOrd_i2p")

set.seed(128)
tierOrd_p2i <- 1:reps |> map(~pre_meta_Guild_use |>
                                   group_by(desc(tiering)) |>
                                   sample_frac() |>
                                   pull(node)) |> 
  set_names("tierOrd_p2i")

# # TEST: species sampled are at opposite ends
# bind_cols(tierOrd_i2p, tierOrd_p2i) |>
#   head()
# bind_cols(tierOrd_i2p, tierOrd_p2i) |>
#   tail()

### Calcification - heavy to light (1) ; light to heavy (2) ----
set.seed(128)
calcOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(calcification) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("calcOrd_h2l")

set.seed(128)
calcOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(desc(calcification)) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("calcOrd_l2h")

# # TEST: species sampled are at opposite ends
# bind_cols(calcOrd_h2l, calcOrd_l2h) |>
#   head()
# bind_cols(calcOrd_h2l, calcOrd_l2h) |>
#   tail()

### Size ----
set.seed(128)
sizeOrd_b2s <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(size) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("sizeOrd_b2s")

set.seed(128)
sizeOrd_s2b <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(desc(size)) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("sizeOrd_s2b")

# # TEST: species sampled are at opposite ends
# bind_cols(sizeOrd_b2s, sizeOrd_s2b) |>
#   head()
# bind_cols(sizeOrd_b2s, sizeOrd_s2b) |>
#   tail()

## NUMERIC VARIABLES: Generality, Vulnerability

### Vulnerability ----

set.seed(128)
vulnOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(vulnerability) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("vulnOrd_l2h")

set.seed(128)
vulnOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |> 
                               group_by(desc(vulnerability)) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("vulnOrd_h2l")

### Vulnerability Category ----

set.seed(128)
vulnCatOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(vulnerability_cat) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("vulnCatOrd_l2h")

set.seed(128)
vulnCatOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                               group_by(desc(vulnerabilit_cat)) |>
                               sample_frac() |>
                               pull(node)) |> 
  set_names("vulnCatOrd_h2l")

### Generality ----
set.seed(128)
genOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                              group_by(generality) |>
                              sample_frac() |>
                              pull(node)) |> 
  set_names("genOrd_l2h")

set.seed(128)
genOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                              group_by(desc(generality)) |>
                              sample_frac() |>
                              pull(node)) |> 
  set_names("genOrd_h2l")


### Generality Category ----

set.seed(128)
genOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
                              group_by(generality_cat) |>
                              sample_frac() |>
                              pull(node)) |> 
  set_names("genCatOrd_l2h")
set.seed(128)
genCatOrd_h2l <- 1:reps |> map(~pre_meta_Guild_use |>
                              group_by(desc(generality_cat)) |>
                              sample_frac() |>
                              pull(node)) |> 
  set_names("genCatOrd_h2l")