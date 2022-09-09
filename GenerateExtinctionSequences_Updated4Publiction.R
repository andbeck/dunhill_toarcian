# Master Code to define and simulate Extinction Sequences ----
# 9 Sept 2022

# Part 1 defines traits that can be targets for primary extinction
# Part 2 generates replicate plausible primary extinction scenario orders
# Part 3 defines function to simulate primary and secondary extinctions using cheddar library
# Part 4 implements function on plausible extinction orders
# Part 5 collects and saves the final 21 species networks and extinction sequences.

## libraries ----

# CRAN R libraries for this script.
# code makes use of purrr::map and dplyr functions throughout
# code makes use of cheddar library functions (RemoveNodes) for secondary extinction analysis
# code makes use of base R Filter() function and others

library(tidyverse)
library(cheddar)

# data and functions ----
load("./Data/ToarcianWebs_Guild_May2021.RData")
source("NewMethod_Functions_update4Pub.R")

# replicates of sequences of extinctions
# for each trait, we generate 50 unique
# stratified random sequences of plausible extinction orders
reps <- 50

### SEQUENCES ####

# define levels of motility and tiering
pre_meta_Guild_use <- pre_meta_Guild %>%
  mutate(motility_fn = factor(motility, levels = c("fast", "facultative", "slow", "nonmotile", NA))) %>%
  mutate(motility_nf = factor(motility, levels = c(NA, "nonmotile", "slow", "facultative", "fast"))) %>%
  # i to p and p to i
  mutate(tiering_i2p = factor(tiering_simple, levels = c("infaunal", "epifaunal", "pelagic", "primary"))) %>%
  mutate(tiering_p2i = factor(tiering_simple, levels = c("primary", "pelagic", "epifaunal","infaunal")))

# # testing arrangements
# pre_meta_Guild_use %>% arrange(tiering_i2p) %>% pull(tiering_i2p)
# pre_meta_Guild_use %>% arrange(tiering_p2i) %>% pull(tiering_p2i)

# random sequence of extinctions for each replicate ----
set.seed(128)

randOrd <- 1:reps %>% map(~sample(NonBasalNodes(preCom_Guild)))

# non random orders ----

# ORDER: Both Orders - fast-> non (1) ; non->fast (2)
motOrd_fast_non <- 1:reps %>% map(~pre_meta_Guild_use %>%
                           arrange(motility_fn) %>%
                           group_by(motility_fn) %>%
                           sample_frac() %>%
                           filter(node != "BASAL NODE") %>%
                           pull(node))

motOrd_non_fast <- 1:reps %>% map(~pre_meta_Guild_use %>%
                                    arrange(motility_nf) %>%
                                    group_by(motility_nf) %>%
                                    sample_frac() %>%
                                    filter(node != "BASAL NODE") %>%
                                    pull(node))


# ORDER: infaunal, epifaunal, pelagic, primary
tierOrd_i2p <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(tiering_i2p) %>%
                            group_by(tiering_i2p) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

tierOrd_p2i <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(tiering_p2i) %>%
                            group_by(tiering_p2i) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

# this is largest to smallest
sizeOrd_b2s <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(size) %>%
                            group_by(size) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

sizeOrd_s2b <- 1:reps %>% map(~pre_meta_Guild_use %>%
                                arrange(desc(size)) %>%
                                group_by(size) %>%
                                sample_frac() %>%
                                filter(node != "BASAL NODE") %>%
                                pull(node))

vulnOrd_l2h <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(vulnerability) %>%
                            group_by(vulnerability) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

vulnOrd_h2l <- 1:reps %>% map(~pre_meta_Guild_use %>%
                                arrange(desc(vulnerability)) %>%
                                group_by(vulnerability) %>%
                                sample_frac() %>%
                                filter(node != "BASAL NODE") %>%
                                pull(node))

genOrd_l2h <- 1:reps %>% map(~pre_meta_Guild_use %>%
                           arrange(generality) %>%
                           group_by(generality) %>%
                           sample_frac() %>%
                           filter(node != "BASAL NODE") %>%
                           pull(node))

genOrd_h2l <- 1:reps %>% map(~pre_meta_Guild_use %>%
                           arrange(desc(generality)) %>%
                           group_by(generality) %>%
                           sample_frac() %>%
                           filter(node != "BASAL NODE") %>%
                           pull(node))

calcOrd_h2l <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(calcification) %>%
                            group_by(calcification) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

calcOrd_l2h <- 1:reps %>% map(~pre_meta_Guild_use %>%
                            arrange(desc(calcification)) %>%
                            group_by(calcification) %>%
                            sample_frac() %>%
                            filter(node != "BASAL NODE") %>%
                            pull(node))

# This function drives primary extinction sequences with secondary allowed ----
# done across all extinction order scenarios

generate_seq <- function(extinctionOrder = randOrd){

  # collection zone for the final webs with 21 species
  SetOfWebs <- list()
  ExtinctionSequence <- list()

  # loop first over reps
  set.seed(128)

  for(j in 1:reps){
    cat(paste(j, "\n"))

    # Sequence TRAIT
    orderExt <- extinctionOrder[[j]]

    # collecting the webs as species go extinct
    collect <- list()
    es <- vector()

    # push to near full extinctoion,
    # go back and get the point at which it is 21.
    for (i in 1:30){

      # web index updating
      if(i == 1)
      {tmp_web <- preCom_Guild}

      # it is possible that the ith species has already gone extinct via secondary
      # so our index is the next TRUE in this match
      next_idx <- which(orderExt %in% NonBasalNodes(tmp_web))[1]
      cat(paste(i, "--", next_idx, "\n"))

      # now do the iterative removal with possible secondary
      # using cascade method in cheddar: a multistep version of ‘secondary’ is applied. 
      # This has the effect of propogating extinctions though the community - 
      # all consumers that are ultimately dependent upon all species in ‘remove’, 
      # and upon no other nodes (except themselves), will be removed.
      out <- RemoveNodes(tmp_web, remove = orderExt[next_idx], method = 'cascade')

      # add the web to the collection list (will generate the node ID below)
      collect[[i]] <- out
      es[i] <- length(NPS(out)$node)
      tmp_web <- out
    }

    # from that sequence, collect the web at 21 species
    outWeb <- Filter(function(x) NumberOfNodes(x) == 21, collect)
    ExtinctionSequence[[j]] <- es

    # deal with situation where sequence of extinctions did not produce a 21 species network
    if(is_empty(outWeb)){outWeb <- list(NULL)}

    # see it
    print(outWeb[[1]])

    # collect webs with 21 species to here
    SetOfWebs[[j]] <- outWeb[[1]]

  }
  # the output from the function is a list of webs and extinction sequences
  return(list(SetOfWebs, ExtinctionSequence))
}

# Apply generate_seq() function to ALL OF THE EXTINCTION ORDERS defined above ----
# all but randomExt are bi-directional.

randomExt <- generate_seq(extinctionOrder = randOrd)

# motility
motOrd_fast_nonExt <- generate_seq(extinctionOrder = motOrd_fast_non)
motOrd_non_fastExt <- generate_seq(extinctionOrder = motOrd_non_fast)

# tiering
tierExt_i2p <- generate_seq(extinctionOrder = tierOrd_i2p)
tierExt_p2i <- generate_seq(extinctionOrder = tierOrd_p2i)

# size
sizeExt_s2b <- generate_seq(extinctionOrder = sizeOrd_s2b)
sizeExt_b2s <- generate_seq(extinctionOrder = sizeOrd_b2s)

# vulnerability
vulnExt_l2h <- generate_seq(extinctionOrder = vulnOrd_l2h)
vulnExt_h2l <- generate_seq(extinctionOrder = vulnOrd_h2l)

# generality
genExt_l2h <- generate_seq(extinctionOrder = genOrd_l2h)
genExt_h2l <- generate_seq(extinctionOrder = genOrd_h2l)

# calcification
calcExt_h2l <- generate_seq(extinctionOrder = calcOrd_h2l)
calcExt_l2h <- generate_seq(extinctionOrder = calcOrd_l2h)

# remove webs where there was not 21 species network (labelled NULL in the sims) ----
# done for all extinction sequences

wrkWebs_rand <- Filter(function(x) !is.null(x), randomExt[[1]])

wrkWebs_motOrd_fast_non <- Filter(function(x) !is.null(x), motOrd_fast_nonExt[[1]])
wrkWebs_motOrd_non_fast <- Filter(function(x) !is.null(x), motOrd_non_fastExt[[1]])

wrkWebs_tier_i2p <- Filter(function(x) !is.null(x), tierExt_i2p[[1]])
wrkWebs_tier_p2i <- Filter(function(x) !is.null(x), tierExt_p2i[[1]])

wrkWebs_size_b2s <- Filter(function(x) !is.null(x), sizeExt_b2s[[1]])
wrkWebs_size_s2b <- Filter(function(x) !is.null(x), sizeExt_s2b[[1]])

wrkWebs_vuln_l2h <- Filter(function(x) !is.null(x), vulnExt_l2h[[1]])
wrkWebs_vuln_h2l <- Filter(function(x) !is.null(x), vulnExt_h2l[[1]])

wrkWebs_gen_l2h <- Filter(function(x) !is.null(x), genExt_l2h[[1]])
wrkWebs_gen_h2l <- Filter(function(x) !is.null(x), genExt_h2l[[1]])

wrkWebs_calc_h2l <- Filter(function(x) !is.null(x), calcExt_h2l[[1]])
wrkWebs_calc_l2h <- Filter(function(x) !is.null(x), calcExt_l2h[[1]])


# collect sequences ----
wrkWebs_allSeqs <- list(rand = wrkWebs_rand,
                        mot_fn = wrkWebs_motOrd_fast_non, mot_nf = wrkWebs_motOrd_non_fast,
                        size_s2b = wrkWebs_size_s2b, size_b2s = wrkWebs_size_b2s,
                        tier_i2p = wrkWebs_tier_i2p, tier_p2i = wrkWebs_tier_p2i,
                        vuln_l2h = wrkWebs_vuln_l2h, vuln_h2l = wrkWebs_vuln_h2l,
                        gen_l2h = wrkWebs_gen_l2h, gen_h2l = wrkWebs_gen_h2l,
                        calc_h2l = wrkWebs_calc_h2l, calc_l2h = wrkWebs_calc_l2h)

# write sequence list for use in "AnalyseToarcianNetworks" Script ---
save(wrkWebs_allSeqs, file = "wrkWebs_allSeqs.RData")


