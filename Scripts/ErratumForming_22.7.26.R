# Erratum - Error - Generate Extinction Sequences ----
# UPDATED July 2026 - Master Code to define and simulate Extinction Sequences ----
# (ORIGINAL 9 Sept 2022)

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

## Notes for update ----
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

library(tidyverse)
library(cheddar)

# data and functions ----
load("Data/ToarcianWebs_Guild_May2021.RData")
source("Scripts/NewMethod_Functions_update4publication.R")

## IMPORTANT: make sure all factors are coded with levels correctly ----
names(pre_meta_Guild)

# Part 1: Define and manage traits ----

## get current levels and adjust ----
pre_meta_Guild |> 
  select(where(is.factor)) |> 
  map(levels)

## adjusting pre_meta_Guild with orders and categories ----

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
  mutate(generality_cat = cut(generality,
                       breaks = c(-Inf, 0.5, 3.0, Inf),
                       labels = c("low", "medium", "high")))


# Part 2: Generate Sequences ----

reps <-10

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
                               group_by(desc(vulnerability_cat)) |>
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
genCatOrd_l2h <- 1:reps |> map(~pre_meta_Guild_use |>
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

# SECTION 3: Define function to drives primary extinction sequences with secondary allowed ----
# done across all extinction order scenarios

generate_seq <- function(extinctionOrder = randOrd){
  
  # collection zone for the final webs with 21 species
  SetOfWebs <- list()
  ExtinctionSequence <- list()
  Primary <- list()
  Secondary <- list()
  
  # loop first over reps
  set.seed(128)
  
  for(j in 1:reps){
    cat(paste(j, "\n"))
    
    # Sequence TRAIT
    orderExt <- extinctionOrder[[j]]
    
    # collecting the webs as species go extinct
    collect <- list()
    es <- vector()
    primaries <- vector()
    
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
      # This has the effect of propagating extinctions though the community - 
      # all consumers that are ultimately dependent upon all species in ‘remove’, 
      # and upon no other nodes (except themselves), will be removed.
      primaries[i] <- orderExt[next_idx] 
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
    
    # primary and secondary listing
    if(is.null(outWeb[[1]])){
      primary = NA
      secondary = NA
    } else
    {
      start <- NPS(preCom_Guild)$node
      primary <- primaries
      intermediate <- start[!(start %in% primary)]
      end <- NPS(outWeb[[1]])$node
      secondary <- start[!(start %in% primary)][!(intermediate %in% end)]
      
    }
    
    Primary[[j]] <- primary
    Secondary[[j]] <- secondary
    
  }
  # the output from the function is a list of webs and extinction sequences
  return(list(SetOfWebs, Primary, Secondary, ExtinctionSequence))
}

# SECTION 4: Apply generate_seq() function to ALL OF THE EXTINCTION ORDERS defined above ----
# NOTE: all but random are bi-directional.
# Note that extinction orders do NOT contain the basal speices

## random ----
randomExt <- generate_seq(extinctionOrder = randOrd)

## motility ----
motOrd_fast_nonExt <- generate_seq(extinctionOrder = motOrd_fast_non)
motOrd_non_fastExt <- generate_seq(extinctionOrder = motOrd_non_fast)

## tiering ----
tierExt_i2p <- generate_seq(extinctionOrder = tierOrd_i2p)
tierExt_p2i <- generate_seq(extinctionOrder = tierOrd_p2i)

## size ----
sizeExt_s2b <- generate_seq(extinctionOrder = sizeOrd_s2b)
sizeExt_b2s <- generate_seq(extinctionOrder = sizeOrd_b2s)

## calcification ----
calcExt_h2l <- generate_seq(extinctionOrder = calcOrd_h2l)
calcExt_l2h <- generate_seq(extinctionOrder = calcOrd_l2h)

## vulnerability ----
vulnExt_l2h <- generate_seq(extinctionOrder = vulnOrd_l2h)
vulnExt_h2l <- generate_seq(extinctionOrder = vulnOrd_h2l)

## vulnerability_cat ----
vulnCatExt_l2h <- generate_seq(extinctionOrder = vulnCatOrd_l2h)
vulnCatExt_h2l <- generate_seq(extinctionOrder = vulnCatOrd_h2l)

## generality ----
genExt_l2h <- generate_seq(extinctionOrder = genOrd_l2h)
genExt_h2l <- generate_seq(extinctionOrder = genOrd_h2l)

## generality categorical ----
genCatExt_l2h <- generate_seq(extinctionOrder = genCatOrd_l2h)
genCatExt_h2l <- generate_seq(extinctionOrder = genCatOrd_h2l)

# ## Check secondary extinction identities
# sort(table(unlist(tierExt_i2p[[3]])), decreasing = TRUE)
# sort(table(unlist(genExt_l2h[[3]])), decreasing = TRUE)
# sort(table(unlist(genExt_h2l[[3]])), decreasing = TRUE)

# par(mar = c(10,4,2,2), mfrow = c(1,3))
# plot(sort(table(unlist(tierExt_i2p[[3]])), decreasing = TRUE), las = 2, ylab = "count")
# title("Tiering Infaunal to Pelagic")
# plot(sort(table(unlist(genExt_l2h[[3]])), decreasing = TRUE), las = 2, ylab = "count")
# title("Generalism Low to High")
# plot(sort(table(unlist(genExt_h2l[[3]])), decreasing = TRUE), las = 2, ylab = "count")
# title("Generalism High to Low")

# remove webs where there was not 21 species network (labelled NULL in the sims) ----
# done for all extinction sequences

wrkWebs_rand <- Filter(function(x) !is.null(x), randomExt[[1]])

wrkWebs_motOrd_fast_non <- Filter(function(x) !is.null(x), motOrd_fast_nonExt[[1]])
wrkWebs_motOrd_non_fast <- Filter(function(x) !is.null(x), motOrd_non_fastExt[[1]])

wrkWebs_tier_i2p <- Filter(function(x) !is.null(x), tierExt_i2p[[1]])
wrkWebs_tier_p2i <- Filter(function(x) !is.null(x), tierExt_p2i[[1]])

wrkWebs_size_b2s <- Filter(function(x) !is.null(x), sizeExt_b2s[[1]])
wrkWebs_size_s2b <- Filter(function(x) !is.null(x), sizeExt_s2b[[1]])

wrkWebs_calc_h2l <- Filter(function(x) !is.null(x), calcExt_h2l[[1]])
wrkWebs_calc_l2h <- Filter(function(x) !is.null(x), calcExt_l2h[[1]])

wrkWebs_vuln_l2h <- Filter(function(x) !is.null(x), vulnExt_l2h[[1]])
wrkWebs_vuln_h2l <- Filter(function(x) !is.null(x), vulnExt_h2l[[1]])

wrkWebs_vulnCat_l2h <- Filter(function(x) !is.null(x), vulnCatExt_l2h[[1]])
wrkWebs_vulnCat_h2l <- Filter(function(x) !is.null(x), vulnCatExt_h2l[[1]])

wrkWebs_gen_l2h <- Filter(function(x) !is.null(x), genExt_l2h[[1]])
wrkWebs_gen_h2l <- Filter(function(x) !is.null(x), genExt_h2l[[1]])

wrkWebs_genCat_l2h <- Filter(function(x) !is.null(x), genCatExt_l2h[[1]])
wrkWebs_genCat_h2l <- Filter(function(x) !is.null(x), genCatExt_h2l[[1]])


# collect sequences to create master analysis data ----
wrkWebs_allSeqs <- list(rand = wrkWebs_rand,
                        mot_fn = wrkWebs_motOrd_fast_non, mot_nf = wrkWebs_motOrd_non_fast,
                        size_s2b = wrkWebs_size_s2b, size_b2s = wrkWebs_size_b2s,
                        tier_i2p = wrkWebs_tier_i2p, tier_p2i = wrkWebs_tier_p2i,
                        calc_h2l = wrkWebs_calc_h2l, calc_l2h = wrkWebs_calc_l2h,
                        vuln_l2h = wrkWebs_vuln_l2h, vuln_h2l = wrkWebs_vuln_h2l,
                        vulnCat_l2h = wrkWebs_vulnCat_l2h, vulnCat_h2l = wrkWebs_vulnCat_h2l, 
                        gen_l2h = wrkWebs_gen_l2h, gen_h2l = wrkWebs_gen_h2l,
                        genCat_l2h = wrkWebs_genCat_l2h, genCat_h2l = wrkWebs_genCat_h2l
                        )

# Final Step - Save for downstream ----
save(wrkWebs_allSeqs, file = "Data/wrkWebs_allSeqs_updateJuly2026.RData")
