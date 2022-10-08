## Assorted Functions used in Dunhill et al 

# function to make true skill stats comparison ----
# Hanssen–Kuipers discriminant
# code from Gupta et al Simultaneously estimating food web connectance and structure with uncertainty
# https://doi.org/10.1002/ece3.8643
# see also https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2006.01214.x

# structure Only (Anubhav and Owen Function)
# ss_sim and ss_real need to be predation matrices

TSS_func_struc <- function(ExtWeb){
  if (is.null(masterSpecies)) {
    stop("masterSpecies List is missing")
  }
  
  ss_sim <- Community(nodes = masterSpecies,
                      properties = list(title = "TestWeb"),
                      trophic.links = ExtWeb$trophic.links) %>% 
    PredationMatrix()
  ss_real <- Community(nodes = masterSpecies,
                       properties = list(title = "RealWeb"),
                       trophic.links = postCom_Guild$trophic.links) %>% 
    PredationMatrix()
  
  
  a <- sum(ss_sim==1 & ss_real==1)
  b <- sum(ss_sim==1 & ss_real==0)
  c <- sum(ss_sim==0 & ss_real==1)
  d <- sum(ss_sim==0 & ss_real==0)
  
  TSS_func <- (a*d-b*c)/((a+c)*(b+d))
  dist_TSS_r <- TSS_func
  if(is.nan(dist_TSS_r)==TRUE){dist_TSS_r <- 1000}
  
  return(dist_TSS_r)
}

# Functions to get Motif and Network Stats ----
# ported from Shaw et al....

# Motifs S1, S2, S4, S5
jackMotifs <- function(web){
  gg <- igraph::graph_from_edgelist(as.matrix(web$trophic.links))
  return(pfwim::calc_motif_stats(gg)[17:20])
}

# Network Stats Jack
jackNetworks <- function(web){
  gg <- igraph::graph_from_edgelist(as.matrix(web$trophic.links))
  return(calc_select_stats(gg))
}

# function to ID who is extinct
whoExtinct <- function(simWeb = web, original = preCom_Guild ){
  simNodes <- NPS(simWeb)$node
  originalNodes <- NPS(original)$node
  extinct <- originalNodes[!originalNodes %in% simNodes]
  return(extinct)
}



