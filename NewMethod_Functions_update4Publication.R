## Assorted Functions used in Dunhill et al 


# # function to get true skill stats after extinctions made
# # Hanssen–Kuipers discriminant
# # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2664.2006.01214.x
# # Anubhav MS with Owen

# structure Only (Anubhav and Owen Function)
# ss_sim and ss_real need to be predation matrices
# 
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

#Function to get proportion of species that are different
speciesChange <- function(truthMat, predictMat){
  
  truth_sp <- dimnames(truthMat)[[1]]
  predict_sp <- dimnames(predictMat)[[1]]
  
  # postCom_Guild$nodes$node == dimnames(PredationMatrix(postCom_Guild))[[1]]
  
  in_truth_in_predict <- predict_sp[predict_sp %in% truth_sp]
  in_truth_not_predict <- predict_sp[!predict_sp %in% truth_sp]
  in_predict_in_truth <- truth_sp[truth_sp %in% predict_sp] # same as 1
  in_predict_not_truth <- truth_sp[!truth_sp %in% predict_sp]
  
  numRetained <- length(in_truth_in_predict)
  numNew <- length(in_predict_not_truth)
  
  propRetained <- numRetained/length(truth_sp)
  propNew <- 1-propRetained
  
  return(list(retainedSpecies = retainedSpecies,
              newSpecies = newSpecies,
              propRetained = propRetained,
              propNew = propNew))
  
}

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

# function to get stats on networks with 21 species
# adjusted to get most of Jacks stats.
# includes use of True Skill stats
# modify with Owens TO DO
cheddarNetworks <- function(web){
  TRmax <-  max(TrophicHeight(web))
  sd_InDeg <-  sd(NormalisedTrophicVulnerability(web))
  med_InDeg <- median(NormalisedTrophicVulnerability(web))
  sd_OutDeg <-  sd(NormalisedTrophicGenerality(web))
  med_OutDeg <- median(NormalisedTrophicGenerality(web))
  
  return(c(TRMax = TRmax,
                    sd_InDeg = sd_InDeg, med_InDeg = med_InDeg, 
                    sd_OutDeg = sd_OutDeg, med_OutDeg = med_OutDeg))
}

# function to ID who is extinct
whoExtinct <- function(simWeb = web, original = preCom_Guild ){
  simNodes <- NPS(simWeb)$node
  originalNodes <- NPS(original)$node
  extinct <- originalNodes[!originalNodes %in% simNodes]
  return(extinct)
}



