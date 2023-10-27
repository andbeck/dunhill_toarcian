##Load other packages that might be needed in this tutorial
library(tidyverse)
library(igraph)
library(rgl)
library(future.apply)

# robustness data 
load("robustnessMeans_500_0.5spread.Rdata")

# plot them
ggplot(all_rob, aes(x = perc_loss, y = meanRob, col = net))+
  geom_line()+geom_smooth()

# # Robustness Analysis ----------
# 
# Takes a few hours to generate the replicate robusteness data
#
# library(future.apply)
# plan(multisession, workers = 4)
# 
# # functions needed to make food webs and plot them
# source("Scripts/pfim_scripts.R")
# source("Scripts/robustness_gradient.R")
# 
# # Data sources ----
# # Data frame listing taxon names and life habits. 
# # Taxon names must correspond to interaction data, and trait column names and
# # levels must correspond to trait rules.
# 
# ex_taxonlist_g1<-read.csv("Data/G1_guilds.csv")
# ex_taxonlist_g2<-read.csv("Data/G2_guilds.csv")
# ex_taxonlist_g3<-read.csv("Data/G3_guilds.csv")
# ex_taxonlist_g4<-read.csv("Data/G4_guilds.csv")
# 
# #Four column dataframe listing rules for feasible resource-consumer trophic interactions. 
# # Follow column naming format when using your own files.
# ex_traitrules<-read.csv("Data/feeding_rules.csv")
# 
# # calculate guild richess (taxa number) ----
# 
# taxa_g1 <- dim(ex_taxonlist_g1)[1]
# taxa_g2 <- dim(ex_taxonlist_g2)[1]
# taxa_g3 <- dim(ex_taxonlist_g3)[1]
# taxa_g4 <- dim(ex_taxonlist_g4)[1]
# 
# # Generating inferred webs----
# 
# ##Infer the food web using g1 ----
# inferred_web_g1 <- infer_edgelist(data=ex_taxonlist_g1,
#                                   col_taxon = "Guild", #indicate column containing taxon names
#                                   #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
#                                   cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
#                                   cat_trait_types = NULL, #if you only want to use a subset of life habit columns enter column names as vector otherwise NULL
#                                   #num_size_rule = function(res_size, con_size) {ifelse(res_size <= con_size, 1, 0)}, #function defining pred-prey size rule, here as long as predator is equal to or larger than prey then an interaction is feasible
#                                   certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
#                                   return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
#                                   print_dropped_taxa=TRUE #Print names of any taxa that have been removed
# )
# 
# ##Infer the food web using g2 ----
# inferred_web_g2 <- infer_edgelist(data=ex_taxonlist_g2,
#                                   col_taxon = "Guild", #indicate column containing taxon names
#                                   #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
#                                   cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
#                                   cat_trait_types = NULL, #if you only want to use a subset of life habit columns enter column names as vector otherwise NULL
#                                   #num_size_rule = function(res_size, con_size) {ifelse(res_size <= con_size, 1, 0)}, #function defining pred-prey size rule, here as long as predator is equal to or larger than prey then an interaction is feasible
#                                   certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
#                                   return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
#                                   print_dropped_taxa=TRUE #Print names of any taxa that have been removed
# )
# 
# ##Infer the food web using g3 ----
# inferred_web_g3 <- infer_edgelist(data=ex_taxonlist_g3,
#                                   col_taxon = "Guild", #indicate column containing taxon names
#                                   #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
#                                   cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
#                                   cat_trait_types = NULL, #if you only want to use a subset of life habit columns enter column names as vector otherwise NULL
#                                   #num_size_rule = function(res_size, con_size) {ifelse(res_size <= con_size, 1, 0)}, #function defining pred-prey size rule, here as long as predator is equal to or larger than prey then an interaction is feasible
#                                   certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
#                                   return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
#                                   print_dropped_taxa=TRUE #Print names of any taxa that have been removed
# )
# 
# ##Infer the food web using g4 ----
# inferred_web_g4 <- infer_edgelist(data=ex_taxonlist_g4,
#                                   col_taxon = "Guild", #indicate column containing taxon names
#                                   #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
#                                   cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
#                                   cat_trait_types = NULL, #if you only want to use a subset of life habit columns enter column names as vector otherwise NULL
#                                   #num_size_rule = function(res_size, con_size) {ifelse(res_size <= con_size, 1, 0)}, #function defining pred-prey size rule, here as long as predator is equal to or larger than prey then an interaction is feasible
#                                   certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
#                                   return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
#                                   print_dropped_taxa=TRUE #Print names of any taxa that have been removed
# )
# 
# 
# # Plotting food webs----
# 
# # Networks can be plotted in 3D, using trophic level to define position. 
# # Currently the easiest way to save an image is to screenshot it.
# 
# ## Make a graph object from inferred web edgelist ----
# graph_inferred_web_g1 <- igraph::graph_from_edgelist(inferred_web_g1,directed=TRUE)
# graph_inferred_web_g2 <- igraph::graph_from_edgelist(inferred_web_g2,directed=TRUE)
# graph_inferred_web_g3 <- igraph::graph_from_edgelist(inferred_web_g3,directed=TRUE)
# graph_inferred_web_g4 <- igraph::graph_from_edgelist(inferred_web_g4,directed=TRUE)
# 
# ## Make 3D food web objects ----
# fw3d_g1 <- make_3dfw(graph_inferred_web_g1)
# fw3d_g2 <- make_3dfw(graph_inferred_web_g2)
# fw3d_g3 <- make_3dfw(graph_inferred_web_g3)
# fw3d_g4 <- make_3dfw(graph_inferred_web_g4)
# 
# ## SECONDARY EXTINCTION ----
# g1_graph <- graph_from_edgelist(inferred_web_g1)
# g2_graph <- graph_from_edgelist(inferred_web_g2)
# g3_graph <- graph_from_edgelist(inferred_web_g3)
# g4_graph <- graph_from_edgelist(inferred_web_g4)
# 
# 
# # define spread of losses
# spread <- seq(from = 1, to = 99, by = 0.5)
# 
# # use replicate on each web to get 100 random starts 
# # and sequences for each
# # perc_loss value
# # run on 6 (or whatever above) cores using future_apply 
# # package (3 times faster)
# 
# rob_g1 <- future_replicate(500, 
#                   robustness_gradient(g1_graph, spread = spread), 
#                   simplify = FALSE)
# 
# rob_g2 <- future_replicate(500, 
#                          robustness_gradient(g2_graph, spread = spread), 
#                          simplify = FALSE)
# 
# rob_g3 <- future_replicate(500, 
#                          robustness_gradient(g3_graph, spread = spread), 
#                          simplify = FALSE)
# 
# rob_g4 <- future_replicate(500, 
#                          robustness_gradient(g4_graph, spread = spread), 
#                          simplify = FALSE)
# 
# 
# # bind the list elements into a dataframe
# # group by perc_loss value
# # take the mean to plot
# 
# out_g1 <- do.call("rbind", rob_g1) |> 
#   data.frame() |> 
#   group_by(perc_loss) |> 
#   summarise(
#     meanRob = mean(robustness)
#   )
# 
# 
# out_g2 <- do.call("rbind", rob_g2) |> 
#   data.frame() |> 
#   group_by(perc_loss) |> 
#   summarise(
#     meanRob = mean(robustness)
#   )
# 
# out_g3 <- do.call("rbind", rob_g3) |> 
#   data.frame() |> 
#   group_by(perc_loss) |> 
#   summarise(
#     meanRob = mean(robustness)
#   )
# 
# out_g4 <- do.call("rbind", rob_g4) |> 
#   data.frame() |> 
#   group_by(perc_loss) |> 
#   summarise(
#     meanRob = mean(robustness)
#   )
# 
# 
# # combine it all
# net = factor(c("pre", "post", "early", "late"),
#              levels = c("pre", "post", "early", "late"))
# 
# all_rob <- bind_rows(out_g1, out_g2,
#                       out_g3, out_g4) |> 
#   mutate(net = rep(net,each = length(spread)))
# 
# ## save the data
# save(out_g1, out_g2, out_g3, out_g4, all_rob,
#      file = "robustnessMeans_500_0.5spread.Rdata")
#   
