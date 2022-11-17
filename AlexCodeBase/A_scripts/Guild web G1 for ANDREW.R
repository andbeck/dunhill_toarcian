##Load other packages that might be needed in this tutorial
library(dplyr)
library(tidyverse)
library(igraph)
library(ggplot2)

#Data frame listing taxon names and life habits. Taxon names must correspond to interaction data, and trait column names and levels must correspond to trait rules.
ex_taxonlist<-read.csv("./data/G1_guilds.csv")

#Four column dataframe listing rules for feasible resource-consumer trophic interactions. Follow column naming format when using your own files.
ex_traitrules<-read.csv("./data/feeding_rules.csv")

####Generating inferred webs----

##Infer the food web using your taxon list
inferred_web <- infer_edgelist(data=ex_taxonlist,
                               col_taxon = "Guild", #indicate column containing taxon names
                               #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
                               cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
                               cat_trait_types = NULL, #if you only want to use a subset of life habit columns enter column names as vector otherwise NULL
                               #num_size_rule = function(res_size, con_size) {ifelse(res_size <= con_size, 1, 0)}, #function defining pred-prey size rule, here as long as predator is equal to or larger than prey then an interaction is feasible
                               certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
                               return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
                               print_dropped_taxa=TRUE #Print names of any taxa that have been removed
)
head(inferred_web)

write.csv(inferred_web, "./outputs/interactions/G1 - pre-extinction guild inferred web.csv")


####Plotting food webs----

#Networks can be plotted in 3D, using trophic level to define position. Currently the easiest way to save an image is to screenshot it.

#Make a graph object from inferred web edgelist
graph_inferred_web<-igraph::graph_from_edgelist(inferred_web,directed=TRUE)

#Make 3D food web objects
fw3d<-make_3dfw(graph_inferred_web)

#The plot might take a few seconds to load
plot_3dfw(fw3d)

####Analyzing inferred webs----

#This package contains individual functions for calculating network-, node-, and motif-level metrics. It also contains functions to perform multiple analyses at once. Motif statistics can be reported separetly, but key motifs are incorporated into the network-level function. Look at the individual summary functions to see which metrics are reported:

#Make graph from inferred web edgelist
graph_inferred_web<-igraph::graph_from_edgelist(inferred_web,directed=TRUE)


#Calculate selected stats for Toarcian study
toarcian_stats<-calc_select_stats(graph_inferred_web)
head(toarcian_stats)

write.csv(toarcian_stats, "./outputs/toarcian stats/G1 - toarcian stats.csv")


