####Load key packages----

##Load the pfwim package

library(pfwim)

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


##If using categorical size data, drop some of the variables
inferred_web_categoricalsizes <- infer_edgelist(data=ex_taxonlist,
                               col_taxon = "Guild", #indicate column containing taxon names
                               cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
                               certainty_req = "all", #state whether an interaction must meet all trait rules (e.g., feeding, motility, tiering, and size) or whether only a subset is require (e.g., if you input "2" then a taxon pair only has to meet feeding+motility or tiering+motility, but not all 4 feeding+motility+tiering+size)
                               return_full_matrix = FALSE, #indicate whether the function should return all taxon pairs with a value indicating how many trait rules were met, or whether the function should just return inferred interactions as an edgelist
                               print_dropped_taxa=TRUE #Print names of any taxa that have been removed
)
head(inferred_web_categoricalsizes)

##Alternatively, you could return all taxon-pairs without selecting only inferred interactions:
inferred_web_longlist <- infer_edgelist(ex_taxonlist,
                                        col_taxon = "Guild", #indicate column containing taxon names
                                        #col_num_size = "size_select", #if using a numerical prey-pred size rule indicate column containing numerical size values
                                        cat_combo_list = ex_traitrules, #point to dataframe containing trait rules
                                        return_full_matrix = TRUE #here indicating to return full matrix
)
head(inferred_web_longlist)
#When returning the full matrix, no taxa or taxon pairs will be removed
#NAs indicate that a rule wasn't met
#"Pres_sum" indicates how many rules were met

####Plotting food webs----

#Networks can be plotted in 3D, using trophic level to define position. Currently the easiest way to save an image is to screenshot it.

#Make a graph object from inferred web edgelist
graph_inferred_web<-igraph::graph_from_edgelist(inferred_web,directed=TRUE)

#Make 3D food web objects
fw3d<-make_3dfw(graph_inferred_web)

#The plot might take a few seconds to load
plot_3dfw(fw3d)

#By default trophic level is calculated using the standard method. Values are rounded into groups and scatter is added to help distinguish nodes
#These variables can be altered to change the 3D food web structure
plot_3dfw(make_3dfw(graph_inferred_web, round_tls = 1,scatter_tls=0,method="short-weighted"))

#You can view the rounded trophic level groups from the 3D food web object
print(head(fw3d[["trophic_levels"]]))

trophiclevels <- fw3d[["trophic_levels"]]

write.csv(trophiclevels, "./outputs/trophic levels/G1 - trophic levels.csv")

####Analyzing inferred webs----

#This package contains individual functions for calculating network-, node-, and motif-level metrics. It also contains functions to perform multiple analyses at once. Motif statistics can be reported separetly, but key motifs are incorporated into the network-level function. Look at the individual summary functions to see which metrics are reported:

#Make graph from inferred web edgelist
graph_inferred_web<-igraph::graph_from_edgelist(inferred_web,directed=TRUE)

#Network-level statistics
stats_inferred_web_network_level<-calc_network_stats(graph_inferred_web)
head(stats_inferred_web_network_level)

write.csv(stats_inferred_web_network_level, "./outputs/network stats/G1 - networks stats.csv")

#Node-level statistics
stats_inferred_web_node_level<-calc_node_stats(graph_inferred_web)
head(stats_inferred_web_node_level)

write.csv(stats_inferred_web_node_level, "./outputs/node stats/G1 - node stats.csv")

#All motif-level statistics (although most important are incorporated into network-level statistics)
stats_inferred_web_motif_level<-calc_motif_stats(graph_inferred_web)
head(stats_inferred_web_motif_level)

write.csv(stats_inferred_web_motif_level, "./outputs/motif stats/G1 - motif stats.csv")

#Calculate selected stats for Toarcian study
toarcian_stats<-calc_select_stats(graph_inferred_web)
head(toarcian_stats)

write.csv(toarcian_stats, "./outputs/toarcian stats/G1 - toarcian stats.csv")


####OPTIONAL: Generating functional diversity metrics----

#This package also contains functions to calculate functional diversity metrics. The function requires the trait list:

stats_functional_data<-calc_functional_stats(ex_taxonlist, cols=c("motility","tiering","feeding")) #Indicate which columns contain relevant functional information
head(stats_functional_data)
