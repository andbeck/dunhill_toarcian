#Need the following packages: igraph, dplyr
#”spread” is a vector of numbers (1-99) that you want to analyze Robustness at (e.g., R50)

robustness_gradient <- function(graph, spread=seq(from=1,to=99,by=5)) {
  
  primprods <- which(igraph::degree(graph, mode = "in") == 0)
  graph_size <- igraph::vcount(graph)
  
  r_val_df <- c()
  for (percent_loss in spread) {
    pfw_df <- c()
    pfw_val1 <- 100
    pfw_val2 <- 1
    pfw_graph <- graph
    pfw_el <- as.data.frame(igraph::as_edgelist(graph))
    colnames(pfw_el) <- c("resource", "consumer")
    
    while (pfw_val1 >= (percent_loss / 100)) {
      
      # name of i (primary extinction)
      tax_nam <- sample(unique(pfw_el[, c("resource")]), 1)
      # name of i's predators (j--secondary extinction)
      tax_prey <- pfw_el %>%
        dplyr::filter(resource == tax_nam) %>%
        pull(consumer)
      tax_prey <- setdiff(tax_prey, primprods)
      
      
      ###### List which predators are only attached to i (i.e., have no other food sources)
      temp <- suppressMessages(pfw_el %>% dplyr::filter(consumer %in% tax_prey) %>% dplyr::group_by(consumer) %>% dplyr::summarise(in_degree = n_distinct(resource)) %>% dplyr::filter(in_degree == 1))
      tax_prey <- temp %>% pull(consumer)
      
      to_remove <- c(tax_nam, tax_prey)
      
      pfw_el2 <- pfw_el %>% dplyr::filter(!resource %in% to_remove & !consumer %in% to_remove)
      
      pfw_graph2_size <- length(unique(unique(c(pfw_el2[, c("resource")], pfw_el2[, c("consumer")]))))
      
      total_percent_loss <- pfw_graph2_size / graph_size
      
      pfw_df <- rbind(pfw_df, cbind(step = pfw_val2, num_sec_ext = length(tax_prey), new_graph_size = pfw_graph2_size, graph_size_as_perc = total_percent_loss))
      
      pfw_val1 <- total_percent_loss
      pfw_val2 <- pfw_val2 + 1
      pfw_el <- pfw_el2
    }
    
    pfw_df2 <- as.data.frame(pfw_df) %>%
      mutate(og_graph_size = graph_size, num_total_ext = num_sec_ext + 1) %>%
      mutate(cumulative_species_loss = og_graph_size - new_graph_size)
    r_val_df <- rbind(r_val_df, cbind(robustness=max(pfw_df2$step / graph_size),perc_loss=percent_loss))
  }
  
  return(r_val_df)
}