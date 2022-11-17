# A selection of scripts written by Jack Shaw
# 1. network, node and motif metrics
# 2. Food Web Plotting
# 3. Helper Scripts.

# network, node and motif metric fucntions ----
# calc_network_stats and calc_select_stats call helper functions below

calc_network_stats <- function(graph) {
  graph <- igraph::simplify(graph, remove.multiple = T, remove.loops = FALSE)
  adj_mat <- as.matrix(igraph::as_adjacency_matrix(graph))
  row <- c(
    diameter = calc_diam(adj_mat),
    connectance = calc_C(adj_mat),
    size = calc_S(adj_mat),
    mean_normalized_degree = calc_degree(adj_mat, normalized = TRUE, return_mean = TRUE),
    cc = calc_cc(adj_mat),
    mean_norm_btw = calc_betweenness(adj_mat, return_mean = TRUE, normalize = TRUE),
    mean_betweenness = calc_betweenness(adj_mat, return_mean = TRUE),
    mean_pathlength = calc_cpl(adj_mat),
    soi_std = calc_system_omnivory(adj_mat, method = "standard"),
    soi_sw = calc_system_omnivory(adj_mat, method = "short-weighted"),
    mean_oi_std = calc_mean_omnivory(adj_mat, method = "standard"),
    mean_oi_sw = calc_mean_omnivory(adj_mat, method = "short-weighted"),
    mean_tl_std = calc_mean_tl(adj_mat, method = "standard"),
    mean_tl_sw = calc_mean_tl(adj_mat, method = "short-weighted")
  )

  mots <- calc_motif_stats(graph)[c("s1", "s2", "s4", "s5", "norm_s1", "norm_s2", "norm_s4", "norm_s5")]
  names(mots) <- c("mot_lin", "mot_omn", "mot_ap_comp", "mot_dir_comp", "norm_mot_lin", "norm_mot_omn", "norm_mot_ap_comp", "norm_mot_dir_comp")
  degs <- degree_tiers(adj_mat)[c("perc_int", "perc_bas", "perc_top")]

  row <- c(row, mots, degs)

  return(row)
}


calc_select_stats <- function(graph) {
  graph <- igraph::simplify(graph, remove.multiple = T, remove.loops = FALSE)

  adj_mat <- as.matrix(igraph::as_adjacency_matrix(graph))

  row <- c(
    connectance = calc_C(adj_mat),
    mean_normalized_degree = calc_degree(adj_mat, normalized = TRUE, return_mean = TRUE),
    sd_normalized_in_degree = sd(calc_degree(adj_mat, normalized = TRUE, method = "in"), na.rm = TRUE),
    sd_normalized_out_degree = sd(calc_degree(adj_mat, normalized = TRUE, method = "out"), na.rm = TRUE),
    mean_tl_std = calc_mean_tl(adj_mat, method = "standard"),
    max_tl_std = max(calc_node_tl_std(adj_mat), na.rm = TRUE),
    soi_std = calc_system_omnivory(adj_mat, method = "standard"),
    size = calc_S(adj_mat),
    diameter = calc_diam(adj_mat),
    mean_norm_btw = calc_betweenness(adj_mat, return_mean = TRUE, normalize = TRUE)
  )

  mots <- calc_motif_stats(graph)[c("norm_s1", "norm_s2", "norm_s4", "norm_s5")]
  names(mots) <- c("norm_mot_lin", "norm_mot_omn", "norm_mot_ap_comp", "norm_mot_dir_comp")

  row <- c(row, mots)

  return(row)
}


calc_C <- function(fw) {
  if (length(unique(c(rownames(fw), colnames(fw)))) == 0) {
    C <- NA
  } else {
    C <- length(fw[fw > 0]) / (length(unique(c(rownames(fw), colnames(fw))))^2)
  }
  return(C)
}

calc_degree <- function(fw, method = NULL, normalized = FALSE, return_mean = FALSE) {
  if (return_mean == FALSE) {
    deg <- degree(graph_from_adjacency_matrix(fw), mode = method, normalized = normalized)
    return(deg)
  } else if (return_mean == TRUE) {
    deg <- degree(graph_from_adjacency_matrix(fw), mode = "in", normalized = normalized)
    return(mean(deg, na.rm = TRUE))
  } else {

  }
}


calc_node_tl_std <- function(fw) {
  TL <- NetIndices::TrophInd(as.matrix(fw))[, 1]
  names(TL) <- rownames(fw)

  return(TL)
}



calc_node_tl_sw <- function(web) {
  if (length(web) > 1) {

 

    web <- as.matrix(web)
    rownames(web) <- colnames(web)

    basal <- rownames(subset(web, apply(web, 2, sum) == 0) & apply(web, 1, sum) != 0)
    edge.list_web <- igraph::graph.adjacency(web, mode = "directed")
    paths_prey <- igraph::shortest.paths(
      graph = edge.list_web, v = igraph::V(edge.list_web),
      to = igraph::V(edge.list_web)[basal], mode = "in", weights = NULL, algorithm = "unweighted"
    )
    paths_prey[is.infinite(paths_prey)] <- NA
    shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm = TRUE)))
    in_deg <- apply(web, 2, sum) 
    out_deg <- apply(web, 1, sum) 


    # Shortest TL
    sTL <- 1 + shortest_paths 
    

    S <- dim(web)[1]

    # Creating the matrix
    short_TL_matrix <- matrix(NA, S, S)
    prey_ave <- rep(NA, S) 
    chain_ave <- rep(NA, S)
    #
    for (j in 1:S) {
      for (i in 1:S) {
        lij <- web[i, j] 
        prey_ave[j] <- 1 / in_deg[j] 
        short_TL_matrix[i, j] <- lij * sTL[i] 
     
      }
    }

    prey_ave[which(prey_ave == Inf)] <- 0
    prey_ave[which(prey_ave == -Inf)] <- 0

    short_TL_matrix[is.nan(short_TL_matrix)] <- 0
   

    short_TL_matrix[which(short_TL_matrix == Inf)] <- 0
    short_TL_matrix[which(short_TL_matrix == -Inf)] <- 0



    sumShortTL <- as.matrix(apply(short_TL_matrix, 2, sum))


    SWTL <- matrix(data = NA, nrow = S, ncol = 1)
    for (i in 1:S) {
      SWTL[i] <- 1 + (prey_ave[i] * sumShortTL[i])
    }

    SWTL <- SWTL[, 1]
    names(SWTL) <- rownames(web)

    return(SWTL)
  } else {
    return(NA)
  }
}

calc_node_tl_all <- function(fw, method = NULL) {
  if (is.null(method)) {
    all_tls <- cbind(
      `standard` = calc_node_tl_std(fw),
      `network` = calc_node_tl_nw(fw),
    
      `short-weighted` = calc_node_tl_sw(fw)
    )
  } else {
    if (method == "standard") {
      all_tls <- calc_node_tl_std(fw)
    } else if (method == "short-weighted") {
      all_tls <- calc_node_tl_sw(fw)

    
    } else if (method == "network") {
      all_tls <- calc_node_tl_nw(fw)
    } else {

    }
  }

  return(all_tls)
}


calc_mean_tl <- function(fw, method = "standard") {
  mean_tl <- mean(calc_node_tl_all(fw, method), na.rm = TRUE)
  return(mean_tl)
}


calc_node_omnivory <- function(fw, method = "standard") {
  TL <- calc_node_tl_all(fw, method)

  Tij <- t(fw)
  p <- diet_get(Tij)
  ncomp <- ncol(Tij)

  OI <- vector(length = nrow(fw))
  for (i in 1:ncomp) OI[i] <- sum((TL - (TL[i] - 1))^2 * p[i, ])

  names(OI) <- rownames(fw)

  return(OI)
}

calc_system_omnivory <- function(fw = NULL, method = "standard") {
  if (length(fw) > 1) {
    Tij <- t(fw)
    OI <- calc_node_omnivory(fw, method)

    ncomp <- ncol(fw)
    Q <- vector(length = ncomp)
    OlogQ <- vector(length = ncomp)

    for (i in 1:ncomp) Q[i] <- sum(fw[, i])

    for (i in 1:ncomp) OlogQ[i] <- OI[i] * log(Q[i])
    ## Don't know if this is ok, to replace NaN with 0
    OlogQ[is.nan(OlogQ)] <- 0

    SOI <- sum(OlogQ) / log(sum(Q))
  } else {
    SOI <- NA
  }
  return(SOI)
}

calc_mean_omnivory <- function(fw, method = "standard") {
  mean_omn <- mean(calc_node_omnivory(fw, method), na.rm = TRUE)
  return(mean_omn)
}


## Functions taken from NetIndices package
InternalNetwork <- function(Tij, # to-from
                            Import, # flow from external (colNr Tij)
                            Export) # flow to external (colNr Tij)
{
  if (is.character(Import)) {
    import <- which(colnames(Tij) %in% Import)
  } else {
    import <- Import
  }
  if (length(import) != length(Import)) {
    stop("Import not recognized")
  }
  if (is.character(Export)) {
    export <- which(rownames(Tij) %in% Export)
  } else {
    export <- Export
  }
  if (length(import) != length(Import)) {
    stop("Import not recognized")
  }
  if (length(export) != length(Export)) {
    stop("Export not recognized")
  }

  ##
  ## CHECK THE INPUT
  ##

  # Flow or Tij should be inputted
  if (is.null(Tij)) {
    stop("cannot calculate indices - Flow or Tij should be inputted")
  }

  ncomp <- ncol(Tij) - length(import)
  if (ncomp != nrow(Tij) - length(export)) {
    stop("cannot calculate indices - internal flow input matrix not square ")
  }

  # indices to elements of T that are internal
  iN <- setdiff(1:nrow(Tij), export) # internal rows    of Tij
  jN <- setdiff(1:ncol(Tij), import) # internal columns of Tij


  # Total internal flows, externals removed.
  Tint <- Tij
  if (!is.null(export)) {
    Tint <- Tint[-export, ]
  }
  if (!is.null(import)) {
    Tint <- Tint[, -import]
  }

  # Total flows, including flow to/from externals
  FlowFrom <- colSums(Tij)
  FlowTo <- rowSums(Tij)
  FlowFromC <- FlowFrom[jN] # just the total from internal compartments
  FlowToC <- FlowTo[iN]

  return(list(
    Tint = Tint, iN = iN, jN = jN,
    import = import, export = export,
    FlowFrom = FlowFrom,
    FlowTo = FlowTo,
    FlowFromC = FlowFromC,
    FlowToC = FlowToC
  ))
}



Diet <- function(Tint, # Calculates diet composition
                 Dead = NULL, # index from Dead to Tint
                 iN = 1:nrow(Tint)) {

  ## p matrix contains the diet composition of predator i
  IntFlowTo <- rowSums(Tint) # Total food ingested

  p <- Tint

  for (i in 1:ncol(Tint)) {
    p[i, ] <- Tint[i, ] / IntFlowTo[i]
  }

  p[is.na(p)] <- 0

  ## take into account dead matter; Dead refers to column/row in Tij
  ## N$iN maps from Tint to Tij

  if (!is.null(Dead)) {
    p[which(iN %in% Dead), ] <- 0
  }
  return(p)
}



diet_get <- function(Tij) {
  IntFlowTo <- rowSums(Tij)
  p <- Tij
  for (i in 1:ncol(Tij)) {
    p[i, ] <- Tij[i, ] / IntFlowTo[i]
  }
  p[is.na(p)] <- 0
  return(p)
}


calc_S <- function(fw) {
  if (length(unique(c(rownames(fw), colnames(fw)))) == 0) {
    S <- NA
  } else {
    S <- vcount(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(S)
}


calc_cc <- function(fw) {
  if (length(unique(c(rownames(fw), colnames(fw)))) == 0) {
    cc <- NA
  } else {
    cc <- igraph::transitivity(graph_from_adjacency_matrix(as.matrix(fw)), "global")
  }
  return(cc)
}


calc_cpl <- function(fw) {
  if (length(unique(c(rownames(fw), colnames(fw)))) == 0) {
    avpath <- NA
  } else {
    avpath <- igraph::average.path.length(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(avpath)
}


calc_diam <- function(fw) {
  if (length(unique(c(rownames(fw), colnames(fw)))) == 0) {
    diam <- NA
  } else {
    diam <- igraph::diameter(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(diam)
}

calc_betweenness <- function(fw, return_mean = FALSE, normalize = FALSE) {
  btwn <- igraph::betweenness(graph_from_adjacency_matrix(as.matrix(fw)), normalized = normalize)

  if (return_mean == FALSE) {
    return(btwn)
  } else {
    return(mean(btwn, na.rm = TRUE))
  }
}


calc_motif_stats <- function(graph) {
  graph <- igraph::simplify(graph, remove.multiple = T, remove.loops = FALSE)
  triad.count <- igraph::triad.census(graph)
  graph_size <- as.numeric(as.character(vcount(graph)))

  names(triad.count) <- c(
    "empty", "single", "mutual", "s5", "s4", "s1", "d4",
    "d3", "s2", "s3", "d8", "d2", "d1", "d5", "d7", "d6"
  )

  triad.count.norm <- triad.count / graph_size^2
  names(triad.count.norm) <- paste("norm", names(triad.count.norm), sep = "_")

  # "norm_mot_lin","norm_mot_omn","norm_mot_ap_comp","norm_mot_dir_comp")
  triad.count.join <- c(triad.count, triad.count.norm[c("norm_s1", "norm_s2", "norm_s4", "norm_s5")])

  return(triad.count.join)
}

niche_model <- function(x, C) {

  # x=vector of s and connectance
  if (length(x) > 1) {
    S <- x[1]
    C <- x[2]
  } else {
    S <- x
    C <- C
  }

  if (S > 1 & C > 0 & is.finite(S) & is.finite(C)) {
    n.i <- sort(runif(S), decreasing = FALSE)
    r.i <- suppressWarnings(rbeta(S, 1, ((1 / (2 * C)) - 1)) * n.i)
    c.i <- suppressWarnings(runif(S, r.i / 2, n.i))

    if (any(is.na(c.i)) == FALSE & any(is.na(r.i)) == FALSE & any(is.na(n.i)) == FALSE) {
      a <- matrix(0, nrow = S, ncol = S)
      for (i in 2:S) {
        for (j in 1:S) {
          if (n.i[j] > (c.i[i] - (0.5 * r.i[i])) & n.i[j] <
            (c.i[i] + 0.5 * r.i[i])) {
            a[j, i] <- 1
          }
        }
      }

      rownames(a) <- 1:nrow(a)
      colnames(a) <- 1:ncol(a)
    } else {
      a <- NULL
    }

    return(a)
  } else {
    return(NULL)
  }
}


get_sc <- function(graph) {
  adj_mat <- as_adjacency_matrix((graph))
  siz <- as.numeric(as.character(calc_S(adj_mat)))
  con <- as.numeric(as.character(calc_C(adj_mat)))

  return(c(siz, con))
}


niche_model_replicates <- function(S, C, reps, return_mean = FALSE) {
  run <- variable <- value <- NULL

  stat_list <- c()
  for (i in 1:reps) {
    mod <- igraph::graph_from_adjacency_matrix(niche_model(S, C), mode = "directed")
    sta <- calc_select_stats(mod)
    stat_list <- rbind(stat_list, sta)
  }
  stat_list <- as.data.frame(stat_list, row.names = FALSE)
  stat_list$run <- 1:nrow(stat_list)

  stat_list2 <- stat_list %>%
    pivot_longer(!run, names_to = "variable", values_to = "value")

  if (return_mean == TRUE) {
    stat_list2 <- stat_list2 %>%
      group_by(variable) %>%
      summarize(
        avg = mean(value, na.rm = TRUE),
        SD = stats::sd(value, na.rm = TRUE),
        min = stats::quantile(value, probs = 0.025, na.rm = TRUE), max = stats::quantile(value, probs = 0.975, na.rm = TRUE), med = stats::median(value, na.rm = TRUE)
      ) %>%
      mutate(run = "averaged across runs")
  } else {

  }

  return(stat_list2)
}



calc_ME <- function(value, model_lower, model_median, model_upper) {
  upper_dif <- model_upper - model_median
  lower_dif <- model_median - model_lower

  model_error <- ifelse(value > model_median,
    ((value - model_median) / upper_dif),
    ((value - model_median) / lower_dif)
  )

  return(model_error)
}


calc_node_stats <- function(graph) {
  graph <- igraph::simplify(graph, remove.multiple = T, remove.loops = FALSE)
  adj_mat <- as.matrix(igraph::as_adjacency_matrix(graph))
  row <- list(
    taxon = rownames(adj_mat),
    betweenness = calc_betweenness(adj_mat),
    norm_btw = calc_betweenness(adj_mat, normalize = TRUE),
    norm_degree_in = calc_degree(adj_mat, method = "in", normalized = TRUE),
    norm_degree_out = calc_degree(adj_mat, method = "out", normalized = TRUE),
    degree_in = calc_degree(adj_mat, method = "in", normalized = FALSE),
    degree_out = calc_degree(adj_mat, method = "out", normalized = FALSE),
    tl_std = calc_node_tl_std(adj_mat),
    tl_sw = calc_node_tl_sw(adj_mat),
    oi_std = calc_node_omnivory(adj_mat, method = "standard"),
    oi_sw = calc_node_omnivory(adj_mat, "short-weighted")
  )
  row2 <- as.data.frame(do.call(cbind, row), row.names = NA) %>% dplyr::mutate_at(dplyr::vars(-("taxon")), anac)

  return(row2)
}

calc_functional_stats <- function(fd_data, cols = c("motility", "tiering", "feeding")) {
  Model <- Param <- feeding <- fw <- motility <- tiering <- NULL

  test <- fd_data[, c(cols)]
  test_join <- test %>%
    mutate(
      motility = as.factor(motility),
      tiering = as.factor(tiering),
      feeding = as.factor(feeding)
    )


  b <- unlist(ecospace::calc_metrics(samples = test_join, increm = FALSE, m = 7) %>% dplyr::select(-Model, -Param))

  names(b) <- c("NumOfSpecies", "NumOfTraits", "MeanPairDist", "MaxDist", "TotalVar", "FRic", "FDiv", "FDisp", "FEve", "qual.FRic")

  # S = number of species
  # H = unique life habits
  # D = mean pairwise distance between species in trait space
  # M = maximum distance between species in trait space
  # V = total variance for each functional trait
  # FRic = functional richness (minimum convex hull)
  # FDiv = functional divergence, mean distance of species from centroid
  # FDis = functional dispersion, total deviance of species from primary circle
  # FEve = functional evenness
  # qual.FRic = proportion (quality) of total PCoA space used

  return(b)
}


infer_edgelist <- function(data,
                           col_taxon = "taxon",
                           col_num_size = NULL,
                           cat_combo_list = NULL,
                           cat_trait_types = NULL,
                           num_size_rule = function(res_size, con_size) {
                             ifelse(res_size <= con_size, 1, 0)
                           },
                           certainty_req = "all",
                           return_full_matrix = FALSE,
                           print_dropped_taxa = FALSE,
                           hide_printout = FALSE,
                           ...) {

  # Hide "no visible binding for global variable" comment
  . <- Var1 <- Var2 <- pres_sum <- size_consumer <- size_resource <- taxon_consumer <- taxon_resource <- trait <- trait_type <- trait_type_consumer <- trait_type_interaction <- trait_type_pres <- trait_value_pres <- NULL

  if (is.null(cat_combo_list)) {
    cat_combo_list <- pfwim::example_data_traitrules
  } else {

  }

  names(data)[names(data) == col_taxon] <- "taxon"
  col_taxon <- "taxon"

  trait_cats <- tolower(unique(c(cat_combo_list[, c("trait_type_resource")], cat_combo_list[, c("trait_type_consumer")])))
  trait_cats_cons <- tolower(unique(c(cat_combo_list[, c("trait_type_consumer")])))
  col_taxon <- tolower(col_taxon)

   if (is.null(cat_trait_types)) {
  } else {
    trait_cats <- cat_trait_types
    trait_cats_cons <- cat_trait_types
  }

  fd <- as.data.frame(data)
  colnames(fd) <- tolower(colnames(fd))


  # Run for all trait types
  fw_match_traits <- c()
  for (i in trait_cats_cons) {
    col_combo_list_mini <- cat_combo_list %>% dplyr::filter(trait_type_consumer == i)

    res_taxa <- fd %>% dplyr::select(c(col_taxon, unique(col_combo_list_mini[, c("trait_type_resource")])))
    con_taxa <- fd %>%
      dplyr::select(c(col_taxon, i)) %>%
      filter_all(any_vars(!is.na(.)))
    con_taxa <- con_taxa[Reduce(`&`, lapply(con_taxa, function(x) !(is.na(x) | x == ""))), ]

    res_taxa_longer <- res_taxa %>%
      dplyr::mutate_all(., as.character) %>%
      tidyr::pivot_longer(-dplyr::all_of(col_taxon), names_to = "trait_type", values_to = "trait") %>%
      dplyr::filter(trait_type %in% col_combo_list_mini[, c("trait_type_resource")])

    con_taxa_longer <- con_taxa %>%
      dplyr::mutate_all(., as.character) %>%
      tidyr::pivot_longer(-dplyr::all_of(col_taxon), names_to = "trait_type", values_to = "trait") %>%
      dplyr::filter(trait != "primary")

    crossing_taxa <- tidyr::crossing(as.data.frame(res_taxa_longer) %>% setNames(paste0(names(.), "_resource")), con_taxa_longer %>% setNames(paste0(names(.), "_consumer"))) %>%
      dplyr::left_join(col_combo_list_mini[, c("trait_type_resource", "trait_type_consumer", "trait_consumer")] %>% dplyr::mutate(trait_type_pres = 1), by = c("trait_type_resource", "trait_type_consumer", "trait_consumer"))
    crossing_taxa <- crossing_taxa %>%
      dplyr::filter(trait_type_pres == 1) %>%
      dplyr::distinct() %>%
      dplyr::left_join(col_combo_list_mini %>% dplyr::mutate(trait_value_pres = 1), by = c("trait_type_resource", "trait_type_consumer", "trait_consumer", "trait_resource"))
    crossing_taxa <- crossing_taxa %>%
      dplyr::filter(trait_value_pres == 1) %>%
      dplyr::mutate(trait_type_interaction = i)
    crossing_taxa <- crossing_taxa[, c("taxon_resource", "taxon_consumer", "trait_type_interaction", "trait_value_pres")]

    fw_match_traits <- rbind(fw_match_traits, crossing_taxa)
    if (hide_printout == FALSE) print(i)
  }

  # If feasibility of interaction also defined by a numerical size rule, then utilize
  if (is.null(col_num_size)) {

  } else {
    fw_size_pastes <- unique(paste(fd[, c(col_taxon)], fd[, c(col_num_size)], sep = "SEP"))
    fw_size_expand <- as.data.frame(expand.grid(fw_size_pastes, fw_size_pastes)) %>%
      tidyr::separate(Var1, c("taxon_resource", "size_resource"), sep = "SEP", convert = TRUE) %>%
      tidyr::separate(Var2, c("taxon_consumer", "size_consumer"), sep = "SEP", convert = TRUE) %>%
      dplyr::mutate(trait_value_pres = num_size_rule(size_resource, size_consumer), trait_type_interaction = col_num_size) %>%
      dplyr::filter(trait_value_pres == 1) %>%
      dplyr::select(taxon_resource, taxon_consumer, trait_value_pres, trait_type_interaction)
    fw_match_traits <- bind_rows(fw_match_traits, fw_size_expand)
    if (hide_printout == FALSE) print(col_num_size)
  }


  fw_match_traits2 <- fw_match_traits %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = trait_type_interaction, values_from = trait_value_pres)
  fw_match_traits2$pres_sum <- fw_match_traits2 %>%
    dplyr::select(all_of(unique(fw_match_traits$trait_type_interaction))) %>%
    rowSums(na.rm = TRUE)

  # How many trait combinations need to be matched in order for an interaction to be viable
  if (certainty_req == "all") {
    certainty_val <- length(unique(fw_match_traits$trait_type_interaction))
  } else {
    certainty_val <- certainty_req
    certainty_val <- ifelse(certainty_val > length(unique(fw_match_traits$trait_type_interaction)), length(unique(fw_match_traits$trait_type_interaction)), certainty_val)
  }


  # Return full matrix with trait feasibility indicated or not
  if (return_full_matrix == TRUE) {
    fw_match_traits3 <- fw_match_traits2 %>% dplyr::distinct()
    # print(paste(length(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer)))), "taxa dropped from web"))



    return(as.matrix(fw_match_traits3))
  } else {
    fw_match_traits3 <- fw_match_traits2 %>%
      dplyr::filter(pres_sum >= certainty_val) %>%
      dplyr::distinct()
    if (hide_printout == FALSE) print(paste(length(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer)))), "taxa dropped from web"))

    if (print_dropped_taxa == TRUE) {
      if (hide_printout == FALSE) print(setdiff(unique(fd$taxon), unique(c(fw_match_traits3$taxon_resource, fw_match_traits3$taxon_consumer))))
    } else {
    }

    return(as.matrix(fw_match_traits3[, c("taxon_resource", "taxon_consumer")]))
  }
}


sample_pdf <- function(M = 100,
                       y = 2.5,
                       func = function(r, M, y) exp(-r / (exp((y - 1) * (log(M) / (y))))),
                       n_samp = 100) {
  row <- c()
  for (i in 1:M) {
    row <- c(row, func(i, M, y))
  }
  row2 <- as.data.frame(row)
  # ggplot(row2)+geom_point(aes(x=1:M,y=row))

  ar <- sum(row)
  row <- row / ar

  return(sample(1:M, n_samp, replace = T, prob = row))
}


powerlaw_prey <- function(el,
                          n_samp = 50,
                          func = function(r, M, y) exp(-r / (exp((y - 1) * (log(M) / (y)))))) {
  con_node_node_name_inferred <- NULL

  edgelist <- as.data.frame(el)
  colnames(edgelist) <- c("res_node_node_name_inferred", "con_node_node_name_inferred")

  web_list <- lapply(1:n_samp, matrix, data = NA, nrow = 0, ncol = 2)

  for (i in unique(edgelist$con_node_node_name_inferred)) {
    min <- edgelist %>% dplyr::filter(con_node_node_name_inferred == i)
    min <- as.data.frame(min)
    rich <- as.numeric(as.character(length(unique(min[, c("res_node_node_name_inferred")]))))

    t <- sample_pdf(M = rich, n_samp = n_samp)

    names(t) <- rep(i, length(t))

    for (j in 1:length(t)) {

      # Currently taking upper bound
      min2 <- min %>% sample_n(min(rich, t[[j]]))



      web_list[[j]] <- rbind(web_list[[j]], min2)
    }
  }

  web_list <- lapply(web_list, as.matrix)
  web_list <- lapply(web_list, na.omit)



  return(web_list)
}


degree_tiers <- function(fw, exclude_basal = FALSE) {
  graph <- graph_from_adjacency_matrix(as.matrix(fw))

  indegrees <- degree(graph, mode = "in")
  outdegrees <- degree(graph, mode = "out")

  numBas <- length(indegrees[which(indegrees == 0)])
  numTop <- length(outdegrees[which(outdegrees == 0)])
  numInt <- vcount(graph) - numBas - numTop

  nonbas <- numTop + numInt

  if (exclude_basal == FALSE) {
    basal <- (numBas / vcount(graph)) * 100
    top <- (numTop / vcount(graph)) * 100
    int <- ((vcount(graph) - (numBas + numTop)) / vcount(graph)) * 100
    indices <- c(n_bas = numBas, n_top = numTop, perc_bas = basal, perc_top = top, perc_int = int)
  } else {
    prop <- (numInt / nonbas) * 100
    indices <- c(n_bas = numBas, n_top = numTop, n_int = numInt, perc_int = prop)
  }

  return(indices)
}


Calc_TSS_new <- function(mat) {
  mat <- as.data.frame(mat)
  colnames(mat) <- c("real", "predicted")
  lis <- as.data.frame(table(mat))


  a <- as.numeric(as.character(subset(lis, real == 1 & predicted == 1)[, c("Freq")]))
  b <- as.numeric(as.character(subset(lis, real == 0 & predicted == 1)[, c("Freq")]))
  c <- as.numeric(as.character(subset(lis, real == 1 & predicted == 0)[, c("Freq")]))
  d <- as.numeric(as.character(subset(lis, real == 0 & predicted == 0)[, c("Freq")]))

  # TSS=(ad−bc) ∕ ((a+c)(b+d)), where a is the proportion of true positives
  # (link both predicted and observed), b is the proportion of false positives
  # (link predicted but not observed), c is the proportion of false negatives
  # (link not predicted, but observed), and d is the proportion of true negatives
  # (link not predicted and not observed)

  Py <- a / (a + b) * 100
  Pn <- d / (c + d) * 100

  # scores range from 1 (perfect prediction) to −1 (inverted prediction)
  TSS <- (a * d - b * c) / ((a + c) * (b + d))

  TSSd <- as.data.frame(data.frame(a = a, b = b, c = c, d = d, TSS = TSS))

  print("a = T+ve (link pred and obs), b=F+ve (link pred not obs), c=F-ve (link not pred but obs), d=T-ve (link not pred not obs)")
  print("")
  print(TSSd)
  return(TSSd)
}








#### UTILS----


round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}


anac <- function(x) {
  a <- as.numeric(as.character(x))
  return(a)
}


mat2list <- function(mat) {
  mat <- as.matrix(mat)
  list <- data.frame(X1 = colnames(mat)[col(mat)], X2 = rownames(mat)[row(mat)], dist = c(mat))
  colnames(list) <- c("col", "row")
  return(list)
}


clean_date <- function(drop_punctuation = TRUE) {
  t <- format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S%Z")

  t <- ifelse(drop_punctuation == FALSE, t, gsub("[^0-9A-Za-z///' ]", "", t, ignore.case = TRUE))

  return(t)
}


jlookup <- function(old_values,
                    replacement_df,
                    matching_col,
                    new_values) {
  temp_df <- as.data.frame(old_values) %>%
    dplyr::left_join(replacement_df[, c(matching_col, new_values)], by = c("old_values" = matching_col))

  return(temp_df[, c(new_values)])
}

make_3dfw_scaled_radius <- function(graph, round_tls = 0.1, scatter_tls = 0.15, method = "standard") {

  # Inspired by #http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization

  Freq <- Var1 <- round_tl <- NULL

  tl <- data.frame(tl = calc_node_tl_all(as_adjacency_matrix(graph), method = method))
  tl$name <- rownames(tl)
  tl$round_tl <- round_any(tl[, 1], round_tls)
  tl <- tl %>% arrange(round_tl)

  tl2 <- tl %>%
    group_by(round_tl) %>%
    dplyr::summarise(Freq = n()) %>%
    as.matrix()

  lay <- c()
  for (i in 1:nrow(tl2)) {
    l <- cbind(circpos(tl2[i, "Freq"], r = tl2[i, "Freq"]),
      # Noise in vertical positioning of nodes
      runif(n = tl2[i, "Freq"], tl2[i, "round_tl"] - scatter_tls, tl2[i, "round_tl"] + scatter_tls),
      tl = as.character(tl2[i, "round_tl"])
    )
    lay <- rbind(lay, l)
  }
  # V3 is the height
  lay <- as.data.frame(lay)
  lay$name <- tl$name
  lay <- lay[match(V(graph)$name, lay$name), ]
  rownames(lay) <- NULL

  # Return list of trophic levels
  rownames(tl) <- NULL
  # Add node colors
  cols <- c(0.9, 1.5, 2, 2.5, 3, 4, 5, 6, 1000)
  values <- as.numeric(as.character(tl$round_tl))
  ii <- cut(values, breaks = cols, include.lowest = TRUE)
  # Green, yellow, dark blue, dark red,light blue, orange, brown, pink
  tl$node_col_val <- c("#7ABF6F", "#EDCF23", "#2F4182", "#AB193B", "#6FBDBC", "#D9912B", "#754D07", "#D48786")[ii]
  tl$node_col_name <- c("green", "yellow", "dark blue", "dark red", "light blue", "orange", "brown", "pink")[ii]


  return(list(fw = graph, layout = lay, trophic_levels = tl))
}

# Food Web Plotting ----
make_3dfw <- function(graph, round_tls = 0.1, scatter_tls = 0.15, method = "standard") {
  
  # Inspired by #http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
  
  Freq <- Var1 <- round_tl <- NULL
  
  tl <- data.frame(tl = calc_node_tl_all(as_adjacency_matrix(graph), method = method))
  tl$name <- rownames(tl)
  tl$round_tl <- round_any(tl[, 1], round_tls)
  tl <- tl %>% arrange(round_tl)
  
  tl2 <- tl %>%
    group_by(round_tl) %>%
    dplyr::summarise(Freq = n()) %>%
    as.matrix()
  
  lay <- c()
  for (i in 1:nrow(tl2)) {
    l <- cbind(circpos(tl2[i, "Freq"], r = tl2[i, "Freq"]),
               # Noise in vertical positioning of nodes
               runif(n = tl2[i, "Freq"], tl2[i, "round_tl"] - scatter_tls, tl2[i, "round_tl"] + scatter_tls),
               tl = as.character(tl2[i, "round_tl"])
    )
    lay <- rbind(lay, l)
  }
  # V3 is the height
  lay <- as.data.frame(lay)
  lay$name <- tl$name
  lay <- lay[match(V(graph)$name, lay$name), ]
  rownames(lay) <- NULL
  
  # Return list of trophic levels
  rownames(tl) <- NULL
  # Add node colors
  cols <- c(0.9, 1.5, 2, 2.5, 3, 4, 5, 6, 1000)
  values <- as.numeric(as.character(tl$round_tl))
  ii <- cut(values, breaks = cols, include.lowest = TRUE)
  # Green, yellow, dark blue, dark red,light blue, orange, brown, pink
  tl$node_col_val <- c("#7ABF6F", "#EDCF23", "#2F4182", "#AB193B", "#6FBDBC", "#D9912B", "#754D07", "#D48786")[ii]
  tl$node_col_name <- c("green", "yellow", "dark blue", "dark red", "light blue", "orange", "brown", "pink")[ii]
  
  
  return(list(fw = graph, layout = lay, trophic_levels = tl))
}


plot_3dfw <- function(list3dfw, plotmany = F, save.rgl = F, bgc = "white") {
  graph <- list3dfw[["fw"]]
  lay <- list3dfw[["layout"]]
  
  
  cols <- c(0.9, 1.5, 2, 2.5, 3, 4, 5, 6, 1000)
  # Fixed colors by trophic level
  values <- as.numeric(as.character(lay$tl))
  ii <- cut(values, breaks = cols, include.lowest = TRUE)
  # Green, yellow, dark blue, dark red,light blue, orange, brown, pink
  nodecols <- c("#7ABF6F", "#EDCF23", "#2F4182", "#AB193B", "#6FBDBC", "#D9912B", "#754D07", "#D48786")[ii]
  
  
  lay$tl <- NULL
  lay$name <- NULL
  lay <- lay %>% mutate_if(is.factor, as.character)
  lay <- lay %>% mutate_if(is.character, as.numeric)
  
  ## Change spacing height of fw (don't know how to change height)
  lay$V3 <- lay$V3^0.5
  lay$V3 <- lay$V3 * 50
  # ggplot(lay)+geom_point(aes(x=V3,y=y))
  
  lay <- as.matrix(lay)
  
  if (plotmany == F) {
    # Open 3d viewer
    # rgl.open()
    open3d()
    par3d(windowRect = c(300, 300, 1200, 1200))
    rgl.bg(color = bgc)
    igraph::rglplot(graph,
                    layout = lay,
                    vertex.size = 10,
                    edge.width = 1,
                    vertex.label = NA,
                    # 0 means no arrows, 1 means backward arrows, 2 is for forward arrows and 3 for both
                    edge.arrow.mode = 0,
                    vertex.color = nodecols
    )
    rgl.viewpoint(theta = 180, phi = 75, fov = 20)
    
    if (save.rgl == "pdf") {
      nam <- paste("3d_net_img_", paste(unlist(str_extract_all(as.character(Sys.time()), "[[:digit:]]+")), collapse = ""), ".pdf", sep = "")
      rgl.postscript(nam, fmt = "pdf")
      print(paste("Saved as", nam))
    } else if (save.rgl == "png") {
      nam <- paste("3d_net_img_", paste(unlist(str_extract_all(as.character(Sys.time()), "[[:digit:]]+")), collapse = ""), ".png", sep = "")
      rgl.snapshot(nam, fmt = "png")
      print(paste("Saved as", nam))
    } else {
      
    }
  } else {
    igraph::rglplot(graph,
                    layout = lay,
                    vertex.size = 10,
                    edge.width = 1,
                    vertex.label = NA,
                    # 0 means no arrows, 1 means backward arrows, 2 is for forward arrows and 3 for both
                    edge.arrow.mode = 0,
                    vertex.color = nodecols
    )
  }
}

## Helper Function for plotting food webs ----
circpos <- function(n, r = 1) { # Coordinates on a circle
  
  ran<-runif(1,min=0,max=1) #Shuffle starting location
  #ran<-0
  rad <- seq(0+ran, ((2 * pi)+ran), length.out = n + 1)[-1]
  x <- cos(rad) * r
  y <- sin(rad) * r
  ret<-cbind(x, y)
  return(ret)
}

