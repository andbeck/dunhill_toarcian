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
