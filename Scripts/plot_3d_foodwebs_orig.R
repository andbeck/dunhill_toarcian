#' Make Components Required For 3D Foodweb
#'
#' @param graph An igraph object
#' @param round_tls How much averaging at each trophic level (indicate decimal places, e.g., 0.001, 0.01, 0.1, 1, 10)
#' @param scatter_tls How much scatter is within each rounded trophic level
#' @param method Method for calculating trophic level (see calc_node_tl_all function for details)
#'
#' @return Returns a list containing food web information and 3D layout information
#' @export
#'
#' @importFrom igraph as_adjacency_matrix
#' @import dplyr
#' @examples
#' make_3dfw(igraph::graph_from_adjacency_matrix(rand_adj_mat()))
make_3dfw <- function(graph, round_tls = 0.1, scatter_tls = 0.15, method = "standard") {

  # Inspired by #http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization

  Freq <- Var1 <- round_tl <- NULL

  tl <- data.frame(tl = calc_node_tl_all(as_adjacency_matrix(graph), method = method))
  tl$name <- rownames(tl)
  tl$round_tl <- round_any(tl[, 1], round_tls)
  tl <- tl %>% arrange(round_tl)

  tl2 <- tl %>%
    group_by(round_tl) %>%
    summarise(Freq = n()) %>%
    as.matrix()

  lay <- c()
  for (i in 1:nrow(tl2)) {
    l <- cbind(circpos(tl2[i, "Freq"], r = tl2[i, "round_tl"]),
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

  # lay[,c("tl")]<-as.numeric(as.character(lay[,c("tl")]))
  # hist(as.numeric(as.character(lay[,c("tl")])))
  # t<-as.data.frame(table(round_any(as.numeric(as.character(lay[,2])),0.1)))
  # ttt<-lay %>% filter(tl>4)
  # hist(as.numeric(as.character(ttt[,2])))

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



#' Plot 3D Food Web
#'
#' @param list3dfw Output from make_3dfw function (combined FW and layout information)
#' @param plotmany NEED TO FIGURE THIS OUT
#' @param save.rgl Indicate whether to save RGL image as PDF ("pdf"; takes a while) or PNG ("png") or not at all (F)
#' @param bgc Background color
#'
#' @return Plots 3D food web
#' @export
#'
#' @importFrom stringr str_extract_all
#' @importFrom igraph rglplot
#' @import dplyr
#' @import rgl


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


#' Plot Multiple 3D Food Web
#'
#' @param many_fw_list List containing numerous igraph objects
#'
#' @return Returns a frame showing many 3D food webs
#' @export
#'
#' @import rgl
plot_many_3dfw <- function(many_fw_list) {
  lenfw <- as.numeric(as.character(length(many_fw_list)))

  mfrow3d(nr = 2, nc = lenfw / 2, sharedMouse = TRUE, byrow = T)

  for (i in 1:lenfw) {
    plot_3dfw(make_3dfw(many_fw_list[[i]]), plotmany = T)
    if (i != lenfw) next3d()
    print(i)
  }
}

#' Layout Network Nodes As A Circle
#'
#' @param n Number of nodes to be around in circle
#' @param r Radius of circle
#'
#' @return Returns a frame showing coordinates for nodes in a circular arrangement
#' @export
#'
#' @importFrom stats runif
#' @examples
#' circpos(10, 1)
circpos <- function(n, r = 1) { # Coordinates on a circle

  ran <- runif(1, min = 0, max = 1) # Shuffle starting location
  # ran<-0
  rad <- seq(0 + ran, ((2 * pi) + ran), length.out = n + 1)[-1]
  x <- cos(rad) * r
  y <- sin(rad) * r
  ret <- cbind(x, y)
  return(ret)
}


# bbb<-make.3dfw(ttt)
# eee<-plot.3dfw(bbb)
#
# nam<-"Ythanjacob"
# many_fw_list<-list(emp=mj_web_list[[nam]][["empirical"]],
#                    nic=mj_web_list[[nam]][["niche"]],
#                    met=mj_web_list[[nam]][["metaweb"]],
#                    spe1=bc_webs[[nam]][[1]],
#                     spe2=bc_webs[[nam]][[1]],
#                   spe3=bc_webs[[nam]][[1]])
#
# ccc<-plot.many.3dfw(many_fw_list)
# rgl.close()
