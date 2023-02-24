# Make Master Figure

library(patchwork)

# run figure generating code ----
## Produces TSS plot and mean difference plot ----
# TSS_graph
# DiffMeans13
source("Scripts/AnalyseToarcianNetworks_update4Publication.R")

# produces structure and motif plot ----
# StructPlot
# motifPlot2
source("Scripts/guild network metric plotting.R")

# combine
(TSS_graph + diffMeans13)/(structPlot+ plot_spacer() + motifPlot2)
