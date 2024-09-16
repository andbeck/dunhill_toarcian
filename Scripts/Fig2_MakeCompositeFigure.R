# Make Master Figure

library(patchwork)

# run figure generating code ----
## Produces TSS plot and mean difference plot ----
# TSS_graph
# DiffMeans13
source("Scripts/Fig2_ab_AnalyseToarcianNetworks_update4Publication.R")

# produces structure and motif plot ----
# StructPlot
# motifPlot2
source("Scripts/Fig2_cd_guild network metric plotting.R")

# combine ----
## Vertical (using) ----
(TSS_graph/diffMeans13)|(structPlot/motifPlot2)

# ## Horizontal ----
# (TSS_graph + diffMeans13)/(structPlot+ plot_spacer() + motifPlot2)