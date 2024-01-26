## dunhill_toarcian
Code to Reproduce Dunhill et al analysis of secondary extinctions in Toarcian food webs, submitted to Nature Ecology and Evolution for review.

The project should be initiated with `dunhill_toarcian.Rproj`.
There is a Scripts and Data folder

### Data Files

    - ToarcianWebs_Guild_May2021.RData - core input data of Toarcian Network
    - wrkWebs_allSeqs.RData - output from GenerateExintctionSequences (can be recreated with set.seed()) and used (loaded) in AnalyseToarcianNetworks.
    - metrics_time2.csv - details among times slices on network metrics
    - metrics_time.csv - details among times slices on network metrics

### Helper Scripts/Function Scripts

1. `pfwim.R` - contains several functions to define food webs using Shaw et al method and estimate network metrics and motifs
2. `GenerateExtinctionSequences.R` - 
    - defines traits; 
    - defines random stratfied samples of traits from core network; 
    - defines function to implement secondary extinction cascades via primary extinctions from b using cheddar RemoveNodes() function; 
    - collects networks and extinction sequences.
3. `plot_3d_foodwebs_orig.R` - Shaw et al functions to plot pfmim based food webs in 3D
4. `robustness_gradient.R` - cacluate R_1-99 robustness following generalised method from Jonsson et al Oikos 124: 446–457, 2015

### Visualisation of networks in 3D - produces Fig 1

1.`MakeCompositeFigure.R`

### Analysis of secondary exinction scenarios - Produces Figure 2a,b.

1.`AnalyseToarcianNetworks.R` - estimates several structural metrics, motifs and True Skill Statistic from simulated networks.  

### Analysis of network structure through time - Produces Figure 2c,d.

1.`guild network metric plotting.R` accesses combined data generated by `GuildWebsBuildPlotAnalyse.R` in Reference-Scratch

### Make Composite Fig 2a-d

1.`MakeCompositeFigure.R` - uses patchwork in R to combine figures

### Robustness Evaluation Fig 4

1.`RobustnessTest_FineScale.R`
