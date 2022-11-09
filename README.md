## dunhill_toarcian
Code to Reproduce Dunhill et al analysis of secondary extinctions in Toarcian food webs.

There are two sets of analysis (2 x scripts folders) and one shared data folder

The project should be initiated with `dunhill_toarcian.Rproj`.

### Data Files
There are two types of data.

1. Raw Network Data to generate metrics and data from scenarios - these are .csv files
2. files for use in plotting/analyses - these are .RData files and can be loaded via the plotting/analysis scripts

    - ToarcianWebs_Guild_May2021.RData - core input data of Toarcian Network
    - wrkWebs_allSeqs.RData - output from GenerateExintctionSequences (can be recreated with set.seed()) and used (loaded) in AnalyseToarcianNetworks.

### Analysis of network structure through time

Produces Figure X, X and X

### Analysis of secondary exinction scenarios

1. `pfwim.R` - contains several functions to estimate network metrics and motifs
2. `GenerateExtinctionSequences.R` - 
    - defines traits; 
    - defines random stratfied samples of traits from core network; 
    - defines function to implement secondary extinction cascades via primary extinctions from b using cheddar RemoveNodes() function; 
    - collects networks and extinction sequences.
3. `AnalyseToarcianNetworks.R` - estimates several structural metrics, motifs and True Skill Statistic from simulated networks.  Produces Figures X, X and X.


