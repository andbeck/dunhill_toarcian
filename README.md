# dunhill_toarcian
Code to Reproduce Dunhill et al analysis of secondary extinctions in Toarcian food webs.

There are 5 files provided.
1. Utility Functions - contains several metric functions and other tools used in other scripts
2. GenerateExtinctionSequences - a. defines traits; b. defines random stratfied samples of traits from core network; c. defines function to implement secondary extinction cascades via primary extinctions from b using cheddar RemoveNodes() function; d. collects networks and extinction sequences.
3. AnalyseToarcianNetworks - estimates several structural metrics, motifs and True Skill Statistic from simulated networks.
4. ToarcianWebs_Guild_May2021.RData - core input data of Toarcian Network
5. wrkWebs_allSeqs.RData - output from GenerateExintctionSequences (can be recreated with set.seed()) and used (loaded) in AnalyseToarcianNetworks.

