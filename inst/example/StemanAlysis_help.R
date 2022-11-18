
library(StemAnalysis)

# Load the data sets
data(d_stem)
data(d_BEF)
data(d_parameters)

# To calculating tree growth profile and trends for an individual tree (for example, the tree number is 8) is needed
 stemanalysism(xtree = 8, stemgrowth = TRUE, stemdata = d_stem)

# To calculating tree carbon storage by allometric models is needed
 stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = d_stem,
               parameterdata = d_parameters)

# To calculating tree carbon storage by volume model is needed
 stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = d_stem,
               BEFdata = d_BEF)

# To fitting the height-diameter relationships
 stemanalysism(xtree = 8, HDmodel = TRUE, stemdata = d_stem)
