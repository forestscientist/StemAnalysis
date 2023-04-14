
library(StemAnalysis)

# Load the data sets
data(stemdata)
data(volumePardata)
data(allomPardata)

# To calculating tree growth profile and trends for an individual tree (for example, the tree number is 8) is needed
 stemanalysism(xtree = 8, stemgrowth = TRUE, stemdata = stemdata)

# To calculating tree carbon storage by allometric models is needed
 stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = stemdata,
               allompardata = allomPardata)

# To calculating tree carbon storage by volume model is needed
 stemanalysism(xtree = 8, treecarbon = TRUE, stemdata = stemdata,
               volumepardata = volumePardata)

# To fitting the height-diameter relationships
 stemanalysism(xtree = 8, HDmodel = TRUE, stemdata = stemdata)
