library(testthat)
library(StemAnalysis)

test_that("stem growth analysis", {

  path <- system.file(stemdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, stemgrowth = TRUE, stemdata = stemdata)
  result2 <- stemanalysism(xtree = 6, stemgrowth = TRUE, stemdata = stemdata)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(stemdata, package = "StemAnalysis")
  path1 <- system.file(allomPardata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = stemdata,
                           allompardata = allomPardata)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = stemdata,
                           allompardata = allomPardata)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(stemdata, package = "StemAnalysis")
  path1 <- system.file(volumePardata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = stemdata,
                           volumepardata = volumePardata)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = stemdata,
                           volumepardata = volumePardata)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("height-diameter model", {

  path <- system.file(stemdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, HDmodel = TRUE, stemdata = stemdata)
  result2 <- stemanalysism(xtree = 6, HDmodel = TRUE, stemdata = stemdata)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("all", {

  path <- system.file(stemdata, package = "StemAnalysis")
  path1 <- system.file(volumePardata, package = "StemAnalysis")
  path2 <- system.file(allomPardata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, stemgrowth = TRUE, HDmodel = TRUE,
                           treecarbon = TRUE, stemdata = stemdata,
                           volumepardata = volumePardata, allompardata = allomPardata)

  expect_type(result1, "list")

})

