library(testthat)
library(StemAnalysis)

test_that("stem growth analysis", {

  path <- system.file(stemdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, stemgrowth = TRUE, stemdata = stemdata)
  result2 <- stemanalysism(xtree = 6, stemgrowth = TRUE, stemdata = stemdata)

  expect_type(result1, "character")
  expect_type(result2, "character")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(stemdata, package = "StemAnalysis")
  path1 <- system.file(parameterdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = stemdata,
                           parameterdata = parameterdata)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = stemdata,
                           parameterdata = parameterdata)

  expect_type(result1, "character")
  expect_type(result2, "character")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(stemdata, package = "StemAnalysis")
  path1 <- system.file(BEFdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = stemdata,
                           BEFdata = BEFdata)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = stemdata,
                           BEFdata = BEFdata)

  expect_type(result1, "character")
  expect_type(result2, "character")

})



test_that("height-diameter model", {

  path <- system.file(stemdata, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, HDmodel = TRUE, stemdata = stemdata)
  result2 <- stemanalysism(xtree = 6, HDmodel = TRUE, stemdata = stemdata)

  expect_type(result1, "list")
  expect_type(result2, "list")

})


