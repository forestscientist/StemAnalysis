library(testthat)
library(StemAnalysis)

test_that("stem growth analysis", {

  path <- system.file(d_stem, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, stemgrowth = TRUE, stemdata = d_stem)
  result2 <- stemanalysism(xtree = 6, stemgrowth = TRUE, stemdata = d_stem)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(d_stem, package = "StemAnalysis")
  path1 <- system.file(d_parameters, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = d_stem,
                           parameterdata = d_parameters)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = d_stem,
                           parameterdata = d_parameters)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("estimation of tree biomass and carbon storage", {

  path <- system.file(d_stem, package = "StemAnalysis")
  path1 <- system.file(d_BEF, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, treecarbon = TRUE, stemdata = d_stem,
                           BEFdata = d_BEF)
  result2 <- stemanalysism(xtree = 6, treecarbon = TRUE, stemdata = d_stem,
                           BEFdata = d_BEF)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("height-diameter model", {

  path <- system.file(d_stem, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, HDmodel = TRUE, stemdata = d_stem)
  result2 <- stemanalysism(xtree = 6, HDmodel = TRUE, stemdata = d_stem)

  expect_type(result1, "list")
  expect_type(result2, "list")

})



test_that("all", {

  path <- system.file(d_stem, package = "StemAnalysis")
  path1 <- system.file(d_BEF, package = "StemAnalysis")
  path2 <- system.file(d_parameters, package = "StemAnalysis")

  result1 <- stemanalysism(xtree = 4, stemgrowth = TRUE, HDmodel = TRUE,
                           treecarbon = TRUE, stemdata = d_stem,
                           BEFdata = d_BEF, parameterdata = d_parameters)

  expect_type(result1, "list")

})

