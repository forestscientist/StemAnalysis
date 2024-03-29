% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stemdata.R
\docType{data}
\name{stemdata}
\alias{stemdata}
\title{Information about input variables in stem analysis data}
\format{
A data frame with 98 rows and 18 variables:
\describe{
  \item{No}{The line number}
  \item{Treeno}{The tree number for the sampled tree, the same number represents the same tree}
  \item{TreeTH}{Tree total height, m}
  \item{stemheight}{The stem height that the cross-sectional discs were obtained, m}
  \item{stemage}{The age namely the number of growth rings of the cross-sectional disc, year}
  \item{Dwithbark}{The maximum diameter of the cross-sectional disc with bark, cm}
  \item{Dnobark0}{The maximum diameter of the cross-sectional disc without bark, cm}
  \item{Dnobark1}{The diameter for the 1th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark2}{The diameter for the 2th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark3}{The diameter for the 3th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark4}{The diameter for the 4th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark5}{The diameter for the 5th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark6}{The diameter for the 6th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark7}{The diameter for the 7th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark8}{The diameter for the 8th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark9}{The diameter for the 9th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark10}{The diameter for the 10th age class inner growth ring of the cross-sectional disc, cm}
  \item{Dnobark11}{The diameter for the 11th age class inner growth ring of the cross-sectional disc, cm}
}
}
\usage{
stemdata
}
\description{
A dataset containing the list of input variables and their description.
    Note: If a user uses the StemAnalysis package to analysis a very big
    tree, the number of inner growth rings that diameter measured for some
    cross-sectional discs may be more than 11, the Dnobark12, Dnobark13, and
    much more variables can be added, which also could successfully run.
}
\examples{
head(stemdata)
}
\keyword{datasets}
