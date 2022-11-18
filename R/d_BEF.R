#' Information about input variables in biomass expansion factor for volume model
#'
#' A dataset containing the list of input variables and their description.
#'     Note: If a user uses the StemAnalysis package to estimate a sampled tree
#'     carbon by volume model, a biomass expansion factor dataset should be
#'     provided. Owing the biomass expansion factor data varies among tree
#'     species, the users should input the biomass expansion factor data for
#'     the correspondingly tree species
#'
#' @format A data frame with 2 rows and 4 variables:
#' \describe{
#'   \item{WD}{The wood density, kg m-3}
#'   \item{BEF}{The Biomass Expasion Factor}
#'   \item{R}{The Root-to-Shoot Ratio}
#'   \item{Cconcentration}{The tree carbon concentration, kg C kg-1}
#' }

#' @examples
#' head(d_BEF)
"d_BEF"
