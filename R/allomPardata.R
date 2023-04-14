#' Information about input variables in parameters of allometric models
#'
#' A dataset containing the list of input variables and their description.
#'     Note: If a user uses the StemAnalysis package to estimate a sampled tree
#'     carbon by allometric models, a parameters dataset should be provided.
#'     Owing the allometric model varies among tree species, the users should
#'     input the parameter data of the allometric models for the correspondingly
#'     tree species.
#'
#' @format A data frame with 3 rows and 7 variables:
#' \describe{
#'   \item{X}{The line number}
#'   \item{DBH}{The tree diameter of breast height, cm}
#'   \item{tissues}{The tree tissues including aboveground and belowground}
#'   \item{a}{The parameter a in the allometric model ln(B)=ln(a)+b×ln(DBH)+c×ln(H)}
#'   \item{b}{The parameter b in the allometric model ln(B)=ln(a)+b×ln(DBH)+c×ln(H)}
#'   \item{c}{The parameter c in the allometric model ln(B)=ln(a)+b×ln(DBH)+c×ln(H)}
#'   \item{Cconcentration}{The carbon concentration in each tree tissues, kg C kg-1}
#' }

#' @examples
#' head(allomPardata)
"allomPardata"
