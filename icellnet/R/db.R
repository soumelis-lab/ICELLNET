#' db
#'
#' @name db
#'
#' @docType data
#'
#' @description Database of manually curated ligand-receptor interactions
#'
#' @format A dataframe with 373 rows and 13 colomuns:
#' ' \describe{
#'   \item{Ligand 1}{First ligand subunit}
#'   \item{Ligand 2}{Second ligand subunit}
#'   \item{Receptor 1}{First receptor subunit}
#'   \item{Receptor 2}{Second receptor subunit}
#'   \item{Receptor 3}{Third receptor subunit}
#'   \item{Alias}{Other known names for the ligand and/or the receptor}
#'   \item{Family}{Family of molecules}
#'   \item{Subfamily}{Further classification in a family, for cytokines.}
#'   \item{Classifications}{All kind of classifications in the same columns}
#'   \item{Source for interaction}{Reference}
#'   \item{PubMed ID}{PMID}
#'   \item{Comments}{Comments}
#' }
#' @source \url{http://biogps.org/dataset/species/human/}
#'
#
"db"
