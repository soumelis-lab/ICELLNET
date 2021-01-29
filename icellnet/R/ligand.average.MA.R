#'
#' ligand.average.MA function
#'
#' @description Compute the average expression for all ligands expressed in dataset provided, only available for Affymetrix Human Genome U133 Plus 2.0 Array.
#'
#' @details This function returns the mean of expression of each ligand (included in the database db) of several transcriptomic profiles that should
#' correspond to the same biological condition. If the ligand is composed of several subunits, the function computes the geometric mean of expression
#' of the mean values obtained for each subunits.
#'
#' @param data A dataframe of transcriptomic profiles with gene names as rownames
#' @param db Ligand/receptor database
#' @export
#' @examples
#' \dontrun{ligand.average.MA(db=db, data = data)}


ligand.average.MA <-
  function(db = db,
           data = data) {
    SYMBOL = rownames(data)
    x.lg = vector(length = dim(db)[1])
    data = as.data.frame(data)
    # for Microarray expression data
    if (is.null(data$Symbol))
      data$Symbol = SYMBOL
    for (mol in seq(1, dim(db)[1])) {
      # In ligand-provider data
      if (db$`Ligand 1`[mol] %in% data$Symbol) {
        if (is.na(db$`Ligand 2`[mol])) {
          x.lg[mol] = mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data) !=
                                                                                              "Symbol")]))
        } else{
          x.lg[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data) !=
                                                                                                               "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Ligand 2`[mol]), which(colnames(data) != "Symbol")]))))
        }
      } else{
        #warning(paste0('Warning: ', db$`Ligand 1`[mol], ' is not found in data matrix'))
        x.lg[mol] = NA
      }
    }
    return(x.lg)
  }
