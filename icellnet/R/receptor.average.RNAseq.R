#'
#' receptor.average.RNAseq function
#'
#' @description Compute the average expression for all receptors expressed in dataset provided, only available for RNA-seq dataset with gene symbol annotation.
#'
#' @details This function returns the mean of expression of each receptor (included in the database db) of several transcriptomic profiles that should
#' correspond to the same biological condition. If the receptor is composed of several subunits, the function computes the geometric mean of expression
#' of the mean values obtained for each subunits.
#'
#' @param data A dataframe of transcriptomic profiles with gene names as rownames
#' @param db Ligand/receptor database
#' @export
#' @examples
#' \dontrun{receptor.average.RNAseq (db=db, data = data)}
#'
receptor.average.RNAseq <-
  function(db = db,
           data = data) {
    SYMBOL = rownames(data)
    # check type of data
    x.rc = vector(length = dim(db)[1])
    data = as.data.frame(data)
    # for RNAseq expression data
    if (is.null(data$Symbol))
      data$Symbol = rownames(data)
    for (mol in seq(1, dim(db)[1])) {
      if (db$`Receptor 1`[mol] %in% data$Symbol) {
        # In receptor-provider data
        if (is.na(db$`Receptor 2`[mol]) & is.na(db$`Receptor 3`[mol])) {
          x.rc[mol] = mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                "Symbol")]))
        } else if (is.na(db$`Receptor 3`[mol])) {
          x.rc[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                                 "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Receptor 2`[mol]), which(colnames(data) != "Symbol")]))))
        } else if (is.na(db$`Receptor 2`[mol])) {
          x.rc[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                                 "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Receptor 3`[mol]), which(colnames(data) != "Symbol")]))))
        } else{
          x.rc[mol] = psych::geometric.mean(c(
            mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                      "Symbol")])),
            mean(data.matrix(data[which(data$Symbol ==
                                          db$`Receptor 3`[mol]), which(colnames(data) != "Symbol")])),
            mean(data.matrix(data[which(data$Symbol ==
                                          db$`Receptor 2`[mol]), which(colnames(data) != "Symbol")]))
          ))
        }
      } else{
        #warning(paste0('Warning: ', db$`Ligand 1`[mol], ' is not found in data matrix'))
        x.rc[mol] = NA
      }
    }
    return(x.rc)
  }
