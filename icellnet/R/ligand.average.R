#'
#' ligand.average function
#'
#' @description Computes the average expression for all ligands expressed in the dataset provided with gene symbol annotation.
#'
#' @details This function returns the mean of expression of each ligand (included in the database db) of several transcriptomic profiles that should
#' correspond to the same biological condition. If the ligand is composed of several subunits, the function computes the geometric mean of expression
#' of the mean values obtained for each subunits.
#'
#' @param data A dataframe of transcriptomic profiles with gene names as rownames
#' @param db Ligand/receptor database
#' @export
#' @examples
#' \dontrun{ ligand.average(db=db, data = data)}

ligand.average <- function(db = db, data = data) {
  data = as.data.frame(data)
  x.lg = vector(length = dim(db)[1])
  if (is.null(data$Symbol)) {
    data$Symbol = rownames(data)
  }
  for (mol in seq(1, dim(db)[1])) {
    # In ligand-provider data
    if (db$`Ligand 1`[mol] %in% data$Symbol) {
      if (!is.na(db$`Ligand 4`[mol])) {
        x.lg[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 2`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 3`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 4`[mol]), which(colnames(data) != "Symbol")]))))
      }else if (!is.na(db$`Ligand 3`[mol])){
        x.lg[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 2`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 3`[mol]), which(colnames(data) != "Symbol")]))))
      }else if (!is.na(db$`Ligand 2`[mol])) {
        x.lg[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data)!= "Symbol")])),
                                            mean(data.matrix(data[which(data$Symbol == db$`Ligand 2`[mol]), which(colnames(data) != "Symbol")]))))
      }else{
        x.lg[mol] = mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data) != "Symbol")]))
      }
    }else{
      x.lg[mol] = NA
    }
  }
  return(x.lg)
}
