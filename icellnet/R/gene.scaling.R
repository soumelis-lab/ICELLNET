
#' Gene scaling
#'
#' @description This function rescale each gene expression across condition to range from 0 to 10.
#'
#' @details Select genes from transcriptional profiles corresponding to the ligands and/or
#' receptors included in the database (db). Then, scale each ligand/receptor
#' gene expression among all the conditions ranging from 0 to 10.
#' For each gene:
#' - the maximum value (10) is defined as the mean expression of the 'n' highest values of expression.
#' - the minimum value (0) is defined as the mean expression of the 'n' lowest values of expression.
#' Default value of n is 1.
#'
#' @param data Matrix of transcriptomic profiles
#' @param n Number of highest expression values selected to define the maximum value of expression for the scaling
#' @param db Ligand/receptor database
#' @export
#' @examples
#' \dontrun{In this example, our matrix (CC.data) contains 80 transcriptional profiles.
#' The n is set to 4 in order to take the mean of the 5% extreme expression values as the maximum/minimum expression value.
#' mean(4 highest values) --> 10
#' mean(4 lowest values) --> 0
#' gene.scaling (data = CC.data, n=4, db = db)
#' }
#'

gene.scaling <-function(data = data,
                        n = n ,
                        db = db) {

  dim2 = dim(data)[2]
  dim1 = dim(data)[1]
  int = intersect(as.matrix(rownames(data)), as.matrix(db[, 1:5]))
  data = data[which(rownames(data) %in% int), ]
  for (i in 1:length(int)) {
    progressBar(i, max = length(int))
    sorted = sort(data[i, ], decreasing = TRUE)
    max = sum(sorted[1, 1:n]) / n
    min = sum(rev(sorted)[1, 1:n]) /n
    data[i, ] = (data[i, ] - min) / (max - min)* 10
    data[i,which(data[i,]>10)] <- 10
    data[i,which(data[i,]<0)] <- 0
  }
  data = (data[complete.cases(data),])

  return(data)
}
