#'
#' Distribution of class of molecules contribution to the communication scores
#'
#' @description Computes the contribution of families or subfamilies of molecules (included in my.selection) to the global communication scores.
#' The function returns a table that can be visualized in a barplot using the LR.family.barplot() function.
#'
#' @param lr Matrix of individual communication scores (sum on each column should be equal to the global communication score)
#' @param my.family Vector selecting family or subfamilies of molecules for the analysis
#' @param db.couple Output of the name.lr.couple() function. name.lr.couple(db, "Family") is set as a default
#
#
#' @export
#' @examples
#' \dontrun{LR.family.score(lr=lr, my.family=my.family, db.couple=db.name.couple)}
#'
LR.family.score = function(lr = lr,
                           my.family = my.family,
                           db.couple = as.data.frame(db.couple)) {
  lr=lr[complete.cases(lr),]
  lr=lr[is.finite(log2(rowSums(lr))),]
  contribution = matrix(nrow = (length(my.family) + 1), ncol = length(lr[1, ]))
  colnames(contribution) = rownames(t(lr))
  rownames(contribution) = c(my.family, c("other"))
  for (family in my.family) {
    list = db.couple[which(db.couple[, 2] %in% family), 1]
    list = intersect(list, rownames(lr))
    if (length(list) > 1) {
      lr.family = lr[list, ]
      contribution[family, ] <-
        colSums(lr.family, na.rm = TRUE, dims = 1)
    } else if (length(list) == 1) {
      contribution[family, ] = lr[list, ]
    } else{
      contribution[family, ] = 0
    }
  }
  family = "other"
  contribution2 = matrix(ncol = length(lr[1, ]))
  contribution2 = colSums(lr, na.rm = TRUE, dims = 1) - colSums(contribution, na.rm =
                                                                  TRUE, dims = 1)
  contribution[family, ] = contribution2
  return(contribution)
}
