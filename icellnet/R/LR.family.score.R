#'
#' Distribution of class of molecules contribution to the communication scores
#'
#' @description Computes the contribution of families or subfamilies of molecules (included in my.selection) to the global communication scores.
#' If plot=FALSE, the function returns a table. If plot=T, returns a ggplot object, barplot representing contribution of families of molecules to the communication scores.
#'
#' @param lr Matrix of individual communication scores (sum on each column should be equal to the global communication score)
#' @param my.family Vector selecting family or subfamilies of molecules for the analysis
#' @param db.couple Output of the name.lr.couple() function. name.lr.couple(db, "Family") is set as a default
#' @param plot Logical, by default plot=F. If plot=T, function returns a ggplot object
#' @param title Title of the barplot (character)
#' @param ymax  Value as limit of the yaxis. If ymax=NULL, ymax takes the maximum value of the table
#'@param family.col Color vector for the family of molecules c("family1"= "color1", "family2"="color2").
#
#
#' @export
#' @examples
#' \dontrun{LR.family.score(lr=lr, my.family=my.family, db.couple=db.name.couple, plot=F)}
#'
LR.family.score = function(lr = lr,
                           my.family = my.family,
                           db.couple = as.data.frame(db.couple), plot=F, title=NULL, ymax=NULL, family.col=NULL) {

  lr=lr[stats::complete.cases(lr),]
  lr=lr[is.finite(log2(rowSums(lr))),]
  contribution = matrix(nrow = (length(my.family) + 1), ncol = length(lr[1, ]))
  colnames(contribution) = rownames(t(lr))
  rownames(contribution) = c(my.family, c("other"))
  for (family in my.family) {
    list = db.couple[which(db.couple[, 2] %in% family), 1]
    list = dplyr::intersect(list, rownames(lr))
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
  if (plot==F){
    return(contribution)
  }else
    table = contribution
  if (is.null(ymax)) {
    ymax = max(colSums(table)) + 1
  }
  melted <- reshape2::melt(table)
  melted$Var1 = as.factor(melted$Var1)
  if (is.null(family.col))
  {note("Colors automatically assigned. If you want to assign particular colors for fmaily of molecules, set it as family.col parameter")
    plot <- ggplot2::ggplot(data = melted,  ggplot2::aes(x = Var2, y = value, fill = Var1)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_brewer()  +
      ggplot2::labs(
        x = NULL,
        y = "score",
        fill = NULL ,
        title = title
      ) +
      ggplot2::guides(fill = guide_legend(reverse = TRUE)) +
      ggplot2::theme_classic() + ggplot2::ylim(0, ymax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  } else
    plot <- ggplot2::ggplot(data = melted,  ggplot2::aes(x = Var2, y = value, fill = Var1)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = family.col)+
    ggplot2::labs(
      x = NULL,
      y = "score",
      fill = NULL ,
      title = title
    ) +
    ggplot2::guides(fill = guide_legend(reverse = TRUE)) +
    ggplot2::theme_classic() + ggplot2::ylim(0, ymax) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  return(plot)
}

