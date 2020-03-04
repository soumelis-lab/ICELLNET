#'
#' Barplot of family of molecules contribution to the communication scores
#'
#' @description Display the contribution of families or subfamilies of molecules (included in my.selection) to the global communication scores.
#' The matrix taken as an argument should be results of the LR.family.score() function.
#'
#' @param table Matrix of the contribution per family of molecules to the  global communication score (sum on each column should be equal to the global communication score)
#' @param title Title of the barplot
#' @param ymax Limits of the yaxis. If ymax=NULL, ymax takes the maximum value of the table.
#
#
#' @export
#' @examples
#' \dontrun{
#' LR.family.barplot(LR.family.score)
#' }
#'
LR.family.barplot = function (table = table,
                              title = NULL,
                              ymax = NULL) {
  if (is.null(ymax)) {
    ymax = max(colSums(table)) + 1
  }
  melted <- melt(table)
  melted$Var1 = as.factor(melted$Var1)
  plot <- ggplot(data = melted,  aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = family.col) +
    labs(
      x = NULL,
      y = "score",
      fill = NULL ,
      title = title
    ) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_classic() + ylim(0, ymax) +
    theme(
      axis.text.x = element_text(angle = 90, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    )
  return(plot)
}
