#'
#' Distribution of class of molecules contribution to the communication scores
#'
#' @description Computes the contribution of families or subfamilies of molecules (included in my.selection) to the global communication scores.
#' If plot=FALSE, the function returns a table. If plot=T, returns a ggplot object, barplot representing contribution of families of molecules to the communication scores.
#'
#' @param lr Matrix of individual communication scores (sum on each column should be equal to the global communication score)
#' @param family Vector selecting family or subfamilies of molecules for the analysis
#' @param db.couple output of the name.lr.couple() function. name.lr.couple(db) is used by default
#' @param plot NULL, "heatmap" or "barplot". If NULL (default), returns a table with numerical values. It returns a ggplot object otherwise.
#' @param title Title of the barplot (character)
#' @param ymax  Value as limit of the yaxis. If ymax=NULL, ymax takes the maximum value of the table
#' @param family.col Color vector for the family of molecules c("family1"= "color1", "family2"="color2"). Used only for the barplot representation, NULL by default.
#' @param title Title of the barplot (character)
#' @param cluster_rows  Argument for ComplexHeatmap heatmap
#' @param cluster_columns Argument for ComplexHeatmap heatmap
#
#
#' @export
#' @examples
#' \dontrun{LR.family.score(lr=lr, plot=F)}
#'
#'

LR.family.score <-function (lr = lr, family = NULL, db.couple = as.data.frame(db.couple),
                            plot = NULL, title = NULL, ymax = NULL, family.col = NULL,
                            cluster_rows = T, cluster_columns = T)
{
  lr = na.omit(lr)
  row_names = rownames(lr)
  coln = colnames(lr)
  if (is.null(family)) {
    family = unique(db.couple[, 2]) %>% na.omit() %>% as.vector()
  }
  colnames(lr) = coln
  contribution = matrix(nrow = (length(family) + 1), ncol = length(lr[1,
  ]), 0)
  colnames(contribution) = rownames(t(lr))
  rownames(contribution) = c(family, c("other"))
  for (f in family) {
    list = db.couple[which(db.couple[, 2] %in% f), 1]
    list = dplyr::intersect(list, rownames(lr))
    if (length(list) > 1) {
      lr_f = as.data.frame(lr[list, ])
      contribution[f, ] <- colSums(lr_f, na.rm = TRUE,
                                   dims = 1)
    }
    else if (length(list) == 1) {
      contribution[f, ] = as.numeric(lr[list, ])
    }
    else {
      contribution[f, ] = 0
    }
  }
  contribution["other", ] = colSums(lr, na.rm = TRUE, dims = 1) -
    colSums(contribution, na.rm = TRUE, dims = 1)
  if (is.null(plot)) {
    return(contribution)
  }
  else {
    contribution=contribution[which(rowSums(contribution)>0),]
    if (is.null(ymax)) {
      ymax = max(colSums(contribution)) + 1
    }
    if (plot == "barplot") {
      melted <- reshape2::melt(contribution)
      melted$Var1 = as.factor(melted$Var1)
      melted$Var2 = as.factor(melted$Var2)
      if (is.null(family.col)) {
        print("Colors automatically assigned. If you want to assign particular colors for family of molecules, set it manually (family.col argument) or add colors as for ggplot objects.")
        plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var2,
                                                            y = value, fill = Var1)) + ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(x = NULL, y = "score", fill = NULL,
                        title = title) + ggplot2::guides(fill = guide_legend(reverse = TRUE)) +
          ggplot2::theme_classic() + ggplot2::ylim(0,
                                                   ymax) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                                                                              hjust = 1, vjust = 0.5))
      }
      else plot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var2,
                                                               y = value, fill = Var1)) + ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = family.col) +
        ggplot2::labs(x = NULL, y = "score", fill = NULL,
                      title = title) + ggplot2::guides(fill = guide_legend(reverse = TRUE)) +
        ggplot2::theme_classic() + ggplot2::ylim(0, ymax) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5))
    }
    else if (plot == "heatmap") {
      contribution_norm = scale(contribution, center = F,
                                scale = colSums(contribution, na.rm = TRUE, dims = 1)) %>%
        as.matrix()
      column_ha = ComplexHeatmap::HeatmapAnnotation(score = ComplexHeatmap::anno_barplot(as.numeric(colSums(lr,
                                                                                                            na.rm = TRUE, dims = 1)), height = unit(2, "cm"),
                                                                                         gap = unit(0.25, "points")))
      plot <- ComplexHeatmap::Heatmap(as.matrix(contribution_norm),
                                      top_annotation = column_ha, cluster_rows = cluster_rows,
                                      cluster_columns = cluster_columns, name = "proportion",
                                      show_column_names = T, show_row_names = T, column_title = title,
                                      col = circlize::colorRamp2(c(0, max(contribution_norm)),
                                                                 c("white", "red")))
    }
    return(plot)
  }
}
