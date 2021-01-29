
#' Comparison of different communication score by displaying pvalues as a correlation plot
#'
#' @description Display pvalue matrix to compare different communication scores between the central cell and different peripheral cells
#'
#' @param pvalue pvalue matrix. See icellnet.score.pvalue() function for more information
#' @param PC Vector selecting a list of peripheral cell for the cell-cell communication analysis
#'
#' @export
#' @examples
#' \dontrun{
#' pvalue.matrix <- icellnet.score.pvalue(...)
#' pvalue.plot(pvalue=pvalue.matrix, PC=PC)
#' }
#'
#'
pvalue.plot <-function(pvalue, PC=PC){
  melted<- reshape2::melt(pvalue, na.rm = TRUE)
  melted$Var1=factor(melted$Var1, levels= rev(PC))
  melted$Var2=factor(melted$Var2, levels= PC)
  # Heatmap
  plot.pvalue = ggplot2::ggplot(data = melted, ggplot2::aes(Var2, Var1, fill = value))+
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradientn(colours = c("red","tomato","white","darkblue"),
                         values = c(0.00,0.05,0.1,1), name="p-value", space = "Lab", na.value = "grey50") +
    ggplot2::theme_minimal()+
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1), panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(vjust = 1, size = 12, hjust = 1))+
    ggplot2::coord_fixed()
  return(plot.pvalue)
}
