
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
  melted<- melt(pvalue, na.rm = TRUE)
  melted$Var1=factor(melted$Var1, levels= rev(PC))
  melted$Var2=factor(melted$Var2, levels= PC)
  # Heatmap
  plot.pvalue =ggplot(data = melted, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradientn(colours = c("red","tomato","white","darkblue"),
                         values = c(0.00,0.05,0.1,1), name="p-value", space = "Lab", na.value = "grey50") +
    theme_minimal()+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), panel.grid = element_blank(), axis.text.y = element_text(vjust = 1, size = 12, hjust = 1))+
    coord_fixed()
  return(plot.pvalue)
}
