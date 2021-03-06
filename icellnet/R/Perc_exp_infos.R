#'
#' Perc_exp_infos function
#'
#' @description Compute the average expression for a gene of interest for one cluster (cell_id) and the percentage of cells expressing the gene in the cluster.
#'
#' @details Starting from a SeuratObject, this function returns as a vector: 1) the percentage of cells expressing the gene (counts >0), 2) the average gene expression for the cluster (slot=data from the SeuratObject)
#'
#' @param object SeuratObject with clusters as Idents(), and with slot counts and data.
#' @param gene Character, name of the gene of interest
#' @param cell_id Name/identity of the cluster, found also in Idents(object).
#' @export
#' @examples
#' \dontrun{Perc_exp_infos(object=object, gene=gene, cell_id=cell_id)}
#'


Perc_exp_infos <-function(object=object, gene=gene, cell_id=cell_id){
  WhichCells=colnames(object)[which(Idents(object)==cell_id)]
  value.pos=sum(object[['RNA']]@counts[gene,WhichCells]>0)/length(WhichCells)
  value.exp=mean(object[['RNA']]@data[gene,WhichCells])
  return(c(value.pos, value.exp))
}
