#'
#' Perc_exp_infos function
#'
#' @description Compute the average expression for a gene of interest for one cluster (cell_id) and the percentage of cells expressing the gene in the cluster.
#'
#' @details Starting from a SeuratObject, this function returns as a vector: 1) the percentage of cells expressing the gene (counts >0), 2) the average gene expression for the cluster (slot=data from the SeuratObject)
#'
#' @param object SeuratObject with clusters as Idents(), and with slot counts and data.
#' @param gene Character, name of the gene of interest
#' @param assay Seurat assay, RNA assay by default.
#' @param cell_id Name/identity of the cluster, found also in Idents(object).
#' @export
#' @examples
#' \dontrun{Perc_exp_infos(object=object, gene=gene, cell_id=cell_id)}
#'


Perc_exp_infos <- function (object = object, assay = "RNA", gene = gene, cell_id = cell_id){
  WhichCells = colnames(object)[which(Idents(object) == cell_id)]
  if (class(object[[assay]]) == "Assay5") {
    if (is.null(object[[assay]]@layers$data)){
      data= object[[assay]]@layers$counts
    }else{data = object[[assay]]@layers$data}
    rownames(data) = rownames(object)
    colnames(data) = colnames(object)
    value.pos = sum(data[gene, WhichCells] > 0)/length(WhichCells)
    value.exp = mean(data[gene, WhichCells])
  }else{
    value.pos = sum(object[[assay]]@counts[gene, WhichCells] >
                      0)/length(WhichCells)
    value.exp = mean(object[[assay]]@data[gene, WhichCells])
  }
  return(c(value.pos, value.exp))
}
