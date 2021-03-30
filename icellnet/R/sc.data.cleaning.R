
#' Filtering scRNAseq data and computing average expression per cluster
#'
#' @description This function computes average expression per cluster from a Seurat object and eventually filters genes to keep only genes expressed by a specific percentage of the cells in their respective cluster. It returns the average the gene expression matrix (filtered) with clusters in columns.
#'
#' @details This function starts from a Seurat object and will used object[['RNA']]@data to compute average expression for each gene and for each cluster.
#' The function will also generate a statistic file (csv format) displaying for each gene in each cluster :
#' 1) the percentage of cells expressing the gene for this specific cluster
#' 2) the gene average expression value, considering all cells from the cluster.
#' If this file is already computed by a previous use of the function, it will be used automatically. To regenerate the file agan, set force.file=TRUE.
#' If filter.perc is not NULL and set to 5, the gene expression matrix will be filtered to keep only gene filtering of the genes expressed in less than 5% of the cells at the cluster level. Other genes (filtered) will be set at 0 for the corresponding cluster(s).
#'
#' @param object Seurat object with cells category (as future partner cells as Idents())
#' @param filter.perc Number between 0 and 100 (percentage), used to filter gene expression in each cluster (see details)
#' @param db Ligand/receptor database
#' @param db Ligand/receptor database
#' @param save_file Logical, TRUE as default. If TRUE, it will save the file scRNAseq_statsInfo_for_ICELLNET.csv, that will be reused automatically if ever needed.
#' @param path path to save the file scRNAseq_statsInfo_for_ICELLNET.csv. Working directory used as default path.
#' @param force.file Logical, FALSE as default. If scRNAseq_statsInfo_for_ICELLNET.csv already exists, set force.file=T will regenerate automatically the intermediate statistics table for data filtering.
#' @export
#' @examples
#' \dontrun{
#' sc.data.cleaning(object=seurat, db = icellnet_db, filter.perc=NULL, save_file=T, path=NULL)
#' }
#'


sc.data.cleaning <- function (object=object, db = db, filter.perc=NULL, save_file=T, path=NULL, force.file=FALSE) {

  if (!file.exists(paste0(path, "scRNAseq_statsInfo_for_ICELLNET.csv")) | force.file==T){
    data=object[['RNA']]@data
    int = dplyr::intersect(as.matrix(rownames(data)), as.matrix(db[, 1:5]))
    #restrict only to gene restricted to matrix
    data=data[which(rownames(data)%in%int),]
    data.int=as.data.frame(expand.grid(rownames(data), unique(Idents(object)))) #create intermediate data frame
    colnames(data.int)=c("Symbol","Cell_ID")
    data.int$Perc_posCell=NA
    data.int$Mean_exp=NA

    #fill in the intermediate dataframe
    print("Filling in intermediate table: percentage of expressing cell per cluster per gene, and mean of expression")
    for (i in seq(1,dim(data.int)[1])){
      data.int[i,3:4]=Perc_exp_infos(object = object, gene=data.int[i,1], cell_id = data.int[i,2])
    }

    #Save matrix and histogramm of Perc_pos_cell
    if (save_file==TRUE){
      if (is.null(path)){
        path=getwd()
        print("Intermediate table were saved in the working directory as scRNAseq_statsInfo_for_ICELLNET.csv. Set path parameter to change the directory to save files.")
      }
      write.csv(data.int, file = paste0(path, "scRNAseq_statsInfo_for_ICELLNET.csv"))
      hist(data.int[,"Perc_posCell"])
      pdf(file = paste0(path, "scRNAseq_statsInfo_Perc_posCell_", Sys.Date(),".pdf"), width = 5, height = 5)
      hist(data.int[,"Perc_posCell"], 100)
      print("Intermediate table were saved as scRNAseq_statsInfo_for_ICELLNET.csv.")
      dev.off()
    }
  }else{
    note( paste0("Following file used as intermediate statistics table: ", path, "scRNAseq_statsInfo_for_ICELLNET.csv. Use force.file=T to regenerate this file"))
    data.int=utils::read.csv(paste0(path,"scRNAseq_statsInfo_for_ICELLNET.csv" ), header = T)
    data.int=data.int[,-1]
  }

  #Filter expression above cell percentage value (ex: 5 -> 5% cells expressing the gene)
  if (!is.null(filter.perc)){
    data.int=dplyr::filter(data.int, Perc_posCell > filter.perc/100)
    print("Filtering done")
  }

  average.cluster=reshape2::dcast(data.int[,-3], formula = Symbol~ Cell_ID, value.var = "Mean_exp", drop = T)
  rownames(average.cluster)=average.cluster$Symbol
  average.cluster[is.na(average.cluster)]<-0
  return (average.cluster)
}
