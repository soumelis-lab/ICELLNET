#'
#' Plot individual communication scores of a specific L/R interactions for all cell pairs.
#'
#' @description  Plot individual communication scores of a specific L/R interactions for all cell pairs. Useful to see specificity of an interaction between two cell types. Relevant only for single cell datasets, deciphering communication between clusters. Returns either a ggplot object (if plot=TRUE) or a matrix.
#'
#'
#' @param data A dataframe of transcriptomic profiles with gene names as rownames
#' @param couple character, name of a interaction "L / R". If several subunits, refer to output of name.lr.couple() function for interaction names.
#' @param db Ligand/receptor database
#' @param plot logical FALSE/TRUE (TRUE by default). If TRUE, the function returns the graphical representation of communication score for each cell pairs as a heatmap. If FALSE, the communication score matrix is returned.
#'
#' @export
#' @examples
#' \dontrun{LR.viz(data=data.icell , couple=couple, db=db , plot=TRUE)}
#'
LR.viz <- function(data=data , couple=couple, db=db , plot=TRUE){

  #check format data
  SYMBOL = rownames(data)
  data = as.data.frame(data)
  if (is.null(data$Symbol)) {
    data$Symbol = SYMBOL
  }
  #format db and select interaction of interest
  rownames(db)=name.lr.couple(db=db, type="Family")[,1]
  db.sel=db[which(rownames(db)==couple),]
  if (dim(db.sel)[1]==0){
    print("ERROR - interaction not in the DB - correct interaction")
  }
  if (is.na(db.sel$`Ligand 2`)) {
    if (db.sel$`Ligand 1`%in% data$Symbol){
      lig=as.vector(data[which(data$Symbol == db.sel$`Ligand 1`) , which(colnames(data) !=  "Symbol")])
    } else {lig = NA}
  }else{
    if (db.sel$`Ligand 1`%in% data$Symbol & db.sel$`Ligand 2`%in% data$Symbol){
      lig=as.vector(apply(data[which(data$Symbol %in% c(db.sel$`Ligand 1`, db.sel$`Ligand 2`)), which(colnames(data) !=
                                                                                                        "Symbol")], 2, psych::geometric.mean) )}
    else {lig = NA}
  }

  if (is.na(db.sel$`Receptor 2`) & is.na(db.sel$`Receptor 3`)){
    if (db.sel$`Receptor 1`%in% data$Symbol){
      rec=as.vector(data[which(data$Symbol == db.sel$`Receptor 1`) , which(colnames(data) !=  "Symbol")])
    }else {rec = NA}

  } else if (is.na(db.sel$`Receptor 3`)){
    if (db.sel$`Receptor 1`%in% data$Symbol & db.sel$`Receptor 2`%in% data$Symbol ){
      rec=as.vector(apply(data[which(data$Symbol %in% c(db.sel$`Receptor 1`, db.sel$`Receptor 2`)), which(colnames(data) !=
                                                                                                            "Symbol")], 2, psych::geometric.mean))
    }else {rec = NA}

  } else{
    if (db.sel$`Receptor 1`%in% data$Symbol & db.sel$`Receptor 2`%in% data$Symbol ){
      rec=as.vector(apply(data[which(data$Symbol %in% c(db.sel$`Receptor 1`, db.sel$`Receptor 2`, db.sel$`Receptor 3`)), which(colnames(data) !=
                                                                                                                                 "Symbol")], 2, psych::geometric.mean))
    }else {rec = NA}
  }
  mat=matrix(as.numeric(lig), ncol=1) %*% matrix(as.numeric(rec), nrow = 1)

  if (all(is.na(mat))){
    # mean that at least one ligand/receptor genes is not present in the data matrix.
    mat = matrix(nrow=dim(data)[2]-1, ncol=dim(data)[2]-1)
    colnames(mat)=colnames(dplyr::select_if(data, is.numeric))
    rownames(mat)=colnames(dplyr::select_if(data, is.numeric))
    print("Cannot display graph : at least one ligand/receptor genes is not present in the data matrix" )
    return(mat)
  }else{

    colnames(mat)=colnames(dplyr::select_if(data, is.numeric))
    rownames(mat)=colnames(dplyr::select_if(data, is.numeric))
    if (plot==T){
      melted=reshape2::melt(mat)
      colnames(melted)=c("Ligand", "Receptor", "value")

      p1 <-ggplot(melted, ggplot2::aes(x = Receptor,
                                       y = Ligand, fill=value)) +  geom_tile(aes(fill = value)) +
        scale_fill_gradient(low="white", high="blue")+
        theme(axis.text.x=ggplot2::element_text(angle=45, hjust = 1)) +
        ggplot2::labs(title = couple)
      return (p1)
    }else
      return(mat)
  }
}
