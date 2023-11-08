#'
#' Plot individual communication scores of a specific L/R interactions for all cell pairs.
#'
#' @description  Plot individual communication scores of a specific L/R interactions for all cell pairs. Useful to see specificity of an interaction between two cell types. Relevant only for single cell datasets, deciphering communication between clusters. Returns either a ggplot object (if plot=TRUE) or a matrix.
#'
#'
#' @param data A dataframe of transcriptomic profiles with gene names as rownames
#' @param int Vector of characters corresponding to  "L / R" interactions. If several subunits, refer to output of name.lr.couple() function for interaction names.
#' @param family Vector of characters corresponding to family names included in ICELLNET DB.
#' @param db Ligand/receptor database
#' @param plot logical FALSE/TRUE (TRUE by default). If TRUE, the function returns the graphical representation of communication score for each cell pairs as a heatmap. If FALSE, the communication score matrix is returned.
#'
#' @export
#' @examples
#' \dontrun{LR.viz(data=data.icell , couple=couple, db=db , plot=TRUE)}
#'
#'
LR.viz <- function(data = data, int = NULL, family=NULL, db = db, plot = TRUE)
{

  SYMBOL = rownames(data)
  data = as.data.frame(data)
  if (is.null(data$Symbol)) { data$Symbol = SYMBOL }
  rownames(db) = name.lr.couple(db = db, type = "Family")[,1]

  if (!is.null(int) & !is.null(family)){
    stop("Set int = NULL or family=NULL ")
  }
  if (!is.null(int) & length(int)==1){
    db.sel = db[which(rownames(db) == int), ]
    if (dim(db.sel)[1] == 0) {
      stop("Interaction not included the DB - correct interaction")
    }
    #ligand profile
    if (all(as.vector(na.omit(c(db.sel$`Ligand 1`, db.sel$`Ligand 2`, db.sel$`Ligand 3`, db.sel$`Ligand 4`))) %in% data$Symbol)) {
      lig = as.vector(apply(data[which(data$Symbol %in% as.vector(na.omit(c(db.sel$`Ligand 1`, db.sel$`Ligand 2`, db.sel$`Ligand 3`, db.sel$`Ligand 4`)))),
                                 which(colnames(data)!= "Symbol")], 2, psych::geometric.mean))
    }else{
      lig = NA}

    #receptor profile
    if (all(as.vector(na.omit(c(db.sel$`Receptor 1`, db.sel$`Receptor 2`, db.sel$`Receptor 3`, db.sel$`Receptor 4`,  db.sel$`Receptor 5`))) %in% data$Symbol)) {
      rec = as.vector(apply(data[which(data$Symbol %in% as.vector(na.omit(c(db.sel$`Receptor 1`, db.sel$`Receptor 2`, db.sel$`Receptor 3`, db.sel$`Receptor 4`,  db.sel$`Receptor 5`))) ),
                                 which(colnames(data)!= "Symbol")], 2, psych::geometric.mean))
    }else{
      rec = NA
    }

    mat = matrix(as.numeric(lig), ncol = 1) %*% matrix(as.numeric(rec), nrow = 1)
    title=int

  }else if(!is.null(family) | length(int) >1){ #same but several interactions to consider
    if (is.null(family)){
      db.sel=db[unique(int),]
      title=paste0(length(unique(int)), " interactions")
    }else{
      db.sel = dplyr::filter(db, Family==family | Subfamily==family)
      title=family
    }

    if (dim(db.sel)[1] == 0) {
      stop("Family not included the DB - correct family name")
    }
    mat = matrix(nrow = dim(data)[2] - 1, ncol = dim(data)[2] - 1, 0)  # create an empty matrix
    for (mol in seq(1, dim(db.sel)[1])) {
      #ligand profile
      if (all(as.vector(na.omit(c(db.sel$`Ligand 1`[mol], db.sel$`Ligand 2`[mol], db.sel$`Ligand 3`[mol], db.sel$`Ligand 4`[mol]))) %in% data$Symbol)) {
        lig = as.vector(apply(data[which(data$Symbol %in% as.vector(na.omit(c(db.sel$`Ligand 1`[mol], db.sel$`Ligand 2`[mol], db.sel$`Ligand 3`[mol], db.sel$`Ligand 4`[mol])))),
                                   which(colnames(data)!= "Symbol")], 2, psych::geometric.mean))
      }else{
        lig = NA}

      #receptor profile
      if (all(as.vector(na.omit(c(db.sel$`Receptor 1`[mol], db.sel$`Receptor 2`[mol], db.sel$`Receptor 3`[mol], db.sel$`Receptor 4`[mol],  db.sel$`Receptor 5`[mol]))) %in% data$Symbol)) {
        rec = as.vector(apply(data[which(data$Symbol %in% as.vector(na.omit(c(db.sel$`Receptor 1`[mol], db.sel$`Receptor 2`[mol], db.sel$`Receptor 3`[mol], db.sel$`Receptor 4`[mol],  db.sel$`Receptor 5`[mol])))),
                                   which(colnames(data)!= "Symbol")], 2, psych::geometric.mean))
      }else{
        rec = NA
      }

      mat_int = matrix(as.numeric(lig), ncol = 1) %*% matrix(as.numeric(rec), nrow = 1)
      if (dim(mat_int)[1]>1 & dim(mat_int)[2]>1){
        mat = mat + mat_int
      }
    }

  }else{stop("Provide an interaction (int argument) or family of molecules (family argument)")}

  if (all(is.na(mat))) {
    mat = matrix(nrow = dim(data)[2] - 1, ncol = dim(data)[2] - 1, 0)  # create an empty matrix instead
    warning("Cannot display graph : at least one ligand/receptor genes is not present in the data matrix")
  }
  colnames(mat) = colnames(dplyr::select_if(data, is.numeric))
  rownames(mat) = colnames(dplyr::select_if(data, is.numeric))
  if (plot == T) {
    melted = reshape2::melt(mat)
    colnames(melted) = c("Ligand", "Receptor", "value")
    p1 <- ggplot(melted, ggplot2::aes(x = Receptor, y = Ligand,
                                      fill = value)) + geom_tile(aes(fill = value)) +
      ggplot2::scale_fill_gradient(low = "white", high = "blue") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                hjust = 1)) + ggplot2::labs(title = title)
    return(p1)
  }
  else return(mat)
}

