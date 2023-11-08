
#' Communication score computation with ICELLNET method
#'
#' @description Returns global communication scores between a central cell (CC) and peripheral cell types (PC).
#'
#' @details From the CC and PC expression matrices (RNA-seq or Microarray), this function computes a communication score using the knowledge-based LR database.
#' It selects ligands and receptors from the database an expressed by the cells of interest. A selection of the PC can be done.
#' The direction of the communication is indicated using the "direction" argument:
#' - out == compute the communication score from the CC to PC
#' - in == compute the communication score from the PC to CC
#'
#' The function returns the global communication score and the individual communication score matrices.
#'
#' @param CC.data Data.frame of transcriptomic profiles corresponding to the central cell
#' @param PC.data Data.frame of transcriptomic profiles corresponding to the peripheral cells
#' @param PC Vector selecting a list of peripheral cell for the cell-cell communication analysis
#' @param PC.target data.frame of information concerning the PC.data samples. Set to NULL by default: use only colnames(PC.data) as meaningful level of information.
#' @param CC.type Type of transcriptomic data for the central cell (either "RNAseq" or "Microarray")
#' @param PC.type Type of transcriptomic data for the peripheral (either "RNAseq" or "Microarray")
#' @param direction Direction of the communication (either "out" or "in")
#' @param db Database of ligand/receptors interactions
#'
#' @export
#' @examples
#' \dontrun{icellnet.score( CC.data= CC.data.selection.S1, PC.data=PC.data, PC = my.selection, ...
#' PC.target = PC.target, CC.type = "RNAseq", PC.type = "Microarray", direction="out")}
#'
#'
#'
icellnet.score <- function (direction = c("out", "in"), CC.data = CC.data, PC.data = PC.data,
                         PC = PC, PC.target = NULL, CC.type = c("RNAseq", "Microarray"),
                         PC.type = c("RNAseq", "Microarray"), db = db) {

  if (dim(CC.data)[2] == 1 | dim(PC.data)[2] == 1) {
    warning("Check that PC.data and/or CC.data contains rownames. Ignore this if it is the case")
  }
  score_global = matrix(nrow = length(PC), ncol = 1)
  rownames(score_global) = PC
  colnames(score_global) = direction
  lr_global = matrix(nrow = dim(db)[1], ncol = length(PC))
  colnames(lr_global) = PC
  rownames(lr_global) = name.lr.couple(db, type = "Family")[, 1]

  #create PC.target file when not provided
  if (is.null(PC.target)){
    if (is.null(colnames(PC.data))){stop("Error:colnames(PC.data) is NULL and should be defined")}
    else{
      PC.target=data.frame("ID"=colnames(PC.data), "Cell_type"=colnames(PC.data), "Class"=colnames(PC.data))
      rownames(PC.target)=PC.target$ID
    }
  }

  #when provided, PC.target format
  if (is.null(PC.target$Class)){
    PC.target$Class="NaN"
    print("PC.target$Class was not defined in input file. Defined in the code as NaN column")}
  if (is.null(PC.target$Cell_type)){
    PC.target$Cell_type="NaN"
    print("PC.target$Cell_type was not defined in input file. Defined in the code as NaN column")}

  if (PC.type == "Microarray") {
    PC.affy.probes = as.data.frame(PC.data[, c(1, 2)])
    PC.affy.probes$ID = rownames(PC.affy.probes)
    PC_Probes_to_symbol = db.hgu133plus2(db, PC.affy.probes)
  }
  if (CC.type == "Microarray") {
    CC.affy.probes = as.data.frame(CC.data[, c(1, 2)])
    CC.affy.probes$ID = rownames(CC.affy.probes)
    CC_Probes_to_symbol = db.hgu133plus2(db, CC.affy.probes)
  }
  if (direction == "out") {
    for (cell in PC) {
      cell.IDs = PC.target$ID[PC.target$Cell_type==cell |  PC.target$Class==cell]
      if (length(cell.IDs) == 0) {
        warning(paste0("Cell type ", cell, " is not found in the PC.target file"))
        score_global[cell, ] = "NaN"
      }
      else {
        if (PC.type == "Microarray") {
          rc.PC = receptor.average(db = PC_Probes_to_symbol,
                                      data = as.data.frame(PC.data[, cell.IDs],
                                                           row.names = rownames(PC.data)))
        }
        else if (PC.type == "RNAseq") {
          rc.PC = receptor.average(db = db, data = as.data.frame(PC.data[, cell.IDs], row.names = rownames(PC.data)))
        }
        else {
          stop("Error: PC.type must be displayed (\"Microarray\" or \"RNAseq\").")
        }
        if (CC.type == "Microarray") {
          lg.CC = ligand.average(db = CC_Probes_to_symbol,
                                    data = CC.data)
        }
        else if (CC.type == "RNAseq") {
          lg.CC = ligand.average(db = db, data = CC.data)
        }
        else {
          stop("Error: CC.type must be displayed (\"Microarray\" or \"RNAseq\").")
        }
        lr_global[, cell] = lg.CC * rc.PC
        score_global[cell, ] = sum(lg.CC * rc.PC, na.rm = T)
      }
    }
  }
  else if (direction == "in") {
    for (cell in PC) {
      cell.IDs = PC.target$ID[PC.target$Cell_type==cell |  PC.target$Class==cell]
      if (length(cell.IDs) == 0) {
        warning(paste0("Cell type ", cell, " is not found in the PC.target file"))
        score_global[cell, ] = "NaN"
      }
      else {
        if (PC.type == "Microarray") {
          lg.PC = ligand.average(db = PC_Probes_to_symbol,
                                    data = as.data.frame(PC.data[, cell.IDs],
                                                         row.names = rownames(PC.data)))
        }
        else if (PC.type == "RNAseq") {
          lg.PC = ligand.average(db = db, data = as.data.frame(PC.data[, cell.IDs], row.names = rownames(PC.data)))
        }
        else {
          stop("Error: PC.type must be displayed (\"Microarray\" or \"RNAseq\").")
        }
        if (CC.type == "Microarray") {
          rc.CC = receptor.average(db = CC_Probes_to_symbol,
                                      data = CC.data)
        }
        else if (CC.type == "RNAseq") {
          rc.CC = receptor.average(db = db, data = CC.data)
        }
        else {
          stop("Error: CC.type must be displayed (\"Microarray\" or \"RNAseq\").")
        }
        lr_global[, cell] = lg.PC * rc.CC
        score_global[cell, ] = sum(lg.PC * rc.CC, na.rm = T)
      }
    }
  }
  else {
    stop("Error: Direction of the communication (\"in\" or \"out\") must be specified ")
  }
  main <- list(score_global, lr_global)
  return(main)
}
