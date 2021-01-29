
#' ICELLNET score computation for replicates, to compute the ovalue (not exported)
#'
#' @description Needed to compute the pvalue in the icellnet.score.pvalue function. Returns communication scores between a central cell (CC) and partner cell types (PC).
#'
#' @details See icellnet.score.pvalue function for more details.
#'
#' The function returns the communication score per replicates.
#'
#' @param CC.data Data.frame of transcriptomic profiles corresponding to the central cell
#' @param PC.data Data.frame of transcriptomic profiles corresponding to the partner cells
#' @param PC Vector selecting a list of partnercell for the cell-cell communication analysis
#' @param PC.target Data.frame of information concerning the PC.data samples
#' @param cell A cell type
#' @param CC.type Type of transcriptomic data for the central cell (either "RNAseq" or "Microarray")
#' @param PC.type Type of transcriptomic data for the partner (either "RNAseq" or "Microarray")
#' @param direction Direction of the communication (either "out" or "in")
#' @param db Database of ligand/receptors interactions
#' @param between Type of pvalue matrix desired : "cells" to compare communication scores of the central cell with different partner cells, "conditions" to compare the communication score of two central cells in different biological conditions with the same partner cell.
#' @examples
#' \dontrun{icellnet.ind.score ( direction = "in", CC.data = CC.data, PC.data = PC.data,...
#' PC = PC, PC.target = PC.target, cell=cell, CC.type = "RNAseq", PC.type="RNAseq" , db = db, ...
#' between="cells")}

icellnet.ind.score = function(direction = c("out", "in"),CC.data = CC.data,
                              PC.data = PC.data,
                              PC = PC,
                              PC.target = PC.target,
                              cell=cell,
                              CC.type = c("RNAseq", "Microarray"),
                              PC.type = c("RNAseq", "Microarray"),
                              db = db,
                              between=c("cells","conditions"))
  {

  if (PC.type == "Microarray") {
    PC.affy.probes = as.data.frame(PC.data[, c(1, 2)])
    PC.affy.probes$ID = rownames(PC.affy.probes)
    PC_Probes_to_symbol = db.hgu133plus2(db, PC.affy.probes)
  }
  if (CC.type == "Microarray") {
    CC.affy.probes = as.data.frame(CC.data[, c(1, 2)])
    CC_Probes_to_symbol = db.hgu133plus2(db, CC.affy.probes)
  }
  #PC data
  cell.IDs=PC.target$ID[grepl(cell, PC.target$Cell_type) | grepl(cell, PC.target$Class)]
  if (length(cell.IDs)==0){
    stop(paste0( "Cell type ", cell, " is not found in the PC.target file"))
  }
  PC.data2= as.data.frame(PC.data[,cell.IDs])

  if (between =="cells"){
    score_cell= matrix(nrow=1, ncol=length(PC.data2))
    rownames(score_cell)=cell

    #Compute the LR product for each replicate of PC.data
    for (i in 1:dim(PC.data2)[2]){
      if (direction =="out"){
        if (PC.type == "Microarray") {
          rc.PC = receptor.average.MA(
            db = PC_Probes_to_symbol, data =as.data.frame(PC.data2[,i], row.names = rownames(PC.data2)))
        } else if (PC.type == "RNAseq") {
          rc.PC = receptor.average.RNAseq(db = db,
                                          data = as.data.frame(PC.data2[,i], row.names = rownames(PC.data2)))
        } else
          stop('Error : PC.type must be displayed ("Microarray" or "RNAseq").')
        if (CC.type == "Microarray") {
          lg.CC = ligand.average.MA(db = CC_Probes_to_symbol,
                                    data = as.data.frame(CC.data,row.names = rownames(CC.data)))
        } else if (CC.type == "RNAseq") {
          lg.CC = ligand.average.RNAseq(db = db,
                                        data = as.data.frame(CC.data,row.names = rownames(CC.data)))
        } else
          stop('Error : CC.type must be displayed ("Microarray" or "RNAseq").')
        score_cell[i]= sum((lg.CC * rc.PC), na.rm=T)

      }else if (direction=="in"){
        if (PC.type == "Microarray") {
          lg.PC = ligand.average.MA(
            db = PC_Probes_to_symbol, data =as.data.frame(PC.data2[,i], row.names = rownames(PC.data2)))
        } else if (PC.type == "RNAseq") {
          lg.PC = ligand.average.RNAseq(db = db,
                                        data = as.data.frame(PC.data2[,i], row.names = rownames(PC.data2)))
        } else
          stop('Error : PC.type must be displayed ("Microarray" or "RNAseq").')

        if (CC.type == "Microarray") {
          rc.CC = receptor.average.MA(db = CC_Probes_to_symbol,
                                      data = as.data.frame(CC.data,row.names = rownames(CC.data)))
        } else if (CC.type == "RNAseq") {
          rc.CC = receptor.average.RNAseq(db = db,
                                          data = as.data.frame(CC.data,row.names = rownames(CC.data)))
        } else
          stop('Error : CC.type must be displayed ("Microarray" or "RNAseq").')
        score_cell[i]= sum((lg.PC * rc.CC) , na.rm=T)
      }else
        stop('Error : Direction of the communication ("in" or "out") must be specified ')

    }
  } else if (between =="conditions"){
    score_cell= matrix(nrow=1, ncol=length(CC.data))
    rownames(score_cell)=cell

    #Compute the LR product for each replicate of PC.data
    for (j in 1:dim(CC.data)[2]){
      if (direction =="out"){
        if (PC.type == "Microarray") {
          rc.PC = receptor.average.MA(
            db = PC_Probes_to_symbol, data =as.data.frame(PC.data2, row.names = rownames(PC.data2)))
        } else if (PC.type == "RNAseq") {
          rc.PC = receptor.average.RNAseq(db = db,
                                          data = as.data.frame(PC.data2, row.names = rownames(PC.data2)))
        } else
          stop('Error : PC.type must be displayed ("Microarray" or "RNAseq").')
        if (CC.type == "Microarray") {
          lg.CC = ligand.average.MA(db = CC_Probes_to_symbol,
                                    data = as.data.frame(CC.data[,j],row.names = rownames(CC.data)))
        } else if (CC.type == "RNAseq") {
          lg.CC = ligand.average.RNAseq(db = db,
                                        data = as.data.frame(CC.data[,j],row.names = rownames(CC.data)))
        } else
          stop('Error : CC.type must be displayed ("Microarray" or "RNAseq").')
        score_cell[j]= sum((lg.CC * rc.PC), na.rm=T)

      }else if (direction=="in"){
        if (PC.type == "Microarray") {
          lg.PC = ligand.average.MA(
            db = PC_Probes_to_symbol, data =as.data.frame(PC.data2, row.names = rownames(PC.data2)))
        } else if (PC.type == "RNAseq") {
          lg.PC = ligand.average.RNAseq(db = db,
                                        data = as.data.frame(PC.data2, row.names = rownames(PC.data2)))
        } else
          stop('Error : PC.type must be displayed ("Microarray" or "RNAseq").')

        if (CC.type == "Microarray") {
          rc.CC = receptor.average.MA(db = CC_Probes_to_symbol,
                                      data = as.data.frame(CC.data[,j],row.names = rownames(CC.data)))
        } else if (CC.type == "RNAseq") {
          rc.CC = receptor.average.RNAseq(db = db,
                                          data = as.data.frame(CC.data[,j],row.names = rownames(CC.data)))
        } else
          stop('Error : CC.type must be displayed ("Microarray" or "RNAseq").')
        score_cell[j]= sum((lg.PC * rc.CC) , na.rm=T)
      }else
        stop('Error : Direction of the communication ("in" or "out") must be specified ')

    }

  }else
    stop ('Error : between must be specified ("cells" or "conditions")')
  return(score_cell)
}
