
#' Comparison of different communication score by computing a pvalue matrix
#'
#' @description Return a pvalue matrix to allow communication score comparison.
#'
#'
#' @details Two types of pvalue can be computed, to compare either the communication scores obtained from the same central cell to different peripheral cells (between="cells"),
#' or to compare communication scores obtained from two different central cells corresponding to different biological conditions with the same peripheral cell (between="conditions").
#' If between="cells", the communication score is computed considering the average expression of ligands for the central cell, and each replicates separately
#' for the receptor expression of the peripheral cells. In this way, for one peripheral cell type, we obtain a distribution of n communication scores, n beeing the number
#' of peripheral cells replicates for this particular cell type. If between="conditions", then, the communication score is computed considering each replicates of the central
#' cell separately, and the average gene expression for the peripheral cells. We obtain a distribution of n communication scores, n beeing the number
#' of central cell replicates in one biological condition.
#' Then, a Wilcoxon statistical test is performed to compare the communication scores distributions. The pvalues are ajusted with p.adjust(), with "BH" method as a default.
#' It returns the pvalue matrix of statistical tests, that can be visualize in a heatmap with the pvalue.plot function.
#'
#' @param CC.data Data.frame of transcriptomic profiles corresponding to the central cell
#' @param CC.data2 Data.frame of transcriptomic profiles corresponding to the central cell in a different biological condition. NULL as a default.
#' @param PC.data Data.frame of transcriptomic profiles corresponding to the peripheral cells
#' @param PC Vector selecting a list of peripheral cell for the cell-cell communication analysis
#' @param PC.target Data.frame of information concerning the PC.data samples
#' @param CC.type Type of transcriptomic data for the central cell (either "RNAseq" or "Microarray")
#' @param PC.type Type of transcriptomic data for the peripheral (either "RNAseq" or "Microarray")
#' @param direction Direction of the communication (either "out" or "in"). "out" represent the communication from the central cell to the peripheral cells, which means that the ligands considered are expressed by the central cell, and the receptors are expressed by the peripheral cells. "in" represents the communication from the peripheral cells to the central cells, which means that the considered ligands are expressed by the peripheral cells and the receptors are expressed by the central cells.
#' @param db Database of ligand/receptors interactions
#' @param between Type of pvalue matrix desired : "cells" to compare communication scores of the central cell with different peripheral cells, "conditions" to compare the communication score of two central cells in different biological conditions with the same peripheral cell.
#' @param method method to adjust the pvalue. See p.adjust() function for the available methods. "BH" method is set as a default.
#'
#'
#' @export
#' @examples
#' \dontrun{
#' icellnet.score.pvalue( CC.data= CC.data.selection.S1, PC.data=PC.data, PC = my.selection,...
#'  PC.target = PC.target, CC.type = "RNAseq", PC.type = "Microarray", direction="out", between="cells")
#'  }
icellnet.score.pvalue = function(direction = c("out", "in"),
                                 CC.data = CC.data,
                                 CC.data2=NULL,
                                 PC.data = PC.data,
                                 PC = PC,
                                 PC.target = PC.target,
                                 CC.type = c("RNAseq", "Microarray"),
                                 PC.type = c("RNAseq", "Microarray"),
                                 db = db,
                                 between=c("cells","conditions"),
                                 method = "BH") {
  scores_rep=list()
  scores_rep1=list()
  scores_rep2=list()

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

  if (between=="cells"){
    pvalue = matrix(nrow = length(PC), ncol = length(PC))
    colnames(pvalue) = PC
    rownames(pvalue) = PC

    for (cell in PC){
      calc=icellnet.ind.score(direction = direction,
                              CC.data = CC.data,
                              PC.data = PC.data,
                              PC.target = PC.target,
                              cell=cell,
                              CC.type =CC.type ,
                              PC.type =PC.type,
                              db = db,
                              between=between)
      scores_rep=rlist::list.append(scores_rep, list(cell,calc))
    }
    #Compute pvalue from the score_rep matrix
    for (i in 1:(length(PC)-1)) {
      for (j in (i+1):length(PC)) {
        score1 = scores_rep[[i]][[2]]
        score2 = scores_rep[[j]][[2]]
        pvalue[PC[j], PC[i]] = stats::wilcox.test(score1, score2)[[3]]

      }
    }
    p.adjust2 <-matrix(stats::p.adjust(pvalue, method = "BH"),
                       nrow = length(PC),
                       ncol = length(PC))
    rownames(p.adjust2) = PC
    colnames(p.adjust2) = PC

    main <- list(p.adjust2, scores_rep)

  } else if (between=="conditions"){
    pvalue = matrix(nrow = length(PC), ncol = 1)
    colnames(pvalue) = "pvalue"
    rownames(pvalue) = PC

    for (cell in PC){
      calc1=icellnet.ind.score(direction = direction,
                               CC.data = CC.data,
                               PC.data = PC.data,
                               PC.target = PC.target,
                               cell=cell,
                               CC.type =CC.type ,
                               PC.type =PC.type,
                               db = db, between=between)
      scores_rep1=rlist::list.append(scores_rep1, list(cell,calc1))

      calc2=icellnet.ind.score(direction = direction,
                               CC.data = CC.data2,
                               PC.data = PC.data,
                               PC.target = PC.target,
                               cell=cell,
                               CC.type =CC.type ,
                               PC.type =PC.type,
                               db = db, between=between)
      scores_rep2=rlist::list.append(scores_rep2, list(cell,calc2))
    }

    for (i in 1:length(PC)){
      score1 = scores_rep1[[i]][[2]]
      score2 = scores_rep2[[i]][[2]]
      pvalue[i] = stats::wilcox.test(score1, score2)[[3]]
    }
    p.adjust2 <-
      matrix(stats::p.adjust(pvalue, method = "BH"),
             nrow = length(PC),
             ncol = 1)
    rownames(p.adjust2) = PC
    colnames(p.adjust2) ="pvalue"
    main <- list(p.adjust2, scores_rep1, scores_rep2)
  }else
    stop ('Error : between must be specified ("cells" or "conditions")')
  return(main)
}
