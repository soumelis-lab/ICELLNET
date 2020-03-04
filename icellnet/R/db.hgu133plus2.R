
#'
#' Conversion of affymetrix probe name into gene symbol
#'
#' @description This function converts LR database gene symbol annotation into Affymetrix Human Genome U133 Plus 2.0 Array annotation.
#'
#' @details This function transform the annotation of the database into probeset annotation, only available for Affymetrix Human Genome U133 Plus 2.0 Array.
#'
#'
#' @param data List of annotation data (ex: rownames of the transcriptomic profiles matrix)
#' @param db Ligand/receptor database
#' @export
#' @examples
#' \dontrun{
#' db.hgu133plus2(db,data.probes)
#' }
#'
#'
db.hgu133plus2 = function(db, data) {
  OUT <-
    AnnotationDbi::select(hgu133plus2.db, data$ID, c("SYMBOL", "ENTREZID", "GENENAME"))
  new.db = db
  #Ligand
  for (mol in seq(1, dim(db)[1])) {
    if (length(which(OUT$SYMBOL == db$`Ligand 1`[mol])) == 1) {
      new.db[mol, "Ligand 1"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Ligand 1`[mol])]
    } else if (length(which(OUT$SYMBOL == db$`Ligand 1`[mol])) > 1) {
      if (is.na(jmap(chip = "hgu133plus2", symbol = db$`Ligand 1`[mol])))
        new.db[mol, "Ligand 1"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Ligand 1`[mol])][2]
      else
        new.db[mol, "Ligand 1"] = jmap(chip = "hgu133plus2", symbol = db$`Ligand 1`[mol])
    }
    if (!is.na(db$`Ligand 2`[mol])) {
      if (length(which(OUT$SYMBOL == db$`Ligand 2`[mol])) == 1)
        new.db[mol, "Ligand 2"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Ligand 2`[mol])]
      else if (length(which(OUT$SYMBOL == db$`Ligand 2`[mol])) > 1) {
        if (is.na(jmap(chip = "hgu133plus2", symbol = db$`Ligand 2`[mol])))
          new.db[mol, "Ligand 2"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Ligand 2`[mol])][2]
        else
          new.db[mol, "Ligand 2"] = jmap(chip = "hgu133plus2", symbol = db$`Ligand 2`[mol])
      }
    }
  }
  #Receptor
  for (mol in seq(1, dim(db)[1])) {
    if (length(which(OUT$SYMBOL == db$`Receptor 1`[mol])) == 1) {
      new.db[mol, "Receptor 1"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 1`[mol])]
    } else if (length(which(OUT$SYMBOL == db$`Receptor 1`[mol])) > 1) {
      if (is.na(jmap(chip = "hgu133plus2", symbol = db$`Receptor 1`[mol])))
        new.db[mol, "Receptor 1"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 1`[mol])][2]
      else
        new.db[mol, "Receptor 1"] = jmap(chip = "hgu133plus2", symbol = db$`Receptor 1`[mol])
    }
    if (!is.na(db$`Receptor 2`[mol])) {
      if (length(which(OUT$SYMBOL == db$`Receptor 2`[mol])) == 1)
        new.db[mol, "Receptor 2"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 2`[mol])]
      else if (length(which(OUT$SYMBOL == db$`Receptor 2`[mol])) > 1) {
        if (is.na(jmap(chip = "hgu133plus2", symbol = db$`Receptor 2`[mol])))
          new.db[mol, "Receptor 2"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 2`[mol])][2]
        else
          new.db[mol, "Receptor 2"] = jmap(chip = "hgu133plus2", symbol = db$`Receptor 2`[mol])
      }
    }
    if (!is.na(db$`Receptor 3`[mol])) {
      if (length(which(OUT$SYMBOL == db$`Receptor 3`[mol])) == 1)
        new.db[mol, "Receptor 3"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 3`[mol])]
      else if (length(which(OUT$SYMBOL == db$`Receptor 3`[mol])) > 1) {
        if (is.na(jmap(chip = "hgu133plus2", symbol = db$`Receptor 3`[mol])))
          new.db[mol, "Receptor 3"] = OUT$PROBEID[which(OUT$SYMBOL == db$`Receptor 3`[mol])][2]
        else
          new.db[mol, "Receptor 3"] = jmap(chip = "hgu133plus2", symbol = db$`Receptor 3`[mol])
      }
    }
  }
  return(new.db)
}
