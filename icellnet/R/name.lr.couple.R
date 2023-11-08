
#' name.lr.couple function
#'
#' @description From the database of ligand/receptors interactions, the function retrieves the corresponding gene symbols.
#' It returns a table containing the name of the couple L/R and the family/subfamily of molecules which the interaction belongs to.
#'
#' @param db Ligand/receptor database
#' @param type Can be set to "Family" or "Subfamily", to display the  family/subfamily of molecules which the L/R interaction belongs to. If you want to focus on a specific subfamily of molecules from the icellnet database, set it as db$Subfamily. (ex: db$Subfamily = db$Checkpoints or db$Subfamily = db$Cytokines"
#' @export
#' @examples
#' \dontrun{name.lr.couple(db = db)}

name.lr.couple <- function (db = db, type ="Family") {
  db.name.lr.couple = matrix(nrow = dim(db)[1], ncol = 2)
  colnames(db.name.lr.couple) = c("Pair", type)
  if (is.null(db$`Ligand 1`) | is.null(db$`Ligand 2`) | is.null(db$`Ligand 3`) |
      is.null(db$`Receptor 1`) | is.null(db$`Receptor 2`) | is.null(db$`Receptor 3`) | is.null(db$`Receptor 4`) | is.null(db$`Receptor 5`) ) {
    warning("Check database columns names : database should contains Ligand 1, Ligand 2, Ligand 3, Receptor 1, Receptor 2, Receptor 3, Receptor 4, and Receptor 5 column names")
  }
  for (mol in seq(1, dim(db)[1])) {
    # ligand name
    if (!is.na(db$`Ligand 4`[mol])){
      ligand_name=paste(db$`Ligand 1`[mol], "+", db$`Ligand 2`[mol], "+", db$`Ligand 3`[mol], "+", db$`Ligand 4`[mol])
    }else{
      if (!is.na(db$`Ligand 3`[mol])){
        ligand_name=paste(db$`Ligand 1`[mol], "+", db$`Ligand 2`[mol], "+", db$`Ligand 3`[mol])
      }else{
        if (!is.na(db$`Ligand 2`[mol])){
          ligand_name=paste(db$`Ligand 1`[mol], "+", db$`Ligand 2`[mol])
        }else{ligand_name = db$`Ligand 1`[mol]}
      }
    }
    # receptor name
    if (!is.na(db$`Receptor 5`[mol])){
      receptor_name=paste(db$`Receptor 1`[mol], "+", db$`Receptor 2`[mol], "+", db$`Receptor 3`[mol], "+", db$`Receptor 4`[mol], "+", db$`Receptor 5`[mol])
    }else{
      if (!is.na(db$`Receptor 4`[mol])){
        receptor_name=paste(db$`Receptor 1`[mol], "+", db$`Receptor 2`[mol], "+", db$`Receptor 3`[mol], "+", db$`Receptor 4`[mol])
      }else{
        if (!is.na(db$`Receptor 3`[mol])){
          receptor_name=paste(db$`Receptor 1`[mol], "+", db$`Receptor 2`[mol], "+", db$`Receptor 3`[mol])
        }else{
          if (!is.na(db$`Receptor 2`[mol])){
            receptor_name=paste(db$`Receptor 1`[mol], "+", db$`Receptor 2`[mol])
          }else{receptor_name = db$`Receptor 1`[mol]}
        }
      }
    }
    # merge both
    db.name.lr.couple[mol] = paste(ligand_name, "/", receptor_name)
  }

  if (is.null(db$Family) | is.null(db$Subfamily)) {
    note("database should contains Family and Subfamily columns names if you want to classify interactions. To focus on a specific subfamily of molecules from the icellnet database, set it to db$Subfamily. (ex: db$Subfamily = db$Checkpoints or db$Subfamily = db$Cytokines")
  }
  if (type == "Family") {
    db.name.lr.couple[, 2] = db$Family
  } else if (type == "Subfamily") {
    db.name.lr.couple[, 2] = db$Subfamily
  } else stop("The type of classification (Family/Subfamily) must be specified")
  return(db.name.lr.couple)
}
