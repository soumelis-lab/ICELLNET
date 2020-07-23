
#' name.lr.couple function
#'
#' @description From the database of ligand/receptors interactions, the function retrieves the corresponding gene symbols.
#' It returns a table containing the name of the couple L/R and the family/subfamily of molecules which the interaction belongs to.
#'
#' @param db Ligand/receptor database
#' @param type Can be set to "Family" or "Subfamily", to display the  family/subfamily of molecules which the L/R interaction belongs to.
#' @export
#' @examples
#' \dontrun{name.lr.couple(db = db, type="Family")}

name.lr.couple <-
  function(db = db,
           type = c("Family", "Subfamily")) {
    # rendre interne
    db.name.lr.couple = matrix(nrow = dim(db)[1], ncol = 2)
    colnames(db.name.lr.couple) = c("Pair", type)
    #c("Pair", "Family","Function")
    for (mol in seq(1, dim(db)[1])) {
      if (is.na(db$`Ligand 2`[mol])) {
        if (is.na(db$`Receptor 2`[mol]) & is.na(db$`Receptor 3`[mol])) {
          db.name.lr.couple[mol] = paste(db$`Ligand 1`[mol], "/", db$`Receptor 1`[mol])
        } else if (is.na(db$`Receptor 3`[mol])) {
          db.name.lr.couple[mol] = paste(db$`Ligand 1`[mol],
                                         "/",
                                         db$`Receptor 1`[mol],
                                         "+",
                                         db$`Receptor 2`[mol])
        } else{
          db.name.lr.couple[mol] = paste(
            db$`Ligand 1`[mol],
            "/",
            db$`Receptor 1`[mol],
            "+",
            db$`Receptor 2`[mol],
            "+",
            db$`Receptor 3`[mol]
          )
        }
      } else{
        if (is.na(db$`Receptor 2`[mol]) & is.na(db$`Receptor 3`[mol])) {
          db.name.lr.couple[mol] = paste(db$`Ligand 1`[mol],
                                         "+",
                                         db$`Ligand 2`[mol],
                                         "/",
                                         db$`Receptor 1`[mol])
        } else if (is.na(db$`Receptor 3`[mol])) {
          db.name.lr.couple[mol] = paste(
            db$`Ligand 1`[mol],
            "+",
            db$`Ligand 2`[mol],
            "/",
            db$`Receptor 1`[mol],
            "+",
            db$`Receptor 2`[mol]
          )
        } else{
          db.name.lr.couple[mol] = paste(
            db$`Ligand 1`[mol],
            "+",
            db$`Ligand 2`[mol],
            "/",
            db$`Receptor 1`[mol],
            "+",
            db$`Receptor 2`[mol],
            "+",
            db$`Receptor 3`[mol]
          )
        }
      }
    }
    if (type == "Family") {
      db.name.lr.couple[, 2] = db$Family
    } else
      if (type == "Subfamily") {
        db.name.lr.couple[, 2] = db$Subfamily
      } else
        stop ("The type of classification (Family/Subfamily) must be specified")
    return(db.name.lr.couple)
  }
