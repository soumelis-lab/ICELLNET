
#' name.lr.couple function
#'
#' @description From the database of ligand/receptors interactions, the function retrieves the corresponding gene symbols.
#' It returns a table containing the name of the couple L/R and the family/subfamily of molecules which the interaction belongs to.
#'
#' @param db Ligand/receptor database
#' @param type Can be set to "Family" or "Subfamily", to display the  family/subfamily of molecules which the L/R interaction belongs to. If you want to focus on a specific subfamily of molecules from the icellnet database, set it as db$Subfamily. (ex: db$Subfamily = db$Checkpoints or db$Subfamily = db$Cytokines"
#' @export
#' @examples
#' \dontrun{name.lr.couple(db = db, type="Family")}

name.lr.couple <-
  function(db = db,
           type = c("Family", "Subfamily")) {
    db.name.lr.couple = matrix(nrow = dim(db)[1], ncol = 2)
    colnames(db.name.lr.couple) = c("Pair", type)

    #Check input format
    if (is.null(db$`Ligand 1`) | is.null(db$`Ligand 2`) | is.null(db$`Receptor 1`) | is.null(db$`Receptor 2`)  | is.null(db$`Receptor 3`)){
      warning("Check database columns names : database should contains Ligand 1, Ligand 2, Receptor 1, Receptor 2 and Receptor 3 column names")
    }
    # Paste interactions pairs
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
    if (is.null(db$`Family`) | is.null(db$`Subfamily`) ){
      note("database should contains Family and Subfamily columns names if you want to classigy interaction. If you want to focus on a specific subfamily of molecules from the icellnet database, set it as db$Subfamily. (ex: db$Subfamily = db$Checkpoints or db$Subfamily = db$Cytokines")
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
