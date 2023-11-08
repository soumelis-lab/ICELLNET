
#' Convert CellPhoneDB database in a compatible format for ICELLNET analysis
#'
#' @description Provide CellPhoneDB ligand-receptor database in a format suitable for ICELLNET analysis
#'
#' @param complex complex_input.csv file from Github (ventolab/cellphonedb-data)
#' @param interaction interaction_input.csv file from Github (ventolab/cellphonedb-data)
#' @param gene_info gene_input.csv file from Github (ventolab/cellphonedb-data)
#' @param ppi logical, to retain only ppi interactions. TRUE by default
#'
#' @export
#' @examples
#' \dontrun{CellPhoneDB_convert()}
#'

CellPhoneDB_convert <- function(complex=utils::read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv"),
                                                        sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""),
                                interaction=utils::read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv"),
                                                            sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""),
                                gene_info=utils::read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv"),
                                                          sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""),
                                ppi=T){

  gene_info=gene_info[!duplicated(gene_info$uniprot),] # remove duplicated rows for hgnc symbol and uniprot (only ensembl differs in these cases)
  if (ppi==T){
    interaction= interaction %>% filter(annotation_strategy=="curated" & is_ppi == "True")
  }else{
    interaction= interaction %>% filter(annotation_strategy=="curated")
  }
  interaction$id_cp_interaction=paste0(interaction$partner_a, "_", interaction$partner_b)

  col=c("Ligand 1","Ligand 2","Ligand 3","Ligand 4","Receptor 1","Receptor 2","Receptor 3","Receptor 4", "Receptor 5",
        "Alias","Family", "Subfamily", "Other family", "Reference")

  db.new= as.data.frame(matrix(ncol=length(col) , nrow=length(interaction$id_cp_interaction)))
  colnames(db.new)=col
  db.new$Interaction_name=interaction$id_cp_interaction
  db.new$Reference=interaction$source # source

  non_ppi_interaction=filter(interaction, is_ppi=="False") %>% pull(id_cp_interaction)
  db.new$Family[which(db.new$Interaction_name %in% non_ppi_interaction)]="Not-protein ligands"


  # Transform cellphone db to be compatible with ICELLNET DB
  q=c() # q <- store the vector of weird interaction that do not fit with ICELLNET format. -> to remove at the end of the process
  for (mol in 1:length(db.new$Interaction_name)){
    #ligand
    if(interaction$partner_a[mol] %in% complex$complex_name){
      complex_info=complex[which(complex$complex_name==interaction$partner_a[mol]),]
      db.new$`Ligand 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_1)]
      if (!is.na(complex_info$uniprot_2)){
        db.new$`Ligand 2`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_2)]
      }
      if (!is.na(complex_info$uniprot_3)){
        db.new$`Ligand 3`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_3)]
      }
      if (!is.na(complex_info$uniprot_4)){
        db.new$`Ligand 4`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_4)]
      }
      if (!is.na(complex_info$uniprot_5)){
        q=c(q, mol)
      }
    }else{ db.new$`Ligand 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==interaction$partner_a[mol])] }

    #receptor
    if(interaction$partner_b[mol] %in% complex$complex_name){
      complex_info=complex[which(complex$complex_name==interaction$partner_b[mol]),]
      db.new$`Receptor 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_1)]
      if (!is.na(complex_info$uniprot_2)){
        db.new$`Receptor 2`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_2)]
      }
      if (!is.na(complex_info$uniprot_3)){
        db.new$`Receptor 3`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_3)]
      }
      if (!is.na(complex_info$uniprot_4)){
        db.new$`Receptor 4`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_4)]
      }
      if (!is.na(complex_info$uniprot_5)){
        db.new$`Receptor 5`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_5)]
      }
    }else{ db.new$`Receptor 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==interaction$partner_b[mol])]}
  }
  if (length(q)>=1){
    warning("The following rows were not included:")
    print(q)
  }
  db.new$Interaction_name=name.lr.couple(db=db.new)[,1]
  #remove duplicates
  db.new=db.new[!duplicated(db.new$Interaction_name),]
  return(db.new)
}

