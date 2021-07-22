
#  Creation of Lupus_Seurat_SingleCell_Landscape.Rds file for tutorial Exemple2_scRNAseq"
*authors: "Melissa Saichi and Lucile Massenet-Regad"*


## PART 1: Data acquisition ##

According to the Nat.Immunol publication, the dataset is publicly available on Immport/ SingleCell Broad Portal, under the accession code: SDY997
Important note: This accession code is for all the AMP consortium : it has other datasets on autoimmune diseases such as Rhumatoid Arthritis (RA) and SLE (equivalent to LN).

### Data acquisition steps:
+ Connect to SingleCell Portal by clicking on : https://portals.broadinstitute.org/single_cell/study/amp-phase-1
This will give you access to the "Download section".
+ The LN dataset was stored as : exprMatrixSleBroad.tsv.gz "Broad dataset" for the Immune part & exprMatrixSleMetro "Metro Dataset" for the "Stromal part".
+ Download the exprMatrixSleBroad.tsv.gz & the metadata file in which were stored the cell metadata of all studied cells from both datasets

## PART 2: Data analysis ##

```{r data loading}
rm(list=ls())
library(dplyr)
library(readr)
library(Seurat)
library(tibble)

meta <- readr::read_delim("~/Desktop/LN_melissa/meta", "\t",    escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()

ln <- readr::read_delim(gzfile("~/Desktop/LN_melissa/exprMatrixSleBroad.tsv.gz"), "\t", escape_double = FALSE, trim_ws = TRUE)
ln <- ln %>% column_to_rownames(var="gene")
stopifnot(identical(length(intersect(colnames(ln), meta$NAME)),ncol(ln)))

sln <- CreateSeuratObject(ln, project = "LN_PubData")

## keep only Metadata of interest:
meta <-  meta %>% dplyr::filter(NAME %in% colnames(sln)) %>% column_to_rownames( var="NAME")
stopifnot(identical(rownames(meta), colnames(sln)))
sln <- AddMetaData(sln, meta)

Idents(sln)=sln$Cluster #important to apply icellnet for desired clusters after. 

saveRDS(sln, "Lupus_Seurat_SingleCell_Landscape.Rds")

```



