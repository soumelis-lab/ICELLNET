
# EXEMPLE of on to run ICELLNET on bulk RNAseq cell-sorted populations. If you only have one population of interest, you can use the internal dataset provided with ICELLNET, 
# see https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple1_CAF.md

# Libraries
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)

# Load and select database
db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
db.name.couple=name.lr.couple(db, type="Family")
head(db.name.couple)

# Load bulk RNAseq data
data=as.data.frame(read.csv("~/Desktop/matrix_gene_expression.csv", sep=",", header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")) 
rownames(data)=data$Symbol # If gene name not already set up as rownames
#make sure colnames(data) give name of replicates

# Load target file with metadata info.
target=as.data.frame(read.csv("~/Desktop/target.csv", sep=",", header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")) 

# Data scaling
data.scaled=gene.scaling(data = data, n=1, db = db)

# data selection
CC.data.selection=target$ID[which(target$Cell_type=="cell_central")] #should give a vector with name of replicates to consider for central cell
PC.data.selection=target$ID[which(target$Cell_type%in% c("cell_1", "cell_2"))] #should give a vector with name of replicates to consider for partner cell                                  
                                   
my_Central_Cell_data=data.scaled[, CC.data.selection]
my_Partner_Cell_data = data.scaled[, PC.data.selection]

#compute communication score
score.computation.1= icellnet.score(direction="out", PC.data=my_Partner_Cell_data, CC.data= my_Central_Cell_data,  PC.target = target, PC=c("cell_1", "cell_2"),  CC.type = "RNAseq",  PC.type = "RNAseq",  db = db) 
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]] 
