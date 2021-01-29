
# Case study 2: IL-10 controls an intercellular communication module in LPS activated dendritic cells <a name="Case-study-2:-IL-10-controls-an-intercellular-communication-module-in-LPS-activated-dendritic-cells"></a>

In this example, we are interested in **studying communication of resting and perturbed immune cells**. To explore the role of autocrine loops, we cultured LPS-activated human monocyte-derived dendritic cells (DCs) in the presence or absence of blocking antibodies (Abs) to the TNF and IL-10 receptors (αTNFR and αIL10R). We want to compare the communication channels that are used by the DCs in the different activation modes.

**Quick experimental information:** Primary cells were extracted from human blood, and the DCs were isolated by negative selection. They were then activated for 8 hours by being cultured either in presence of LPS, LPS+aTNFR, LPS+aIL10. The control condition corresponds to dendritic cells cultured with medium only. Transcriptomic profiles of dendritic cells in each condition were generated using **Affymetrix technology** (hgu133plus2 platform).

### Load database

```{r,echo=T}
db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/database.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
db.name.couple=name.lr.couple(db, type="Family")
head(db.name.couple)
```

### Load partner cell types: Human Primary Cell Atlas dataset 

```{r,echo=T}
my.selection = c("Epith", "Fblast","Endoth","Mono","Macroph","NK","Neutrop","B cell")
PC.target = PC.target.all[which(PC.target.all$Class%in%my.selection),c("ID","Class","Cell_type")]
PC.data = PC.data.all[,PC.target$ID]
```

To use Human Primary Cell Atlas dataset, we have to : 

1. Create a conversion chart between the AffyID and the gene symbol that are used in the database, using the hgu133plus2.db() function. 

2.  Perform the gene.scaling() function, that will a) select genes corresponding to the ligands and/or receptors included in the database (db). b) scale each ligand/receptor gene expression among all the conditions ranging from 0 to 10. For each gene: - the maximum value (10) is defined as the mean expression of the 'n' highest values of expression. - the minimum value (0) is defined as the mean expression of the 'n' lowest values of expression. Default value of n is 1. Outliers are rescaled at either 0 (if below minimum value) or 10 (if above maximum value).

In this example, n is set to take the mean of the 1% extreme expression values as the maximum/minimum.

```{r, warning=F ,echo=T}
### Convert the gene symbol to affy ID 
PC.affy.probes = as.data.frame(PC.data[,c(1,2)])
PC.affy.probes$ID = rownames(PC.affy.probes) # for format purpose
transform = db.hgu133plus2(db,PC.affy.probes) # creation of a new db2 database with AffyID instead of gene symbol

##Gene scaling of the partner celldataset
PC.data=gene.scaling(data = PC.data[, 2: dim(PC.data)[2]], n=0.01*dim(PC.data)[2], db = transform) 
PC.data$ID=rownames(PC.data) # for format purpose
PC.data$Symbol=rownames(PC.data) # for format purpose
```


### Load central cell: dendritic cell transcriptomic profiles
```{r,echo=T}
# Central cell data file (processed gene expression matrix)
data=read.table("data_DC.txt", sep="", header = T)
CC.data=data[,-dim(data)[2]]

#Target central cell file (description of the different samples)
CC.target = as.data.frame(read.table("target_DC.txt",sep = "\t",header=T))
```

Same as for the PC.data, the gene expression matrix is rescaled ranging from 0 to 10 considering all the CC.data samples. Here, the microarray data for the central cell are already annotated with gene symbol so we can consider them as "RNAseq" data for the next steps. 

```{r,echo=T}
CC.data= as.data.frame(gene.scaling(data = CC.data[,2:dim(CC.data)[2]], n=round(dim(CC.data)[2]*0.05), db = db)) #to not consider the SYMBOL column
CC.data$Symbol=rownames(CC.data) #for format purpose
```

### Selection of the different biological conditions for the central cell
```{r,echo=T}
CC.selection.S1 = CC.target[which(CC.target$Condition=="M+IgG_8h"),"Nomenclature"]
CC.selection.S2 = CC.target[which(CC.target$Condition=="L+IgG_8h"),"Nomenclature"]
CC.selection.S3 = CC.target[which(CC.target$Condition=="L+a-T_8h"),"Nomenclature"]
CC.selection.S4 = CC.target[which(CC.target$Condition=="L+ a-10_8h"),"Nomenclature"]

CC.data.selection.S1 = CC.data[,which(colnames(CC.data)%in%CC.selection.S1)]
CC.data.selection.S2 = CC.data[,which(colnames(CC.data)%in%CC.selection.S2)]
CC.data.selection.S3 = CC.data[,which(colnames(CC.data)%in%CC.selection.S3)]
CC.data.selection.S4 = CC.data[,which(colnames(CC.data)%in%CC.selection.S4)]
```

### Computation of ICELLNET intercellular communication scores

```{r,echo=T}
score.computation.1 = icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S1,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db, family.type = "Family")
score1=as.data.frame(score.computation.1[[1]])
colnames(score1)="M+IgG_8h"
lr1=score.computation.1[[2]]

score.computation.2 = icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S2,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db, family.type = "Family")
score2=as.data.frame(score.computation.2[[1]])
colnames(score2)="L+IgG_8h"
lr2=score.computation.2[[2]]

score.computation.3 = icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S3,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db, family.type = "Family")
score3=as.data.frame(score.computation.3[[1]])
colnames(score3)="L+a-TNFR"
lr3=score.computation.3[[2]]

score.computation.4 = icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S4,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db, family.type = "Family")
score4=as.data.frame(score.computation.4[[1]])
colnames(score4)="L+ a-IL10"
lr4=score.computation.4[[2]]

Scores= cbind(score1,score2,score3,score4)
colnames(Scores)=c("M+IgG","L+IgG","L+a-TNFR","L+ a-IL10")
```

score1, score2, score3 and score4 correspond to global scores, that are just the sum of the individual scores. Matrix of scores (Scores) corresponds to a summary of the global communication scores computed with ICELLNET between all peripheral cells and the central cell. lr1,lr2, lr3, lr4 correspond to individual score matrix, and will be useful for the further visualisation steps.

### Normalisation and rescaling of the global score 

In this example, we want to see the global variation of communication compared to the medium condition, so we are going to divide each communication score of perturbed condition (activated DCs) by the communication score in the medium condition. If you want to study the difference of communication score between different cell types and the central cell, you do not want to normalise (See case study 2). 
The Scores.norm matrix is then rescaled ranging from 0 to 10 to facilitate the visualisation of the intercellular network after.

```{r, echo=T, warning=FALSE, fig.align='center'}
#Score normalisation by the medium condition
Scores.norm=Scores
for (i in 1:length(my.selection)){
  Scores.norm[i,]=Scores[i,]/Scores[i,1]
}
Scores.norm
#Score scaling
Scores.norm2=(Scores.norm-min(Scores.norm))/(max(Scores.norm)-min(Scores.norm))*9+1
```


### Intercellular communication network representation

```{r, echo=T, warning=FALSE, fig.height=10, fig.width=12}

# Color label
PC.col = c("Epith"="#C37B90", "Muscle_cell"="#c100b9","Fblast_B"="#88b04b", "Fblast"="#88b04b","Endoth"="#88b04b",
           "Mono"="#ff962c","Macroph"="#ff962c","moDC"="#ff962c","DC1"="#ff962c","DC2"="#ff962c","pDC"="#ff962c","NK"="#ff962c","Neutrop"="#ff962c",
           "CD4 T cell"="#5EA9C3","CD8 T cell"="#5EA9C3","Treg"="#5EA9C3","B cell"="#5EA9C3")

# Display intercellular communication networks
network.plot1 = network.create(icn.score = Scores.norm2[1], scale = c(round(min(Scores.norm2)),round(max(Scores.norm2))), direction = "out", PC.col)
network.plot2 = network.create(icn.score =Scores.norm2[2], scale = c(round(min(Scores.norm2)),round(max(Scores.norm2))), direction = "out",PC.col)
network.plot3 = network.create(icn.score = Scores.norm2[3], scale = c(round(min(Scores.norm2)),round(max(Scores.norm2))), direction = "out", PC.col)
network.plot4 = network.create(icn.score =Scores.norm2[4], scale = c(round(min(Scores.norm2)),round(max(Scores.norm2))), direction = "out",PC.col)
grid.arrange(network.plot1, network.plot2, network.plot3, network.plot4, nrow=2, ncol=2)

```
Here, we observe a general increase of communication by blocking the IL10 communication channel, which suggests that autocrine IL10 secretion controls the communication in LPS-activated DCs.
To assess the differences between scores in a quantitative manner, a statistical test can be performed (see "Compute pvalue to compare communication scores" section). 
To see the use of the other visualisation modes, see case study 2. 

