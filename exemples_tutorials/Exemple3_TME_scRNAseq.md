
# Case study 3: dissect intercellular commmunication within clear cell renal carcinoma TME to identify cancer-cell specific interaction

Data used in this tutorial are coming from [Massenet-Regad et al., 2023](https://pubmed.ncbi.nlm.nih.gov/38025776/) study. Seurat object used here can be retrived from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222703)

### Load libraries and ICELLNET database and compute ligand-receptor pair names (db.name.couple function)

```{r,echo=T}
library(Seurat)
library(icellnet)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(gridExtra)

db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
db.name.couple=name.lr.couple(db, type="Family")
```
### 1 - Load Seurat object

```{r, warning=F, echo=T}

#Load data
seurat <- readRDS("GSE222703_SeuratObj_ALL_integrated.rds")
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
DimPlot(seurat, reduction = 'umap', group.by = 'Celltype_Harmony2', split.by="Tissue, label = T)
```
<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/TME_UMAP.png" width="1000" height="1000">

### 2 - Keep cells from the TME and preprocess data before ICELLNET score.
 
In this tutorial, we will focus only on intercellular communication between cells from tumoral samples (and not from juxtatumors). Once cells of interest are selected, `sc.data.cleaning()` function is used to: 
- remove all genes not present in icellnet database (or database filled with db argument)
- remove all the genes that are expressed by a very low number of cells (less than a percentage of cells per cluster, here 10%). This is chosen by the user with the filt.perc argument.If you are applying ICELLNET for the first time on your dataset, we advice to apply first ICELLNET without filtering (filter.perc = 0), and then with filtering at 5 or 10% to see the differences in filtered genes. This will help in the analysis and biological interpretation of the results. 

*Note: This step takes about 30 min with the current version of the database.*

```{r, warning=F, echo=T}
#select only cells from tumor samples
Idents(seurat)=seurat$Celltype_Harmony2
seurat.tum = subset(seurat, cells=which(seurat$Tissue=="Tumor"))

#data preprocessing
data <- sc.data.cleaning(object = seurat.tum, db = db, filter.perc = 10, save_file = T, force.file = F, path= "~/Desktop/")
data.icell= as.data.frame(gene.scaling(data, n=1, db=db))
```
### 3 - Apply icellnet pipeline on cluster of interest

In this exemple, we investigate tumor cell communication from ccRCC1 and ccRCC2 cluster towards the other cells in the TME. This means that we consider ligands expressed by the tumor cells and receptors expressed by the other cells from the TME.
Outward communication -> direction = "out"

**Define central (CC) and partner cells (PC)**

```{r, warning=F, echo=T}
CC1="ccRCC1"
CC2="ccRCC2"

# target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), #"Class"=colnames(data.icell))
# rownames(target)=target$ID

PC=c("B_cell","PlasmaC", "Treg","CD4T", "Prolif", "CD8T", "NK_T", "NK", "NK_CD160",
     "Mono_nc", "Mono_c", "Mac_C1QA","cDC2","cDC1", "pDC", 
     "Neutrop","Mast",
     "PT_GPX3", "PT_MT1G","ccRCC1","ccRCC2", "Epith", "Endoth","Fibro")

#direction communication score
direction="out"

# data selection 
CC1.data=as.data.frame(data.icell[,c(CC1, "Symbol")], row.names = rownames(data.icell))
CC2.data=as.data.frame(data.icell[,c(CC2, "Symbol")], row.names = rownames(data.icell))
PC.data=as.data.frame(data.icell[,c(PC, "Symbol")], row.names = rownames(data.icell))

```
**Compute intercellular communication scores**

```{r, warning=F, echo=T}
score.computation.1= icellnet.score(direction=direction, PC.data=PC.data, 
                                    CC.data= CC1.data, PC=PC,  db = db)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]

score.computation.2= icellnet.score(direction=direction, PC.data=PC.data, 
                                    CC.data= CC2.data, PC=PC,  db = db)
score2=as.data.frame(score.computation.2[[1]])
lr2=score.computation.2[[2]]

Scores=cbind(score1, score2)
colnames(Scores)=c(paste0(CC1,"_", direction), paste0(CC2,"_", direction))

````

**Visualization of contribution of family of molecules to communication scores**

```{r, warning=F, echo=T}
ymax=max(Scores)+1 #for later visualisations

# table of contribution of each family of molecule to the scores
LR.family.score(lr=lr1, db.couple=db.name.couple, plot=NULL) 

#display heatmap
LR.family.score(lr=lr1, db.couple=db.name.couple, plot="heatmap")

#display barplot
colors=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(length(unique(db$Family)))
LR.family.score(lr=lr1, db.couple=db.name.couple, plot="barplot", family.col=colors)  
```
<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/TME_heatmap_all.png" width="1200" height="1200">

**Remarks on biological interpretation:**

- **ICELLNET will always set, for each gene, maximum gene expression value at 10**. Then, the maximum score that you can obtain for an individual interaction is 100 (10 for the ligand, 10 for the receptor). 

- This means that high interaction scores does not mean high expression. You **should come back to the initial SeuratObject to look at individual gene expression**, and that the ligand/receptor of interest if effectively expressed by the cluster.

- **Filtering of genes expressed by each cluster according to cell percentage expressing the gene (= with counts >0) for each cluster can be an option to remove false-negative interactions scores.**  This can be done with the `sc.data.cleaning()` function, by setting filter.perc to a defined value (2 for 2%, 5 for 5% etc...). Filtered genes (expressed by a number of cells among the cluster below the percentage) will be set to 0. 


### 4 - Identify ccRCC2-specific communication channels within the TME
In this section, we focus on ccRCC2 cancer cells. The aim is to identify ligand-receptor specifically used by ccRCC2 to interact with the other cell types of the TME. 

We will : 
1) use `LR.viz()` to compute pairwise communication scores (outward communication between all possible cell pairs in the TM, i.e direction="out") of individual L/R interactions. The ouput of this step is a table gathering, for each cell type, its global contribution to each LR. 

2) filter the output to highlight the interactions specifically used by ccRCC2 (using `int.spe()` function provided below). This step will provide the list of interactions for which ccRCC2 communication scores are 1.5 times higher than for other cell types in the TME (thresh=1.5, can be tuned by the user).

```{r, warning=F, echo=T, eval=T}

#compute pairwise communication scores 
global_LR=matrix(nrow=length(rownames(db)), ncol=length(PC))
rownames(global_LR)=db.name.couple[,1]
colnames(global_LR)=colnames(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], int=as.character(db.name.couple[1,1]), db=db, plot=F))

direction="out"
for (i in db.name.couple[,1]){
  if (direction =="out"){
    global_LR[i,]=rowSums(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], int=i, db=db, plot=F)) # direction OUT
  }
  if (direction =="in"){
    global_LR[i,]=colSums(LR.viz(data = data.icell[,colnames(data.icell)%in%c(PC,"Symbol")], int=i, db=db, plot=F)) # direction IN
  }
}


#function to assess cell specificity of an interaction

int.spe <- function(LR.mat, cell, thresh){  # function to determine specific interactions for a cell type
  colMax= apply(LR.mat, 1, which.max)
  a=which(colnames(LR.mat)==cell)
  LR.mat2=LR.mat[which(colMax==a),]
  #create dataframe to store maximum value of "cell""
  value=data.frame("vmax"=LR.mat2[,a])
  rownames(value)=rownames(LR.mat2)
  value$vmax2=apply(LR.mat2[,-a], 1, function(x) (max(x)+0.1))
  value$vmax2_nb=apply(LR.mat2[,-a], 1, which.max)
  value$vmax2_cell=colnames(LR.mat2[,-a])[value$vmax2_nb]
  value=value%>% dplyr::select(-c("vmax2_nb"))
  value$vmean_pos=apply( LR.mat2[,-a], 1, function(x) mean(x[x!=0]))
  value$ratio=value$vmax/value$vmax2
  value = value %>% filter(ratio >thresh)
  return(value)
}

thresh=1.5

#give the list of interactions for which ccRCC2 contribution is 1.5 times higher than other cell types
specific_ccRCC2=int.spe(global_LR, cell="ccRCC2", thresh=thresh)

#plot the ccRCC2 specific interactions
LR.heatmap(lr = lr2[rownames(specific_ccRCC2),], thresh = 0 , topn= dim(specific_ccRCC2)[1], sort.by="var",  title="ccRCC2 specific interactions")
```

<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/TME_spe_ccRCC2_heatmap.png" width="600" height="900">


The plot above show all interactions specifically used by ccRCC2 to interact with the other cell types of the TME. To check the specificity of the interactions, LR.viz function can be used:

```{r, warning=F, echo=T}
spe1=LR.viz(data=data.icell, db = db, int="HHLA2 / TMIGD2", plot=T)
spe2=LR.viz(data=data.icell, db = db, int="CD70 / CD27", plot=T)
spe3=LR.viz(data=data.icell, db = db, int="VEGFA / NRP1", plot=T)
```

<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/TME_spe_single.png" width="1000" height="1000">



