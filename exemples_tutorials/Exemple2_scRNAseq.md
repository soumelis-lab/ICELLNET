
# Case study 2: dissect intercellular commmunication in single cell data from a Seurat object

Data used in this tutorial are coming from [Arazi et al.2019](https://pubmed.ncbi.nlm.nih.gov/31209404/) study. 
The code used to create the SeuratObject used in the tutorial is provided [here](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Code_SeuratObject_Exemple2.md). 

### Load libraries and ICELLNET database and compute ligand-receptor pair names (db.name.couple function)
```{r,echo=T}

library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)

library(Seurat)

db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))

db.name.couple=name.lr.couple(db, type="Family")
head(db.name.couple)
```
### 1 - Load Seurat object

```{r, warning=F, echo=T}

#Load data
seurat <- readRDS(file = "Lupus_Seurat_SingleCell_Landscape.Rds")
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)

#only for UMAP visualization, not for ICELLNET purpose
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:50)
DimPlot(seurat, reduction = 'umap', group.by = 'Cluster', label = T)
```
<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/ICELLNET_scRNAseq_UMAP.png" width=50% height=50%>

### 2 - Retrieve gene expression matrix 
 
#### a - Compute manually average gene expression per cluster without filtering 

First, you need to group the cells according to your classification of interest (here, group cells by annotated cluster)
```{r, warning=F, echo=T}
# Taking into account the total nb of cells in each cluster
filter.perc=0
Idents(seurat)=seurat$Cluster
average.clean= sc.data.cleaning(object = seurat, db = db, filter.perc = filter.perc, save_file = T, path="path/", force.file = F)
```

#### b - Compute manually average gene expression per cluster with filtering for gene expression by a defined cell percentage at a cluster level

ICELLNET offers the possibility to filter the initial gene expression matrix to keep genes at least expressed by defined percentage of cell in their respective cluster (below 10%): 

```{r, warning=F, echo=T}
filter.perc=10
average.clean= sc.data.cleaning(object = seurat, db = db, filter.perc = filter.perc, save_file = T, path="path/", force.file = F)
```
This filtering allows to remove all the genes that are expressed by a very low number of cells in some clusters, to avoid false negative cell-cell interactions scores. If you are applying ICELLNET for the first time on your dataset, we advice to apply first ICELLNET without filtering, and then with filtering at 10% to see the differences and filtered genes. This will help in the analysis and biological interpretation of the results (see steps 3 and 4).

#### c - Compute manually average gene expression per cluster without starting from a SeuratObject.

If your are starting from a SeuratObject, continue directly on step 3.

If you are not starting from a SeuratObject but from a matrix of count, check that the matrix is normalized by library size and follow the steps below (explained below with Seurat formatting as reference for clarity).

```{r, warning=F, echo=T}
data <- as.data.frame(GetAssayData(seurat, slot = "data")) #or other matrix for which features expression values # are scaled by the total expression in each cell

target <- seurat@meta.data
target$Class=target$author_annotation
target$Cell=rownames(target)

average.manual=matrix(ncol=length(unique(target$author_annotation)), nrow=length(rownames(data)))
colnames(average.manual)=unique(target$author_annotation)
rownames(average.manual)=rownames(data)
dim(average.manual)
for (cell in unique(target$author_annotation)){
  cells.clust=target$Cell[which(target$author_annotation==cell)]
  average.manual[,cell]=apply(data[,which(colnames(data)%in%cells.clust)], 1, mean)
}

average.clean=average.manual
```

### 3 - Apply icellnet pipeline on cluster of interest

In this exemple, we investigate cDC to T cell communication from CM3 cluster (= conventional dendritic cells, 82 cells), to either CT3b or CT0a clusters (CT3b=TFH-like cells, 50 cells ; CT0a = effector memory CD4+ T cells, 220 cells). 


Format CC.data and PC.data and PC.target
```{r, warning=F, echo=T}

data.icell=as.data.frame(gene.scaling(as.data.frame(average.clean), n=1, db=db))

PC.data=as.data.frame(data.icell[,c("CT3b","CT0a", "Symbol")], row.names = rownames(data.icell))
my.selection=c("CT3b","CT0a")
```

**Compute intercellular communication scores**

We investigate conventional dendritic cells (cDCs, CM3 cluster) to T cell (either CT3b or CT0a clusters) outward communication, so this means that we consider ligands expressed by cDCs and receptors expressed by T cells to compute intercellular communication scores. 
Outward communication -> direction = "out"

```{r, warning=F, echo=T}
score.computation.1= icellnet.score(direction="out", PC.data=PC.data, 
                                    CC.data= as.data.frame(data.icell[,c("CM3")], row.names = rownames(data.icell)),  
                                    PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]

```
**Visualization of contribution of family of molecules to communication scores**

```{r, warning=F, echo=T}
# label and color label if you are working families of molecules already present in the database
my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding", "ECM")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  , "ECM"="slateblue4",
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "slateblue1",  "other" = "#908F90",  "NA"="#908F90")

ymax=round(max(score1))+1 #to define the y axis range of the barplot

LR.family.score(lr=lr1, db.couple=db.name.couple, plot=F) # table of contribution of each family of molecule to the scores

LR.family.score(lr=lr1, db.couple=db.name.couple, plot=T, title="DC-T comm") #display heatmap
```

<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/ICELLNET_scRNAseq_barplot.png" width=50% height=50%>

**Visualization of highest and most different interactions between the two conditions (selection of topn=20 interactions):** 
Can be displayed by `LR.heatmap()` and `LR.balloon.plot()` functions. 
`
**20 first most contributing interactions (sort.by="sum")**
```{r, warning=F, echo=T}
colnames(lr1)=c("CM3_to_CT3b", "CM3_to_CT0a")
LR.heatmap(lr = lr1, thresh = 0 , topn=20 , sort.by="sum",  title="Most contributing interactions")
```

**20 first most different interactions between the conditions (sort.by="var")**
```{r, warning=F, echo=T}
colnames(lr1)=c("CM3_to_CT3b", "CM3_to_CT0a")
LR.heatmap(lr = lr1, thresh = 0 , topn=20 , sort.by="var",  title="Most different interactions")
```
<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/ICELLNET_scRNAseq_heatmap_spe.png" width=50% height=50%>

### Remarks on biological interpretation: 

- **ICELLNET will always set, for each gene, maximum gene expression value at 10**. Then, the maximum score that you can obtain for an individual interaction is 100 (10 for the ligand, 10 for the receptor). 

- This means that high interaction scores does not mean high expression. You **should come back to the initial SeuratObject to look at individual gene expression**, and that the ligand/receptor of interest if effectively expressed by the cluster.

- **Filtering of genes expressed by each cluster according to cell percentage expressing the gene (= with counts >0) for each cluster can be an option to remove false-negative interactions scores.**  This can be done with the sc.data.clean function, by setting filter.perc to a defined value (2 for 2%, 5 for 5% etc...). Filtered genes (expressed by a number of cells among the cluster below the percentage) will be set to 0. 


# Broader cell-cell communication analysis and prioritisation strategies

Let's say we now want to decipher communication between conventional dendritic cells (cDCs, CM3 cluster) to all T cell clusters  (CT0a, CT0b, CT1, CT2, CT3a, CT3b, CT4, CT5a, CT5b, CT6 clusters) outward communication.

```{r, warning=F, echo=T}

PC=c("CT0a", "CT0b", "CT1", "CT2", "CT3a", "CT3b", "CT4", "CT5a", "CT5b", "CT6") #vectors including name of partner cells
PC.data=as.data.frame(data.icell[,c(PC,"Symbol")], row.names = rownames(data.icell))

PC.target=data.frame("ID"=colnames(data.icell), "Cell_type"=colnames(data.icell), "Class"=colnames(data.icell))
rownames(PC.target)=PC.target$ID

score.computation.1= icellnet.score(direction="out", PC.data=PC.data, 
                                    CC.data= as.data.frame(data.icell[,c("CM3")], row.names = rownames(data.icell)),  
                                    PC=PC, CC.type = "RNAseq", 
                                    PC.type = "RNAseq",  db = db)
score1=as.data.frame(score.computation.1[[1]])
lr1=score.computation.1[[2]]

```

## Heatmap representation 
Similar to `LR.balloon.plot`, the goal is to easily visualize specific interactions (most different between conditions or most contributing to the scores).

```{r, echo=T}
LR.heatmap(lr = lr1, thresh = 0 , topn=20 , sort.by="var",  title="Most different interactions")
LR.heatmap(lr = lr1, thresh = 0 , topn=20 , sort.by="sum",  title="Most contributing interactions")
```
![](exemples_tutorials/pictures/ICELLNET_scRNAseq_heatmap.png)

It is also possible to select the same interactions with the `LR.selection()` function and use other package for visualization such as *ComplexHeatmap package* that provides clustering of interaction and cell types. See exemple below.

```{r, echo=T}
pairs=LR.selection(lr = lr1, thresh = 0 , topn=20 , sort.by="var")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 100), c("white", "red"))
ComplexHeatmap::Heatmap(as.matrix(pairs), cluster_rows = T, cluster_columns = T, clustering_method_rows = "ward.D", name="Score", 
            clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names = T , show_row_names = T,   
            column_title ="Communication from cDC to T cell clusters", col = col_fun)
```

<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/ICELLNET_scRNAseq_CHeatmap.png" width=50% height=50%>

## Interaction specificity plot

This feature allows to represent graphically the communication score of a specific interaction for all combinations of cell pairs included in the dataset. *This representation is useful only when studying communication between cells all coming from the same dataset (ex: different clusters of a single cell RNAseq dataset* This representation allows to easily identify :
+ which clusters are expressing the ligand (in y axis) among all cell types included in the dataset
+ which clusters are expressing the receptor (x axis) among all cell types included in the dataset
+ the intensity of communication score for this specific interaction for all cell pairs
+ the specificity of the studied interaction

If plot=T, the function returns a graphical representation (ggplot object) of communication score for each cell pairs as a heatmap. If FALSE, the communication score matrix is returned.

A single interaction or a whole family of interactions can be displayed. Several interactions can be considered in the same plot when added as a vector

```{r, echo=T}
LR.viz(data=data.icell, db = db, int="CD86 / CD28", plot=T)
LR.viz(data=data.icell, db = db, int=c("CD86 / CD28", "TNFSF15 / TNFRSF25"), plot=T)
LR.viz(data=data.icell, db = db, family = "Checkpoint", plot=T)
```

![](exemples_tutorials/pictures/ICELLNET_scRNAseq_spe_3.png)






