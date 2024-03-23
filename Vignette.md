This vignette explains the use of the ICELLNET package and demonstrates typical workflows to dissect intercellular communication between multiple cell types, based on transcriptomic profiles.

----
# Table of content

- [Introduction to ICELLNET R package](#Introduction-to-ICELLNET-R-package)
- [ICELLNET ligand/receptor interaction database](#ICELLNET-ligand/receptor-interaction-database) 
- [Input data](#Input-data)

- [Typical workflow](#Typical-workflow)
- [How is the intercellular communication score computed?](#How-is-the-intercellular-communication-score-computed?)
- [Visualization modes](#Visualization-modes)

  
- [How to install ICELLNET package?](#How-to-install-ICELLNET-package?)
- [How to format your own data to ue ICELLNET package?](#How-to-format-your-own-data-to-use-ICELLNET-package?)
- [Use cases exemples](#Use-cases-exemples)

- [Software information](#Software-information)

<!-- toc -->
----


# Introduction to ICELLNET R package  <a name="Introduction-to-ICELLNET-R-package"></a>

Cell-to-cell communication is at the basis of the higher-order organization observed in tissues and organisms, at steady state and in response to stress. The availability of large-scale transcriptomics datasets from several cell types has opened the possibility of **reconstructing cell-cell interactions based on co-expression of ligand-receptor pairs**.

We developed **ICELLNET**, a transcriptomic-based framework to **dissect cell communication in a global manner**. It integrates an original expert-curated **database of ligand-receptor interactions** taking into account multiple subunits expression. Based on transcriptomic profiles, ICELLNET package allows to compute **communication scores** between cells and provides **several visualization modes** that are helpful to dig into cell-cell interaction mechanism and extend biological knowledge. 

# ICELLNET ligand/receptor interaction database <a name="ICELLNET-ligand/receptor-interaction-database"></a>

We curated a comprehensive database of **1669 ligand-receptor interactions** from the literature and public databases. This database takes into account the **multiple subunits** of the ligands and the receptors. Interactions have been classified into 10 families of communication molecules: **Cytokine, Chemokine, Checkpoint, Notch pathway, HLA recognition, Innate immune, Wnt pathway, Growth factor, Cell adhesion and ECM interaction**. 

Cytokines have been further classified into 7 subfamilies according to reference classifications essentially based on structural protein motifs: type 1 cytokines, type 2 cytokines, IL-1 family, IL-17 family, TNF family, TGFb family and RTK cytokines. 

Other interactions and classifications of molecules will be implemented. The most recent version of ligand-receptor interaction database can always be downloaded [here](https://github.com/soumelis-lab/ICELLNET/blob/master/data/ICELLNETdb.tsv).
In R, you can visualize ICELLNET database and its structure: 

```{r db, echo=T}
db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
head(db)
```

You can use either all the database or restrict it by selecting some specific class of molecules (Cytokines, Growth factor etc..). Below, we show you how to restrict the study to cytokines, chemokines, and checkpoints, and how you can take into consideration subfamily of molecules.

```{r, echo=T}
#list of the different family of molecules considered in the database
table(db$Family)

#list of the different subfamily of cytokines considered in the database
table(db$Subfamily[which(db$Family=="Cytokine")])

#Restrict the database to some family of molecules 
my.selection.LR=c("Cytokine", "Chemokine", "Checkpoint")
db2 <- db[grepl(paste(my.selection.LR, collapse="|"),db$Family),] 
db.name.couple=name.lr.couple(db2, type="Family")
head(db.name.couple)

#Restrict the database to cytokines and consider the subfamilies of cytokines
my.selection.LR=c("Cytokine")
db3 <- db[grepl(paste(my.selection.LR, collapse="|"),db$Family),] #if you want to use all the database, do instead : db2=db
db.name.couple=name.lr.couple(db3, type="Subfamily")
head(db.name.couple)
```
Instead of using the ICELLNET database, it is also possible to use its own database as long as it is correctly formatted with specific columns as below. The Family and Subfamily colums correspond to two independant classifications (per family of molecules, or other) of your choice, but each interaction should fit only in one category of the classification (for exemple, an interaction cannot be classified in "type 1" and also "type 2" cytokines in the ICELLNET database).

|  Ligand 1 | Ligand 2  | Ligand 3  | Ligand 4  | Receptor 1  | Receptor 2  | Receptor 3  |  Receptor 4  |  Receptor 5  | Family | Subfamily | Other family |
|---|---|---|---|---|---|---|---|---|---|---|---|
|   |   |   |   |   |   |   |   |   |   |   |   |
|   |   |   |   |   |   |   |   |   |   |   |   |


# Input data <a name="Input-data"></a>

ICELLNET pipeline considers transcriptomic profiles of a defined "central cell", that can correspond for exemple to a cell type in several biological conditions. ICELLNET will then allow to compare the communication channels used by the central cells in these different conditions with partner cells.
As partner cells, we can use Human Primary Cell Atlas, a public datasets of 745 transcriptomic profiles among 31 cell types generated with the same technology (Affymetrix microarray, hgu133plus2 platform), already processed. 

As partner cells, the user can also use other transcriptomic profiles instead of Human Primary Cell Atlas (see section  [How to format your own data to use ICELLNET package?](#How-to-format-your-own-data-to-use-ICELLNET-package?) for format details).

```{r,echo=T}
#download PC.data.all and PC.target.all objects from the github and open them on your Rstudio session - adapt path if needed
PC.data.all=as.data.frame(read.csv("~/Downloads/PC.data.all.csv", sep=",", header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
PC.target.all=as.data.frame(read.csv("~/Downloads/PC.target.all.csv", sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))

head(PC.data.all[1:5,1:5])
```

It is possible to select different cell types to connect with the central cell. The different options are visible by running `table(PC.target.all$Class_broad)` and `table(PC.target.all$Class)`, thus offering different level of granularity. 

List of possible partner cell types :  Adipo, Adr_med, B cell, chondro, DC (or if desired: CD16 DC, DC1, DC2, moDC, pDC), Endoth, Epith, ESM, gameto, Hepato, HSC, Kerat, lps_cells, Macroph, Mono, Neutrop, NK, ostblast, Shwann, SMC, T cell (or if desired: CD4 T cell, CD8 T cell, Treg), TSC: 

```{r,echo=T}
table(PC.target.all$Class_broad)
table(PC.target.all$Class)
```

# Typical workflow <a name="Typical-workflow"></a>
 
Here we describe the different steps of the ICELLNET package to analyze cell-cell communication: 

1. Selection of the genes coding for the ligands and the receptors in our database from the transcriptomic profiles of the central cell and the partner cells. 

2. Rescale gene expression to avoid communication score to be driven only by highly expressed genes.

3. Compute ICELLNET communication scores (`direction="out"` for outward communication, `direction="in"` for inward communication, see [How is the intercellular communication score computed?](##How-is-the-intercellular-communication-score-computed?) for more details)

4. Display different visualization modes to dissect intercellular communication scores



# How is the intercellular communication score computed? <a name="How-is-the-intercellular-communication-score-computed?"></a>

The quantification of intercellular communication consist of scoring the intensity of each ligand/receptor interaction between two cell types with known expression profiles. No filtering threshold is applied on the L/R expression. If the communication molecule (ligand or receptor or both) is not expressed by a cell, the score will be zero. By default, all interactions of the database are considered to compute the score. It is also possible to reduce the number of interactions by manually selecting specific families of molecules in the database or considering DEG to compute the score, depending on the biological question. Whenever needed, we take into account multiple ligand units, or receptor chains, using logical rules.

The score of an individual ligand/receptor interaction is computed as the product of their expression levels respectively by the source (central) and by the target (partner) cell. These individual scores are then combined into a global metric assessing the overall exchange of information between the cell types of interest

Since cell-to-cell communication is directional, we consider ligand expression from the central cell and receptor expression from the partner cells to assess outward communication (`direction="out"`). On the other way, we select receptor expression from the central cell and ligand expression from partner cells to assess inward communication (`direction="in"`). This is controlled by the  *direction* argument ("in" or "out") in the `icellnet.score()` function. 

```{r,echo=T}
#not run - exemple
score.computation.1= icellnet.score(direction="out", PC.data=PC.data, CC.data= CC.data.selection.S1,  
                                    PC.target = PC.target, PC=my.selection, CC.type = "RNAseq", 
                                    PC.type = "Microarray",  db = db2)
score1=as.data.frame(score.computation.1[[1]]) #communication scores
lr1=score.computation.1[[2]] # detail of the ligand/receptor interactions scores matrix

```


# Visualization modes <a name="Visualization-modes"></a>

<img src="https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/pictures/ICELLNET_visualisations.png" width=100% height=100%>


### Intercellular communication network representation

This allows to visualize intercellular communication networks in a global manner through the function `network.create()`. In these directed graphs, nodes represent cell types, the width of the edges connecting two cell types is proportional to a global measure of the intensity of the communication between them and the arrows indicate the direction of communication. 


### Communication molecules distribution

`LR.family.score()` allows to dissect the global scores at a level of class of molecules, and allows to identify patterns of co-expressed molecules from the same family. This layer of analysis helps the interpretation on a qualitative level. 
The contribution of each family of molecules to the communication scores can be displayed as a heatmap (when setting plot="heatmap"), or as a barplot (when setting plot="barplot") 

### Individual communication scores distribution

The balloon plot (`LR.balloon.plot()` function) and heatmap (`LR.heatmap()` function)  are the deepest level of representation of the communication, displaying the most contributing ligand/receptor pairs to the communication score. This allows to identify specific individual interactions that can drive the intercellular communication and should be confirmed experimentally. 
*topn* parameter set the top n highest interactions to display (ex: `topn=10` means display the 10 highest contribution of specific interactions to the communication scores)
*thresh* parameter allows to set a threshold, in order to display interaction contributing to communication scores only above this value. (ex: thresh=20, only interactions with an individual score above 20 will be visualized). Looking at individual scores distribution is helpful to identify relevant threshold values for you case.
*sort.by* parameter ("sum" as default, can be set to "var") to ran the interactions either as the most contributing interactions in all conditions ("sum"), or the most different interactions among conditions ("var", ranking by computing the variance among scores for each interaction)


### Specificity of an interaction (for single cell dataset only)

This feature (`LR.viz()` function) allows to represent graphically the communication score of a specific interaction for all combinations of cell pairs included in the dataset. *This representation is useful only when studying communication between cells all coming from the same dataset (ex: different clusters of a single cell RNAseq dataset* This representation allows to easily identify :
+ which clusters are expressing the ligand (in y axis) among all cell types included in the dataset
+ which clusters are expressing the receptor (x axis) among all cell types included in the dataset
+ the intensity of communication score for this specific interaction for all cell pairs
+ the specificity of the studied interaction

If plot=T, the function returns a graphical representation (ggplot object) of communication score for each cell pairs as a heatmap. If FALSE, the communication score matrix is returned.

For usage of this function, see end of [Exemple_2](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple2_scRNAseq.md)


### Compute pvalue to compare communication scores

Two types of pvalue can be computed (`icellnet.score.pvalue()` function), to compare either the communication scores obtained from the same central cell to different partner cells (between="cells"), or to compare the communication scores obtained from two different central cells corresponding to different biological conditions with the same partner cell (between="conditions").If between="cells", the communication score is computed considering the average expression of ligands for the central cell, and each replicates separately for the receptor expression of the partner cells. In this way, for one partner cell, we obtain a distribution of n communication scores, n beeing the number of partner cells replicates for this particular cell type. If between="conditions", then, the communication score is computed considering each replicates of the central cell separately, and the average gene expression for the partner cells. We obtain a distribution of n communication scores, n beeing the number of central cell replicates in one biological condition. Then, a Wilcoxon statistical test is performed to compare the communication scores distributions. The pvalues are ajusted with `stats::p.adjust()`, with "BH" method as a default. 

It returns the pvalue matrix of statistical tests, that can be visualize as a heatmap with the `pvalue.plot()` function. This allows to interpret the difference of communication score in a quantitative manner.
 
# How to install ICELLNET package? <a name="How-to-install-ICELLNET-package?"></a>

To install `icellnet` package, the easiest way is to use the `R` package `devtools` and its function `devtools::install_github()`. If you don't have all the dependancies needed to use ICELLNET package, run the commands below:  

    install.packages(c("devtools", "curl", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist")) ##Installs devtools and the icellnet CRAN dependancies

    if (!requireNamespace("BiocManager", quietly = TRUE)) # Installs Bioconductor dependancies 
    install.packages("BiocManager")
    BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))
    
Then you just have to load `devtools` package and run the command below:

    library(devtools)
    install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
    
Then load the dependancies below and `icellnet` package.

    library(BiocGenerics)
    library("hgu133plus2.db")
    library("org.Hs.eg.db")
    library(jetset)
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
    library(icellnet)


To load the most recent version of the ligand/receptor database, you should run the command below: 

    db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))


# How to format your own data to use ICELLNET package? <a name="How-to-format-your-own-data-to-use-ICELLNET-package?"></a>

## Data files format
**For the central cell:** It can be any transcriptomic profile data of one cell type. For **RNA-seq data**, the dataset should be a datframe annotated with gene symbol as rownames. For **microarray data**, the ICELLNET functions are adapted to handle Affymetrix Human Genome U133 Plus 2.0 Array annotation. Nevertherless, if the dataset have been generated with an other Affymetrix technology, you have 2 possibilities to adapt the tool : a) Annotate your data with gene symbol before using ICELLNET and then consider your data as "RNA-Seq" for CC.type argument. b) adapt the R code of the `db.hgu133plus2()` function to have the adapted annotation conversion when using ICELLNET. Gene annotations should be set as rownames. 


**For the partner cell (if you don't want to use Human Primary Cell Atlas as partner cells):** This can be interesting for exemple if you possess transcriptomic data of several cell types of the same sample, to see how they interact together. As for the central cell, the transcriptomic profiles should be correctly formated (see previous paragraph above for more information). If your transcriptomic profiles are annotated with gene symbol, PC.type should be set to "RNA-seq" (even if your data come from microarray technology). 

## Target files format
You can define a target file. This dataframe usually describes the different samples with known metadata. If not provided, a target file will be automatically generated containing the sample names as the only metadata information available.

**PC.target** should contains at least an 'ID' column including the name of the samples (usually rownames(PC.target) or colnames(PC.data) ), and a 'Class' column corresponding to a classification of your different samples included in PC.target, such as a cell type classification. The different categories included in the 'Class' column will define the different partner cells in the graphs.

If not provided for `icellnet.score()` function, the column names of the expression matrix will be used as default to create a PC.target file.

# Use cases exemples <a name="Use-cases-exemples"></a>

- [Case study 1](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Exemple1_CAF.md): dissect intercellular commmunication of Cancer Associated Fibroblasts subsets. Show how to apply ICELLNET pipeline on transcriptomic profiles from 2 CAF-subsets, and how to restrict use of icellnet database to cytokines only (or other family of molecules).

- [Case study 2](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Exemple2_scRNAseq.md): Application of ICELLNET pipeline to scRNAseq from a Seurat object to infer intercellular communication between clusters.

- [Case study 3](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Exemple3_scRNAseq_TME.md): Application of ICELLNET to identify cancer cell-specific communication channels within TME from scRNAseq data.
