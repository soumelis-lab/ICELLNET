 
This vignette explains the use of the ICELLNET package and demonstrates typical workflows to dissect intercellular communication between multiple cell types, based on transcriptomic profiles.

---
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
---


# Introduction to ICELLNET R package  <a name="Introduction-to-ICELLNET-R-package"></a>

Cell-to-cell communication is at the basis of the higher-order organization observed in tissues and organisms, at steady state and in response to stress. The availability of large-scale transcriptomics datasets from several cell types has opened the possibility of **reconstructing cell-cell interactions based on co-expression of ligand-receptor pairs**.

We developed **ICELLNET**, a transcriptomic-based framework to **dissect cell communication in a global manner**. It integrates an original expert-curated **database of ligand-receptor interactions** taking into account multiple subunits expression. Based on transcriptomic profiles, ICELLNET package allows to compute **communication scores** between cells and provides **several visualization modes** that are helpful to dig into cell-cell interaction mechanism and extend biological knowledge. 

# ICELLNET ligand/receptor interaction database <a name="ICELLNET-ligand/receptor-interaction-database"></a>

We curated a comprehensive database of ligand-receptor interactions from the literature and public databases. This database takes into account the **multiple subunits** of the ligands and the receptors. Interactions have been classified into 6 families of communication molecules, with strong implication in inflammatory and immune processes: **Growth factors, Cytokines, Chemokines, Checkpoints, Notch family, Antigen binding**. Cytokines have been further classified into 7 subfamilies according to reference classifications essentially based on structural protein motifs: **type 1 cytokines, type 2 cytokines, IL-1 family, IL-17 family, TNF family, TGFb family and RTK cytokines**. 

Other interactions and classifications of molecules will be implemented.
The most recent version of ligand-receptor interaction database can always be downloaded [here](https://github.com/soumelis-lab/ICELLNET/blob/master/data/ICELLNETdb.tsv).
In R, you can visualize ICELLNET database and its structure: 

```{r db, echo=T}
db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
head(db)
```

You can use either all the database or restrict it by selecting some specific class of molecules (Cytokines, Growth factor etc..). Below, we show you how to restrict the study to cytokines, chemokines, and checkpoints, and how you can take into consideration subfamily of molecules.

```{r, echo=T}
summary(as.factor(db$Family)) # list of the different family of molecules considered in the database
db$Subfamily=db$Cytokine
summary(as.factor(db$Subfamily[which(db$Family=="Cytokine")])) # list of the different subfamily of cytokines considered in the database

#Restrict the database to some family of molecules 
my.selection.LR=c("Cytokine", "Chemokine", "Checkpoint")
db2 <- db[grepl(paste(my.selection.LR, collapse="|"),db$Classifications),] 
db.name.couple=name.lr.couple(db2, type="Family")
head(db.name.couple)

#Restrict the database to cytokines and consider the subfamilies of cytokines
my.selection.LR=c("Cytokine")
db3 <- db[grepl(paste(my.selection.LR, collapse="|"),db$Classifications),] #if you want to use all the database, do instead : db2=db
db.name.couple=name.lr.couple(db3, type="Subfamily")
head(db.name.couple)
```
Instead of using the ICELLNET database, it is also possible to use its own database as long as it is correctly formatted with specific columns as below. The Family and Subfamily colums correspond to two independant classifications (per family of molecules, or other) of your choice, but each interaction should fit only in one category of the classification (for example, an interaction cannot be classified in "type 1" and also "type 2" cytokines in the ICELLNET database). In the Classifications category, you should add all the terms used to classify the interaction : the one of Family, Subfamily, but also other words that can be used to select some specific interactions (for example "interleukin" in ICELLNET database).

|  Ligand 1 | Ligand 2  | Receptor 1  | Receptor 2  | Receptor 3  | Family | Subfamily | Classifications |
|---|---|---|---|---|---|---|---|
|   |   |   |   |   |   |   |   |
|   |   |   |   |   |   |   |   |


# Input data <a name="Input-data"></a>

ICELLNET pipeline considers transcriptomic profiles of a defined "central cell", that can correspond for exemple to a cell type in several biological conditions. ICELLNET will then allow to compare the communication channels used by the central cells in these different conditions with partner cells.
As partner cells, we can use Human Primary Cell Atlas, a public datasets of 745 transcriptomic profiles among 31 cell types generated with the same technology (Affymetrix microarray, hgu133plus2 platform), already processed. 

As partner cells, the user can also use other transcriptomic profiles instead of Human Primary Cell Atlas (see section  [How to format your own data to use ICELLNET package?](#How-to-format-your-own-data-to-use-ICELLNET-package?)for format details).

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
 
Here we describe the different stages of the ICELLNET package to compute intercellular communication scores: 

1. Selection of the genes coding for the ligands and the receptors in our database from the transcriptomic profiles of the central cell and the partner cells. 

2. Rescale gene expression to avoid communication score to be driven only by highly expressed genes.

3. Compute ICELLNET communication scores (`direction="out"` for outward communication, `direction="out"` for inward communication, see [How is the intercellular communication score computed?](##How-is-the-intercellular-communication-score-computed?) for more details)

4. Display different visualization modes to dissect intercellular communication scores 

![](pictures/ICELLNET_Figure2_V10.png)


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

### Intercellular communication network representation

This allows to visualize intercellular communication networks in a global manner through the function `network.create()`. In these directed graphs, nodes represent cell types, the width of the edges connecting two cell types is proportional to a global measure of the intensity of the communication between them and the arrows indicate the direction of communication. 


### Communication molecules distribution

The barplot representation (`LR.family.score()` function (with plot=T) allows to dissect the global scores at a level of class of molecules, and allows to identify patterns of co-expressed molecules from the same family. This layer of analysis helps the interpretation on a qualitative level.

Color legends for these functions can be found below (easy to adapt to study other family of molecules):

```{r, echo=T, warning=F,  fig.height=5, fig.width=12 }
    ## label and color label if you are working families of molecules already present in the database
# my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding") 
# family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  ,
#             "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "#12A039",  "other" = "#908F90",  "NA"="#908F90")
     
    ## label and color label if you are working with subfamilies of cytokines
my.family=c("type 1", "type 2", "IL1.", "IL17", "TNF","TGF","RTK")
family.col = c( "type 1"=  "#A8CF90", "type 2"= "#676766", "IL1."= "#1D1D18" ,
            "IL17" ="#66ABDF", "TNF" ="#12A039", "TGF" = "#156399", "RTK"="#AECBE3", "other" = "#908F90","NA"="#908F90")

```

### Individual communication scores distribution

The balloon plot (`LR.balloon.plot()` function) and heatmap (`LR.heatmap()` function)  are the deepest level of representation of the communication, displaying the most contributing ligand/receptor pairs to the communication score. This allows to identify specific individual interactions that can drive the intercellular communication and should be confirmed experimentally. 
*topn* parameter set the top n highest interactions to display (ex: `topn=10` means display the 10 highest contribution of specific interactions to the communication scores)
*thresh* parameter allows to set a threshold, in order to display interaction contributing to communication scores only above this value. (ex: thresh=20, only interactions with an individual score above 20 will be visualized). Looking at individual scores distribution is helpful to identify relevant threshold values for you case.
*sort.by* parameter ("sum" as default, can be set to "var") to ran the interactions either as the most contributing interactions in all conditions ("sum"), or the most different interactions among conditions ("var", ranking by computing the variance among scores for each interaction)

Same color legend used as `LR.family.score()` function above.



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


To load the most recent version of the ligand/receptor database, you should run the command below (last updated:  28/01/2021): 

    db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))


# How to format your own data to use ICELLNET package? <a name="How-to-format-your-own-data-to-use-ICELLNET-package?"></a>

## Data files format
**For the central cell:** It can be any transcriptomic profile data of one cell type. For **RNA-seq data**, the dataset should be a datframe annotated with gene symbol as rownames. For **microarray data**, the ICELLNET functions are adapted to handle Affymetrix Human Genome U133 Plus 2.0 Array annotation. Nevertherless, if the dataset have been generated with an other Affymetrix technology, you have 2 possibilities to adapt the tool : a) Annotate your data with gene symbol before using ICELLNET and then consider your data as "RNA-Seq" for CC.type argument. b) adapt the R code of the `db.hgu133plus2()` function to have the adapted annotation conversion when using ICELLNET. Gene annotations should be set as rownames. 


**For the partner cell (if you don't want to use Human Primary Cell Atlas as partner cells):** This can be interesting for example if you possess transcriptomic data of several cell types of the same sample, to see how they interact together. As for the central cell, the transcriptomic profiles should be correctly formated (see previous paragraph above for more information). If your transcriptomic profiles are annotated with gene symbol, PC.type should be set to "RNA-seq" (even if your data come from microarray technology). 

## Target files format
You should define two dataframes as target files, one corresponding to the central cell, and the other one corresponding to the partner cells. These dataframe usually describes the different samples. 

**PC.target** should contains at least an 'ID' column including the name of the samples (usually rownames(PC.target) or colnames(PC.data) ), and a 'Class' column corresponding to a classification of your different samples included in PC.target, such as a cell type classification. The different categories included in the 'Class' column will define the different partner cells in the graphs.


# Use cases exemples <a name="Use-cases-exemples"></a>


- [Case study 1](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple1_CAF.md): dissect intercellular commmunication of Cancer Associated Fibroblasts subsets. Show how to apply ICELLNET pipeline on transcriptomic profiles from 2 CAF-subsets, and how to restrict use of icellnet database to cytokines only (or other family of molecules).

- [Case study 2](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple2_scRNAseq.md): Application of ICELLNET pipeline to scRNAseq from a Seurat object to infer intercellular communication between clusters.



# Software information <a name="Software-information"></a>

```
session.Info()
```
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.4
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] circlize_0.4.12       ComplexHeatmap_2.8.0  SeuratObject_4.0.2    Seurat_4.0.2         
 [5] gridExtra_2.3         icellnet_1.00         dplyr_1.0.7           ggplot2_3.3.5        
 [9] jetset_3.4.0          hgu133plus2.db_3.13.0 org.Hs.eg.db_3.13.0   AnnotationDbi_1.54.1 
[13] IRanges_2.26.0        S4Vectors_0.30.0      Biobase_2.52.0        BiocGenerics_0.38.0  

loaded via a namespace (and not attached):
  [1] readxl_1.3.1           plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2        
  [5] splines_4.1.0          listenv_0.8.0          scattermore_0.7        GenomeInfoDb_1.28.1   
  [9] digest_0.6.27          foreach_1.5.1          htmltools_0.5.1.1      fansi_0.5.0           
 [13] magrittr_2.0.1         memoise_2.0.0          doParallel_1.0.16      tensor_1.5            
 [17] cluster_2.1.2          ROCR_1.0-11            globals_0.14.0         Biostrings_2.60.1     
 [21] matrixStats_0.59.0     spatstat.sparse_2.0-0  colorspace_2.0-2       blob_1.2.1            
 [25] ggrepel_0.9.1          crayon_1.4.1           RCurl_1.98-1.3         jsonlite_1.7.2        
 [29] spatstat.data_2.1-0    iterators_1.0.13       survival_3.2-11        zoo_1.8-9             
 [33] glue_1.4.2             polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.38.0       
 [37] XVector_0.32.0         leiden_0.3.8           GetoptLong_1.0.5       shape_1.4.6           
 [41] future.apply_1.7.0     abind_1.4-5            scales_1.1.1           DBI_1.1.1             
 [45] miniUI_0.1.1.1         Rcpp_1.0.7             viridisLite_0.4.0      xtable_1.8-4          
 [49] clue_0.3-59            tmvnsim_1.0-2          reticulate_1.20        spatstat.core_2.1-2   
 [53] bit_4.0.4              htmlwidgets_1.5.3      httr_1.4.2             RColorBrewer_1.1-2    
 [57] ellipsis_0.3.2         ica_1.0-2              pkgconfig_2.0.3        farver_2.1.0          
 [61] uwot_0.1.10            deldir_0.2-10          utf8_1.2.1             tidyselect_1.1.1      
 [65] labeling_0.4.2         rlang_0.4.11           reshape2_1.4.4         later_1.2.0           
 [69] munsell_0.5.0          cellranger_1.1.0       tools_4.1.0            cachem_1.0.5          
 [73] cli_3.0.1              generics_0.1.0         RSQLite_2.2.7          ggridges_0.5.3        
 [77] stringr_1.4.0          fastmap_1.1.0          goftest_1.2-2          bit64_4.0.5           
 [81] fitdistrplus_1.1-5     purrr_0.3.4            RANN_2.6.1             KEGGREST_1.32.0       
 [85] pbapply_1.4-3          future_1.21.0          nlme_3.1-152           mime_0.11             
 [89] compiler_4.1.0         rstudioapi_0.13        plotly_4.9.4           curl_4.3.2            
 [93] png_0.1-7              spatstat.utils_2.1-0   tibble_3.1.2           stringi_1.7.3         
 [97] RSpectra_0.16-0        lattice_0.20-44        Matrix_1.3-4           psych_2.1.6           
[101] vctrs_0.3.8            pillar_1.6.1           lifecycle_1.0.0        GlobalOptions_0.1.2   
[105] spatstat.geom_2.1-0    lmtest_0.9-38          RcppAnnoy_0.0.18       data.table_1.14.0     
[109] cowplot_1.1.1          bitops_1.0-7           irlba_2.3.3            httpuv_1.6.1          
[113] patchwork_1.1.1        R6_2.5.0               promises_1.2.0.1       KernSmooth_2.23-20    
[117] parallelly_1.25.0      codetools_0.2-18       MASS_7.3-54            assertthat_0.2.1      
[121] rjson_0.2.20           withr_2.4.2            sctransform_0.3.2      mnormt_2.0.2          
[125] GenomeInfoDbData_1.2.6 mgcv_1.8-36            rpart_4.1-15           tidyr_1.1.3           
[129] Cairo_1.5-12.2         Rtsne_0.15             shiny_1.6.0 