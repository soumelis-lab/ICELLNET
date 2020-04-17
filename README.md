# ICELLNET: A transcriptome-based framework to dissect intercellular communication
---
Floriane Noël and Lucile Massenet-Regad 
---

This repository hosts the source code corresponding to the method described in [Noël et al, 2020](https://www.biorxiv.org/content/10.1101/2020.03.05.976878v1) to infer intercellular communication networks and dissect intercellular communication between multiples cell types based on their transcriptomic profiles.

To install it, the easiest way is to use the `R` package `devtools` and its function `install_github`. If you don't have all the dependancies needed to use ICELLNET package, run the commands below:  

    install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist")) ##Installs devtools and the icellnet CRAN dependancies

    if (!requireNamespace("BiocManager", quietly = TRUE)) # Installs icellnet Bioconductor dependancies 
        install.packages("BiocManager")
    BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))
    
Then you just have to load `devtools` package and run the command below:

    library(devtools)
    install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")

Once all the dependencies are downloaded and loaded, you can load the ‘icellnet’ package.    
    
    library(BiocGenerics)
    library("org.Hs.eg.db")
    library("hgu133plus2.db")
    library("annotate")
    library(jetset)
    library(readxl)
    library(psych)
    library(GGally)
    library(gplots)
    library(ggplot2)
    library(RColorBrewer)
    library(data.table)
    library(grid)
    library(gridExtra)
    library(ggthemes)
    library(scales)
    library(rlist)
    library(icellnet)
    
Examples on how to use `icellnet` package functions can be found in the Vignette.
