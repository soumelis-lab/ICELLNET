# Dissection of intercellular communication using the transcriptome-based framework ICELLNET
---
Floriane Noël and Lucile Massenet-Regad 
---

This repository hosts the source code corresponding to the method described in [Noël, F., Massenet-Regad, L., Carmi-Levy, I. et al. ](https://www.nature.com/articles/s41467-021-21244-x) to infer intercellular communication networks and dissect intercellular communication between multiples cell types based on their transcriptomic profiles.

## Package installation

To install it, the easiest way is to use the `R` package `devtools` and its function `install_github`. If you don't have all the dependancies needed to use ICELLNET package, run the commands below:  

    install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist")) ##Installs devtools and the icellnet CRAN dependancies

    if (!requireNamespace("BiocManager", quietly = TRUE)) # Installs icellnet Bioconductor dependancies 
        install.packages("BiocManager")
    BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))
    
Then you just have to load `devtools` package and run the command below:

    library(devtools)
    install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")

Once all the dependencies are downloaded and loaded, you can load the ‘icellnet’ package.    
Examples on how to use `icellnet` package functions can be found in the Vignette.

## News - ICELLNET
Last package update: 2021-03-31

**New on ICELLNET:**
- ICELLNET ligand /receptor database include now **543 human ligand/receptor interactions** manually curated, with a particular interest on immune checkpoints, cytokines, and chemokines.
- ICELLNET includes now new features to handle cell-cell communication studies from single-cell RNAseq datasets. See [Vignette](https://github.com/soumelis-lab/ICELLNET/blob/master/Vignette.md) and use case [Exemple2_scRNAseq.md](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple2_scRNAseq.md) for details.



