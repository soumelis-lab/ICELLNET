# ICELLNET: A transcriptome-based framework to dissect intercellular communication
---
Lucile Massenet-Regad and Floriane Noël
---

This repository hosts the source code corresponding to the method described in [Noël et al, 2020](https://www.biorxiv.org/content/10.1101/2020.03.05.976878v1) to infer intercellular communication networks and dissect intercellular communication between multiples cell types based on their transcriptomic profiles.

To install it, the easiest way is to use the `R` package `devtools` and its function `install_github`:

    install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist") ##Installs devtools and the icellnet dependancies

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate")
    
    library(devtools)
    install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
    
Examples on how to use `icellnet` package functions can be found in the Vignette.
