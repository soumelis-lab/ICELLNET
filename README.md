
# Dissection of intercellular communication using the transcriptome-based framework ICELLNET

This repository hosts the source code corresponding to the method described in [Noël, F., Massenet-Regad, L., Carmi-Levy, I. et al. ](https://www.nature.com/articles/s41467-021-21244-x) to infer intercellular communication networks and dissect intercellular communication between multiples cell types based on their transcriptomic profiles.

---
## ICELLNET key features: 
- Contains a **manually curated ligand-receptor interactions database**: curated exclusively from **human** studies, ICELLNET database takes into account the **multiple** subunits of ligand and receptor complexes (more than 1600 interactions)
- **Versatile** tool applicable on a wide range of transcriptomic technologies (microarray, bulk RNAseq, scRNAseq, spatial transcriptomics)
- Several visualization modes of cell-cell communication analysis results
---

## New on ICELLNET v2 - corresponding paper [link](https://pubmed.ncbi.nlm.nih.gov/38490248/)
- New structure and extension of ICELLNET database up to 1669 interactions.
- Possibility to ICELLNET analysis using other databases such as CellPhoneDB as input.
- Updates of several functions to match the new structure of ICELLNET database.

Additional details on updates and previous versions can be found [here](https://github.com/soumelis-lab/ICELLNET/blob/master/UPDATES.md).


*Remark: ICELLNET v2.0.0 requires R version < 4.1. For latest R versions (> 4.1.x), please use ICELLNET v2.2.x.*


## Information related to ICELLNET and tutorials

- [Detailed documentation](https://github.com/soumelis-lab/ICELLNET/blob/master/Vignette.md)
- ICELLNET analysis on scRNAseq data: [Tutorial](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Exemple2_scRNAseq.md) 
- ICELLNET analysis on sorted cell-population transcriptomic data: 
  - [Script](https://github.com/soumelis-lab/ICELLNET/issues/12): analyze cell-cell communication between different cell populations from the same sample
  - [Tutorial](https://github.com/soumelis-lab/ICELLNET/blob/master/exemples_tutorials/Exemple1_CAF.md):  infer cell-cell communication from the transcriptome of a single cell type of interest with the Primary Cell Atlas (Database of reference transcriptome of 33 cell types). 
  

## Installation

### Installation of dependencies:

    install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist")) ##Installs devtools and the icellnet CRAN dependancies

    if (!requireNamespace("BiocManager", quietly = TRUE)) # Installs icellnet Bioconductor dependancies 
        install.packages("BiocManager")
    BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))
    
### Installation of ICELLNET package:

    library(devtools)
    install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")

Once all the dependencies are downloaded and loaded, you can load the ‘icellnet’ package.    


## Citation
If you use our human curated ligand-receptor interaction database or the ICELLNET method, please cite our papers: 

Noël, F., Massenet-Regad, L., Carmi-Levy, I. et al. Dissection of intercellular communication using the transcriptome-based framework ICELLNET. Nat Commun 12, 1089 (2021). https://doi.org/10.1038/s41467-021-21244-x [link](https://www.nature.com/articles/s41467-021-21244-x)

Massenet-Regad L, Soumelis V. ICELLNET v2: a versatile method for cell-cell communication analysis from human transcriptomic data. Bioinformatics. 2024 Mar 4;40(3):btae089. doi: 10.1093/bioinformatics/btae089. PMID: 38490248. [link](https://pubmed.ncbi.nlm.nih.gov/38490248/)


