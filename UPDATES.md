#ICELLNET news:
Last package update: 2024-03-21
Last ligand-receptor interactions database update: 2023-11-08


## ICELLNET v2.1.0 (2023-11-08):
- Compatibility of ICELLNET functions with Seurat v5 and latest R versions (R> 4.1.x).
- Update of `LR.family.score()` to reintroduce the "barplot" vizualisation option.
- New use case study, showing how to identify cell type-specific interactions based on scRNAseq data.

- Specific functions modifications:
  - `sc.data.cleaning()` and  `Perc_exp_infos()` functions were adapted to work with Seurat v5 and latest R versions. 
  - `LR.family.score()`: plot argument (previously TRUE or FALSE) can take NULL (no plot, returns a table), "heatmap" or "barplot" values, in order to provide corresponding vizualisations.
  - n=1 argument was set as default in `gene.scaling()` function. 
  - CC.type and PC.type were set to "RNAseq"" as default in `icellnet.score()` function


## ICELLNET v2.0.0 (2023-11-08): 
- New structure and extension of ICELLNET database to consider more subunits of ligand and receptors (enable 4 ligand subunits and 5 receptor subunits)  - now contains 1669 interactions;
- New classification of ICELLNET database: "Antigen binding" family renamed "HLA recognition"; "Adhesion molecules" renamed "Cell adhesion". New families includes "Wnt pathway" and "Innate immune".
- Formatting functions provided to use CellphoneDB, matching the structure of ICELLNET databases: `CellPhoneDB_convert()` 
- Update of several functions to match the new structure of ICELLNET database: `name.lr.couple()`,`ligand.average()` (replacing previous ligand.average.RNAseq() and ligand.average.MA() functions), `receptor.average()` (replacing previous receptor.average.RNAseq() and receptor.average.MA() functions ), `LR.viz()`, `icellnet.score()`, `icellnet.ind.score()`
- Update of several visualisation functions: `LR.viz()`, `LR.family.score()`, `LR.family.score()`, `LR.heatmap()`, `LR.selection()`

- Specific functions modifications: 
    - `LR.viz()`: can now display the pairwise interaction score of several interactions at a time, provided in "int" arguments (replace "couple" argument). Also possible to visualise directly a family of interactions when using "family" input argument. This function requires gene symbols as inputs and cannot be used with microarray probe names.
    - `LR.family.score()`: family (replacing my.family argument) and db.couple (replacing db.name.couple) arguments are not required (defaults values provided). Arguments "cluster_rows" and "cluster_columns" are added (used by ComplexHeatmap for heatmap visualisation) 
    - `name.lr.couple()`: LR database becomes the only required input argument.
    - `network.create()`: not necessary anymore to rescale communication score from 1 to 10. 


## ICELLNET v1.3.0 (2023-06-29):
Updated functions with minor changes:  Correction in `gene.scaling()` function. note() removed in sc_data_cleaning function. Correction of LR.family.score function when considering only one cell type. Vignette update.
New package release. 

## ICELLNET v1.2.0 (2023-02-05):
ICELLNET can now handle spatial transcriptomic data!
Allow to use other assay than "RNA" assay in Seurat objects, such as "SCT" or "Spatial" assays.

## Update of ICELLNET database (2022-06-15):
ICELLNET database has been expanded to **1034 interactions**, focusing on adding **tumor microenvironment** related interactions.

## ICELLNET v1.0.0 - 2021-07-21:
- ICELLNET ligand/receptor database include now **752 human ligand/receptor interactions**, manually curated. In addition to immune checkpoints, cytokines, and chemokines, it now includes additional interactions from classified in an additional family "extracellular matrix", among others added interactions.
- New release (version 1.0.0): New graphical representations to handle cell-cell communication studies from single-cell RNAseq datasets. See [Exemple2_scRNAseq.md](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple2_scRNAseq.md) for details.

## First ICELLNET version - 2021-03-30:
- ICELLNET ligand/receptor database include now **543 human ligand/receptor interactions** manually curated, with a particular interest on immune checkpoints, cytokines, and chemokines.
- ICELLNET includes now new features to handle cell-cell communication studies from single-cell RNAseq datasets. See [Vignette](https://github.com/soumelis-lab/ICELLNET/blob/master/Vignette.md) and use case [Exemple2_scRNAseq.md](https://github.com/soumelis-lab/ICELLNET/blob/master/Exemple2_scRNAseq.md) for details.



