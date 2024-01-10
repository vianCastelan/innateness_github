# Mouse Innateness defined by transcriptomical gradient

> Gabriel Ascui<sup>1,2,3 *</sup>, Viankail Cedillo-Castelan<sup>1 *</sup>, Alba Mendis<sup>1</sup>, Eleni Phung<sup>1</sup>, Hsin-Yu Liu<sup>1</sup>, Greet Verstichel<sup>1</sup>, Shilpi Chandra<sup>1</sup>, Mallory P. Murray<sup>1,3</sup>, Michael Croft<sup>1</sup>, Hilde Cheroutre<sup>1</sup>, Mitchell Kronenberg<sup>1</sup>.

<sup>1</sup> La Jolla Institute for Immunology, La Jolla, California, US; <sup>2</sup> University of California San Diego, La Jolla, California, US; <sup>3</sup> Immunological Genome Project Consortium.

<sup>*</sup>: Equal contribution.

Corresponding author: <mitch@lji.org> 

> questions/issues/comments about code : <vcastelan@lji.org> or <gascui@lji.org>

## Abstract

Innate T cells, such as NKT cells, MAIT cells, &gamma;&delta; T cells and some intraepithelial T cells, are populations with diverse developmental pathways, antigen specificities and functional capacities, but they all share the ability to respond rapidly in TCR-dependent and cytokine-dependent but TCR-independent activation. Recently, a transcriptional program that explains a gradient of innateness has been described in human blood lymphoid populations. Here, using the Immunological Genome Project Consortium publicly available bulk RNA-seq and ATAC-seq datasets of several mouse lymphocyte populations, we constructed linear-mixed models of innateness for mouse lymphoid populations.  Natural Killer (NK) cells mark the highest end of the scale, as germline-encoded fully differentiated innate lymphocytes, whereas the other end is marked by naive CD4 and CD8 T cells, as the most adaptive populations. Pathway analysis shows the resulting innateness gradient to contain transcriptional programs related to NK cell functionality, chemotaxis and motility, all traits of innate T cells.  Applying our models to conventional CD4 or CD8 T cell transcriptional data assigned higher innateness scores to effector and effector memory populations over central memory T cells. A picture emerges, which indicates that for T cells innateness is acquired with some types of antigen-experience and parallels with a loss in expansion capacity and a gain in functional maturation ultimately leading to terminal differentiation. Our results also correlate higher innateness scores with lower levels of calcium-dependent T-cell activation, which we confirmed experimentally, and a higher dependence on protein kinase C phosphorylation pathways.  Therefore, these cells have a higher threshold or different requirements for antigen receptor-dependent activation.

## GitHub Repository

This [repository](github.com/viancastelan/innateness) has all you will need to reproduce the results of Ascui & Cedillo-Castelan et al. 2023. 

## ImmGEN datasets

[ImmGEN](immgen.org) bulk RNA-seq and ATAC-seq datasets are available online here: <https://www.immgen.org/Databrowser19/DatabrowserPage.html>

MAIT cell datasets will be uploaded shortly. ATAC-seq count matrix can be shared upon request. 

# Code 
-----
## Requirements

You will require `R` version `3.5.6` or higher and the following installed  packages:

| Package | Version |
| ------- | -------- |
| `tidyverse`| tidyverse_1.3.2  |
| `lme4`  | lme4_1.1-30 |
| `lmerTest` | lmerTest_3.1-3 |
| `DESeq2` | DESeq2_1.28.1 |
| `yaml` | yaml_2.3.5 |
| `scales` | scales_1.2.1 |

You will require `python3` and the the following installed libraries:

| Package | Version |
| ------- | -------- |
| yaml | version|
| pandas | version|
| sys | version |
| os | version |

### YAML config

The `config.yml` file contains the read and write directories. Make sure you are using has the correct ones assigned here.

### Running Linear-mixed models

Use the `lmm_job.sh` to run linear-mixed models. This script will call the `rna_m_lmm.R` script to generate linear mixed models based on RNA-seq transcript expression and PC1. Here you can also filter for how lower or higher levels of `beta` and `p-value`. This script is meant for a Torque-based job submission in a high-computing cluster, but can be adapted if necessary.

```
> bash lmm_job.sh -s mouse -p 0.05 -b 10
```

This script will require a count matrix with the ImmGEN dataset to generate linear models. The `R` script will call `DESeq2` to generate a PCA analysis, from which PC1 values will be evaluated against gene expression for each gene. Next, using the `lme4` package, linear models will be evaluated. 

Linear models will follow the following formula:

$$X_i \sim PC1_j + (1|celltype_j)$$
    
>Where *X* is the expression of a particular gene *i* in `celltype:` *j*. In this formula, `celltype` is considered as a fixed effect of the model.

This will generate a `results_mouse_beta_table.csv` file containing the following values: 

| gene | beta | pval | var | 
| -- | -- | -- | -- |
| *Ccl5* | 456.34 | 3.45e-40 | 45 |
| ... | ... | ... | ... |
| ... | ... | ... | ... |
| ... | ... | ... | ... |


## Calculating Innateness scores

Use the `innate_score.py` script to generate an innateness score per each celltype of a transcriptome count table. This could either be a microarray dataset or bulk/single cell RNA-sequencing datasets. 

This script will require the following arguments: 

`-c` or `--countmatrix` : Transcriptome countmatrix to calculate innateness score for each sample.  
`-b` or `--betatable` : Beta levels table file  
`-o` or `--output` : Output directory  

To run this script over the same ImmGEN bulk RNA-seq data, run the following code:

```
> python3 python/innate_score.py -c data/countmatrix/immgen_ULI_RNAseq.csv \
       -b output/b_levels/results_mouse_beta_table.tsv \
       -o output/b_scores/
```


Innateness score per cell type are calculated as for the following equation:

$$ \sum_{k=1}^n \beta_i\times X_i = innateness\: score $$ 

----
## Visualization

To visualize the results, you can use many functions written in `R` in this repository. 

### Principal Component gradient of innateness

Run the following script to generate Principal Component coordinates and plot where each cell type samples is on PC1 for the transcriptomical datasets. 

```
> Rscript R/PCA_plots.R
```

![PCA plot](output/figures/pca_plot.png?raw=true "PCA plot for RNA-seq data at 1000 MVG")


[//]: <> (This is a comment)


### Volcano of $\beta$-levels

run the following script: 

```
>Rscript R/volcano_beta.R
```


&beta;-levels are calculated according to each mixed linear model (first formula).

### barplot of $\beta$-levels

run the following script: 

```
> Rscript R/barplot.R
```

### Individual gene expression per tissue

Ues the following script to look at the expression of individual genes in the ImmGen dataset: `R/plot_gene.R`. 

This script will ask for the following variables: 

- gene 
- tissue


```
> Rscript R/plot_gene.R --help
> Rscript R/plot_gene.R --gene=Cxcr6 --tissue=spleen
```


### ShinyApp

We have deployed a ShinyApp here: <https://shinnyapp.io>.

### ATAC-seq motif enrichement: ChromVAR 

```
> chromVAR.Rmd
```

### T helper innateness separation and hierarchical clustering

First, run the models with the specific script to generate the models using T helper ranks instead of PC1 coordinates for all relevant cell types. 

```
> Rscript R/rna_m_thelper_lmm.R
```

Afterwards, run the `thelper.R` script to generate homologous genes based on this [report](https://doi.org/10.1261%2Frna.075929.120). 

```
> Rscript R/thelper.R
```

The following script will generate Venn Diagram comparing resulting Th models between each other and the original innateness models generated above. 

```
> Rscript R/thelper_comparison.R
```
>Modify the filtering variables in this script to be more or less stringent on significance:

```{r}
## filter
beta_th1 <- na.omit(beta_th1[beta_th1$pval < 0.01 & abs(beta_th1$beta) > 50,])
beta_th2 <- na.omit(beta_th2[beta_th2$pval < 0.01 & abs(beta_th2$beta) > 50,])
beta_th17 <- na.omit(beta_th17[beta_th17$pval < 0.01 & abs(beta_th17$beta) > 10,])
```

Finally, to generate heatmaps with heirarchical clustering run this final script:

```
> Rscript R/heatmaps.R
```

### ATAC-seq tracks

Tracks were generated with `Gviz` R package from Bioconductor. 

Running these will require some additional packages:



```
> Rscript R/gviz_tracks.R
```

## Issues

Please report any issues to gascui@lji.org or vcastelan@lji.org, or preferably using the Github issues tab here: <https://github.com/vcastelan/innateness>

## Citation

## Data

<https://geo.ncbi.org> 

<https://zenodo.org>