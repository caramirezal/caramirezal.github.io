
## Dependencies
library(rmarkdown)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)

set.seed(333)


## Loading example data
pbmc_url  <- 'https://github.com/caramirezal/caramirezal.github.io/blob/master/bookdown-minimal/data/pbmc_10X_250_cells.seu.rds?raw=true'
seurat <- readRDS(url(pbmc_url))
seurat

