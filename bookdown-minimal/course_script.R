
## Dependencies
library(rmarkdown)
library(dplyr)

library(Seurat)

set.seed(333)

## Loading example data
pbmc_url  <- 'https://github.com/caramirezal/caramirezal.github.io/raw/master/bookdown-minimal/data/pbmc_10X_500_cells.mtx.rds'
pbmc.mtx <- readRDS(url(pbmc_url))

dim(pbmc.mtx)


pbmc.seurat <- CreateSeuratObject(
  counts = pbmc.mtx,
  project = 'PBMC',
  assay = 'RNA',
  min.cells = 1,
  min.features = 1
)

pbmc.seurat[["percent.mt"]] <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")

pbmc.filtered <- subset(pbmc.seurat,
                        nFeature_RNA < 1250 &
                          nCount_RNA < 4000 &
                          percent.mt < 5)

pbmc.filtered <- FindVariableFeatures(pbmc.filtered, nfeatures = 1000)

pbmc.filtered <- NormalizeData(pbmc.filtered)

pbmc.filtered <- ScaleData(pbmc.filtered)

pbmc.filtered <- RunPCA(pbmc.filtered,
                        features = VariableFeatures(pbmc.filtered))

pbmc.filtered <- FindNeighbors(pbmc.filtered)
pbmc.filtered <- FindClusters(pbmc.filtered,
                              resolution = 0.1,
                              verbose = FALSE)

pbmc.filtered <- RunUMAP(pbmc.filtered,
                         dims = 1:7,
                         verbose = FALSE)

Idents(pbmc.filtered) <- pbmc.filtered$seurat_clusters

pbmc.degs <- FindAllMarkers(pbmc.filtered,
                            logfc.threshold = 1,
                            min.pct = 0.05,
                            min.cells.feature = 10,
                            verbose = FALSE)
