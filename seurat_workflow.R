### Installing Packages----
install.packages("Seurat")
install.packages("tidyverse")

### Loading Libraries----
library(Seurat)
library(SeuratObject)
library(tidyverse)

### Reading Data----
NML1 <- Read10X(data.dir = "GSM3891621/NML1")
NML2 <- Read10X(data.dir = "GSM3891621/NML2")
NML3 <- Read10X(data.dir = "GSM3891621/NML3")

### Creating Seurat Object----
NML1 <- CreateSeuratObject(counts = NML1, project = "NML1", min.cells = 3, min.features = 200)
NML2 <- CreateSeuratObject(counts = NML2, project = "NML2", min.cells = 3, min.features = 200)
NML3 <- CreateSeuratObject(counts = NML3, project = "NML3", min.cells = 3, min.features = 200)

# View Metadata
view(NML1@meta.data)
view(NML2@meta.data)
view(NML3@meta.data)

# Saving .RDS
saveRDS(NML1, file = "NML1_seurat.rds")
saveRDS(NML2, file = "NML2_seurat.rds")
saveRDS(NML3, file = "NML3_seurat.rds")

# Merge Seurat Objects 
merge_seurat <- merge(
  x = NML1,
  y = list(NML2, NML3),
  add.cell.ids = c("NML1","NML2", "NML3"),
  project = "NML_Merge"
)

# Save merged filed as .RDS

saveRDS(merge_seurat, file = "merge_seurat")
saveRDS(merge_seurat, file = "merge_NML")

### Standard Pre-processing Workflow----

# Adding MT% in MetaData
merge_seurat <- PercentageFeatureSet(merge_seurat, pattern = "^MT-", col.name = "percent.mt")
merge_seurat <- PercentageFeatureSet(merge_seurat, pattern = "^MT-", col.name = "percent.mt")
view(merge_seurat@meta.data)
range(merge_seurat$percent.mt) #checking the range of MT Percentage

### Selecting Cell for Further Analysis----

# Visualizing QC Metric as Vlnplot
VlnPlot(merge_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
merge_seurat <- subset(merge_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)


VlnPlot(merge_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


### Normalization----

merge_seurat <- NormalizeData(merge_seurat,
              normalization.method = "LogNormalize", scale.factor = 10000)


# Finding Highly Variable Features
merge_seurat <- FindVariableFeatures(merge_seurat, selection.method = "vst", nfeatures = 2000)

# Top Feature and Feature Plot
VariableFeatures(merge_seurat)
VariableFeaturePlot(merge_seurat)

# Data Scaling

merge_seurat <- ScaleData(merge_seurat)
all.gene <- row.names(merge_seurat)
merge_seurat <- ScaleData(merge_seurat, features = all.gene)
merge_seurat[["RNA"]]$scale.data

names(merge_seurat@commands) # View all the commands performed on data

# Save RDS file
saveRDS(merge_seurat, file = "merge_seurat.rds")
view(merge_seurat@meta.data)

### Perform PCA on Scaled Data (Linear Dimensional Reduction)----

merge_seurat <- RunPCA(merge_seurat)
DimPlot(merge_seurat, reduction = "pca", dims = c(1,2))

# Dimensions of DataSet
ElbowPlot(merge_seurat)
ElbowPlot(merge_seurat, ndims = 50)

### Cluster the Cells----

merge_seurat <- FindNeighbors(merge_seurat, dims = 1:20)
merge_seurat <- FindClusters(merge_seurat, resolution = 0.1)
merge_seurat <- FindClusters(merge_seurat, resolution = 0.3)

### Run UMAP----
merge_seurat <- RunUMAP(merge_seurat, dims = 1:20)
Plot1 <- DimPlot(merge_seurat, reduction = "umap", label = T, repel = T, group.by = "orig.ident")

view(merge_seurat@meta.data)


### Splitting Object into Lists----

merge_NML <- readRDS("merge_NML")
merge_seurat_list <- SplitObject(merge_NML, split.by = "orig.ident")

# Normalize Each Seurat Object
merge_seurat_list <- lapply(merge_seurat_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nFeatures = 2000) } )
# Selecting Features
features <- SelectIntegrationFeatures(object.list = merge_seurat_list)

# Finding Integration Anchor
merge_seurat_anchors <- FindIntegrationAnchors(object.list = merge_seurat_list, anchor.features = features)

# Integrating Data
merge_seurat_integ <- IntegrateData(anchorset = merge_seurat_anchors)

# Setting Default Assay for Downstream Analysis
DefaultAssay(merge_seurat_integ) <- "integrated"

# Running the Standard Workflow for Integrated Data
merge_seurat_integ <- ScaleData(merge_seurat_integ, verbose = FALSE)
merge_seurat_integ <- RunPCA(merge_seurat_integ, npcs = 50, verbose = FALSE)
merge_seurat_integ <- FindNeighbors(merge_seurat_integ, reduction = "pca", dims = 1:30)
merge_seurat_integ <- FindClusters(merge_seurat_integ, resolution = 0.3)
merge_seurat_integ <- RunUMAP(merge_seurat_integ, reduction = "pca", dims = 1:30)

# Visualization 
Plot2 <- DimPlot(merge_seurat_integ, reduction = "umap", label = T, repel = T, group.by = "orig.ident")
Plot1 + Plot2

# Save RDS file
saveRDS(merge_seurat_integ, file = "merge_seurat_integ.rds")
