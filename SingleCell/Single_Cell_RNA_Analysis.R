#1. Setting up your workspace
#Installing the packages necessary for the Single Cell Analysis

install.packages("dplyr")
install.packages("SeuratObject")
install.packages("Seurat")
install.packages("patchwork")

#Loading the package libraries

library(dplyr)
library(SeuratObject)
library(Seurat)
library(patchwork)

# Load the PBMC dataset from the link "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz" 
#Make sure you have the file unzipped and you have the correct data directory to read the file
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#2. Pre-processing workflow

#Quality Control and Filtering Cells
#Mitrochondrial RNAs and low counts indicate problems with sample processing and cell health


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#After looking at the data distribution, we filter out cells expressing less than 200 RNA types
#or more than 2500 as well as cells having more than 5% mitochondrial contributions. 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalizing the Data-After Quality Control data is normalized to log Scale.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify the highly variable Features(i.e. top 2000 variable genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data-Before dimensional reduction we scale the data so that highly expressed genes do not dominate

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#3. Dimensional Reduction to simplify the dataset while retaining its essential components

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#Determine the dimensionality of the dataset
#What dimensions would allow the dataset to retain its complexity while still simplifying it
ElbowPlot(pbmc)
#There is little to no variation around the dimensions 9-10 for Principle Component Analysis

#4. Clustering the Cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Nonlinear Dimensional Reduction with UMAP and t-SNE
pbmc <- RunUMAP(pbmc, dims = 1:10)
# individual clusters
DimPlot(pbmc, reduction = "umap")

#5. Finding differentially expressed features as cluster biomarkers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#Measuring the specificity of the markers
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Visualizing differential gene expression between clusters for genes of interest
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

#Picking specific features for visualization across the clusters
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

#6. Assigning cell type identity to clusters

#Generating a heatmap for top10 markers
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#Check which markers correspond to which cell type in the literature. 
#There are additional programs that also aid in this.
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)