#!/bin/R
#R script to run sc-RNA data analysis using seurat package in R

#load seurat package:
library(Seurat)

#Step 1. Create a Seurat object
#Seurat implements a new data type which is named 'Seurat'. the first step is to read in the data and create a Seurat object. 
counts <- Read10X(data.dir = "data/DS1/")
seurat <- CreateSeuratObject(counts, project="DS1")

#The Read10X function does is to read in the matrix and rename its row names and col names by gene symbols and cell barcodes, respectively.

#manually loading the data and creating seurat object:
library(Matrix)
counts <- readMM("data/DS1/matrix.mtx.gz")
barcodes <- read.table("data/DS1/barcodes.tsv.gz", stringsAsFactors=F)[,1]
features <- read.csv("data/DS1/features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes
seurat <- CreateSeuratObject(counts, project="DS1")


#Step 2. Quality control: After creating the Seurat object, the next step is to do quality control on the data. The most common quality control is to filter out

#Need to filter out cells based on these conditions: 
#Cells with too few genes detected. They usually represent cells which are not sequenced deep enough for reliable characterization.
#Cells with too many genes detected. They may represent doublets or multiplets (i.e. two or more cells in the same droplet, therefore sharing the same cell barcode).
#Cells with high mitochondrial transcript percentage.

#calculate mitochondial transcript percentages:

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

# look at the distribution by creating a violin plot for each of the metrics.
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Or if you don't like the dots (individual cells)

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

#Due to the correlation of gene number and transcript number, we only need to set a cutoff to either one of these metrics, combined with an upper threshold of mitochondrial transcript percentage, for the QC. For instance, for this data set, a detected gene number between 500 and 5000, and a mitochondrial transcript percentage lower than 5% would be quite reasonable, but it is fine to use different thresholds.
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
#It is worth to mention that sometimes more QC may need to be applied.

#Normalization:
#A normalization step, aiming to make gene expression levels between different cells comparable, is therefore necessary. 

seurat <- NormalizeData(seurat)

#Feature selection for following heterogeneity analysis:
#In Seurat, or more general in scRNA-seq data analysis, this step usually refers to the identification of highly variable features/genes, which are genes with the most varied expression levels across cells.

seurat <- FindVariableFeatures(seurat, nfeatures = 3000)

#One can visualize the result in a variable feature plot, but this is optional.
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

# A scaling is applied to the data using the selected features, just like one usually does in any data science field.
seurat <- ScaleData(seurat)

#At this step, one can also remove unwanted sources of variation from the data set by setting the parameter var.to.regress. For instance,
seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

#For scRNA-seq data, the linear dimension reduction mostly refers to principal component analysis, short PCA.
seurat <- RunPCA(seurat, npcs = 50)

#cluster the cells:
#First of all, a k-nearest neighbor network of cells is generated. Every cells is firstly connected to cells with the shortest distances, based on their corresponding PC values. Only cell pairs which are neighbors of each other are considered as connected. Proportion of shared neighbors between every cell pairs is then calculated and used to describe the strength of the connection between two cells. Weak connections are trimmed. This gives the resulted Shared Nearest Neighbor (SNN) network. In practice, this is very simple in Seurat.

seurat <- FindNeighbors(seurat, dims = 1:20)
#With the network constructed, the louvain community identification algorithm is applied to the netowkr to look for communities in the network, i.e. cell groups that cells in the same group tend to connect with each other, while connections between cells in different groups are sparse.
seurat <- FindClusters(seurat, resolution = 1)

#visualize the clustering result using the tSNE and UMAP embeddings that are generated before.

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2

#annotate the cell clusters:
#The easiest to visualize expression of marker genes of interest across cell clusters is probably by a heatmap.

ct_markers <- c("MKI67","NES","DCX","FOXG1", # G2M, NPC, neuron, telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon related
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related
DoHeatmap(seurat, features = ct_markers) + NoLegend()

#we should firstly identify cluster markers for each of the cell cluster identified. In Seurat, this can be done using the FindAllMarkers function

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#First of all, cells of interest are extracted. Afterwards, we re-identify highly variable genes for the subset cells, as genes representing differences between dorsal telencephalic cells and other cells are no longer informative

seurat_dorsal <- subset(seurat, subset = RNA_snn_res.1 %in% c(0,2,5,6,10))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)

#We can try to reduce cell cycle effect by excluding cell cycle related genes from the identified highly variable gene list.

VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unlist(cc.genes))

#We can then check how the data look like, by creating a new UMAP embedding and do some feature plots

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)

#the ScaleData function has the option to include variables representing sources of unwanted variations. We can try to use that to further reduce the cell cycle influence; but before that, we need to generate cell-cycle-related scores for every cell to describe their cell cycle status.

seurat_dorsal <- CellCycleScoring(seurat_dorsal,
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score", "G2M.Score"))
#We can then check how the data look like, by creating a new UMAP embedding and do some feature plots

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)

#Now let's try to run diffusion map to get the cells ordered.

library(destiny)
dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:20])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(dpt$dpt)
FeaturePlot(seurat_dorsal, c("dpt","GLI3","EOMES","NEUROD6"), ncol=4)
