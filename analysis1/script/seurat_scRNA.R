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

#
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
