This directory contains the part 1 of the analysis: using Seurat package in R to do single cell RNA- Seq analysis

Table of Content
Introduction
Preparation
Import Seurat package
Step 1. Create a Seurat object
Step 2. Quality control
Step 3. Normalization
Step 4. Feature selection for following heterogeneity analysis
Step 5. Data scaling
(Optional and advanced) Alternative step 3-5: to use SCTransform
Step 6. Linear dimension reduction using principal component analysis (PCA)
Step 7. Non-linear dimension reduction for visualization
Step 8. Cluster the cells
Step 9. Annotate cell clusters
Step 10. Pseudotemporal cell ordering
Step 11. Save the result


Introduction
After getting the scRNA-seq data of your samples, you will want to analyze it properly.

Multiple toolkits and analytic frameworks have been developed to facilitate scRNA-seq data analysis. These options include but are not limit to Seurat, developed by Rahul Satija's Lab in R, and scanpy, developed by Fabian Theis's Lab in Python. Both toolkits provide functions and rich parameter sets that serve most of the routine analysis that one usually does on scRNA-seq data. However, one should be aware that these analytic frameworks do not cover all the interesting analyses that one can do when analyzing data. It is also important to get to know other tools for scRNA-seq data analysis.

Since this is a tutorial for beginners, we will mostly introduce how to use Seurat to analyze your scRNA-seq data in R. At the end, we will also mention some other additional tools (e.g. presto, destiny, Harmony, simspec, etc.), which provide additional functionalities that you may miss if you only use Seurat. In the most recent update, we also provide the briefly example of some commonly used advanced analysis, such as RNA velocity.

Preparation
This tutorial assumes that the sequencing data preprocessing steps, including base calling, mapping and read counting, have been done. 10x Genomics has its own analysis pipeline Cell Ranger for data generated with the 10x Genomics Chromium Single Cell Gene Expression Solution. At the end of the Cell Ranger pipeline, a count matrix is generated. If your scRNA-seq data is generated using another technology (e.g. well-based experiments using Smart-Seq2 and others), the Cell Ranger pipeline is likely unapplicable, and you will have to find another solution to generate the count matrix.

As part of this tutorial, we are providing two data sets (DS1 and DS2), both generated using 10x Genomics and preprocessed using Cell Ranger. They are both public scRNA-seq data of human cerebral organoids and are part of the data presented in this paper. The first part of this tutorial, which includes most of the general analysis pipeline, is based on DS1, while the second part, which focuses on data integration and batch effect correction, is based on both data sets.

As a test for yourself, please try to apply what is learned in the first part to DS2 and only then continue with part 2 of the vignette. This will also give you an idea which types of cells are in DS2 and how comparable it is to DS1, before doing any data integration of both data sets.

Now let's start Part 1
Step 0. Import Seurat package
First of all, please make sure that Seurat is installed in your R.

library(Seurat)
This imports your installed Seurat package into your current R session. No error should be seen but some verbose information is likely. If it warns you that the package is unavailable, please install Seurat first

install.packages("Seurat")
library(Seurat)
