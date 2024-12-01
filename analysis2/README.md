**`Introduction: jointly analyze multiple scRNA-seq data sets`**

Although some experimental strategy, e.g. cell hashing, as well as computational demultiplexing methods such as demuxlet and scSplit to some extend allow pooling multiple samples together for the scRNA-seq library preparation and sequencing, it is unavoidable that certain steps, e.g. tissue dissociation, would have to be done separately for diffent samples. Therefore, just like when dealing with bulk RNA-seq data, batch effect is usually a critical confounder of the result that one has to resolve.

In this part of the tutorial, several scRNA-seq integration methods would be introduced. We will use DS1 which has been described in the first part of the tutorial, together with DS2 which you should have analyzed following this vignette. 
``The datasets are available in scRNA-Seq-Seurat-/datasets/``

***`Step 0. Load data`***
Let's start with importing Seurat and load the saved Seurat object.

```r
library(Seurat)
library(dplyr)
library(patchwork)
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
seurat_DS2 <- readRDS("DS2/seurat_obj_all.rds")
```

***`Step 1. Merge the two data sets`***
First of all, there is some chances that batch effect is small so that no integration is necessary. Therefore, we should firstly take a look at the two data sets by simply merging them together.
```r
seurat <- merge(seurat_DS1, seurat_DS2) %>%
FindVariableFeatures(nfeatures = 3000) %>%
ScaleData() %>%
RunPCA(npcs = 50) %>%
RunUMAP(dims = 1:20)
plot1 <- DimPlot(seurat, group.by="orig.ident")
plot2 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
plot1 + plot2 + plot_layout(widths = c(1.5, 2))
```
![umap_merged_datasets](https://github.com/user-attachments/assets/2cfb8dab-91dc-465b-8791-5f455e8faeea)



Obviously, the two data sets separate from each other on the embedding. However, the marker expression patterns suggest that the two data sets indeed share quite many cell types. Ideally, cells of the same cell type in the two data sets should be mixed with each other. However, because of the batch effect, this is not happening. So we need to do data integration. What we hope is that after the integration, cells of the same cell type in the two data sets intermix, while cells of different cell types/states still separate.

Here we will try different methods, including

Seurat

Harmony

LIGER

MNN

RSS to BrainSpan

CSS

***`Step 2-1. Data integration using Seurat`***

Seurat has its own data integration procedure implemented. In brief, it firstly applies canonical correlation analaysis (CCA) to the data sets that need to be integrated, rotating them separately so that the covariance of the two data sets is maximized. In other words, Seurat uses CCA to find the way maximizing the similarities between data sets. Next, Seurat introduces an anchoring mechanism, looking for cell anchors in the two data sets. Cell anchors are cell pairs with each cell in a different data set. The two cells are one of the nearest neighbors of each other in the CCA space, while the nearest neighbors of one cell in its own data set also tend to be neighbors of the nearest neighbors of the other cell of the cell pair. The two anchored cells are seen as corresponding cells from one data set to the other, and an integration procedure is then applied by subtracting expression of one data set by the transformation matrix calculated by comparing the anchoring cell pairs in the two data sets. People interested in its detailed methodology can read its paper.

To do integration using Seurat, one needs to firstly normalize and identify highly variable genes for each of data set to be integrated (which should have been done). If it hasn't been done, do it first:
```r
seurat_DS1 <- NormalizeData(seurat_DS1) %>% FindVariableFeatures(nfeatures = 3000)
seurat_DS2 <- NormalizeData(seurat_DS2) %>% FindVariableFeatures(nfeatures = 3000)
```


Next, we identify anchors of data sets. At this step, Seurat takes a list of Seurat objects as the input. Please note that Seurat allows integration of more than two samples. One just needs to put them into a list.
```r
seurat_objs <- list(DS1 = seurat_DS1, DS2 = seurat_DS2)
anchors <- FindIntegrationAnchors(object.list = seurat_objs, dims = 1:30)
```

P.S. The dims parameter determines the number of CC components to take into account, and one should try different values to fine-tune the results.

Next, the identified anchor set is passed to the the IntegrateData function to do the expression level correction.

```r
seurat <- IntegrateData(anchors, dims = 1:30)
```


Running the ```IntegrateData``` function creates a new ```Assay``` object (by default it is called integrated), where the batch-corrected expression matrix is stored. The uncorrected values are not lost, but store in the original Assay object (called RNA by default). The default assay of the resulted Seurat object is automatically set to integrated, but one can switch to the other one by using ```e.g. DefaultAssay(seurat) <- "RNA".```

Next, we just take the corrected Seurat object and re-run the procedure in Part 1, except for the first two steps (normalization and highly variable gene identification) which should be skipped here.
```r
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, dims = 1:20) %>% FindClusters(resolution = 0.6)
```


# You may also want to save the object
```r
saveRDS(seurat, file="integrated_seurat.rds")
```

Please be aware, that while the tSNE/UMAP embedding and clustering should be done with the integrated assay, the corrected values are no longer very reliable as the quantitative measure of gene expression. It is recommended that for the other analysis such as cluster marker identification and visualization, to use the uncorrected expression values instead, by setting the DefaultAssay back to RNA
```r
DefaultAssay(seurat) <- "RNA"
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```

![image](https://github.com/user-attachments/assets/af0e4f93-910f-40dd-b6f8-517e5a3e1c68)



It is not perfect but it does help to make the two data sets more comparable.

If you want to further improve the result, there are several parameters that one may consider to tune (all parameters above are either default or by gut feeling so there should be space for improvement). First of all, the FindIntegrationAnchors function chooses genes for integration based on their frequencies being identified as highly variable genes in individual data sets. Therefore, the nfeatures parameter when doing FindVariableFeatures on the two data sets definitely influence the gene set for integration. Next, since the anchoring step is the crucial step in Seurat integration, any parameter substantially affect the anchoring procedure can change the final integration. For instance, the FindIntegrationAnchors function chooses 2000 genes with the highest frequencies of being highly variable genes in individual data sets for integration by default, and this number of genes for integration can be changed by setting the anchor.features parameter in the FindIntegrationAnchors function. Similar to the issue of how many PCs to use for making tSNE/UMAP and clustering, one needs to decide which CCs to use to define cross-data-set neighbors, as set in the dims parameter. This is another parameter which can influence the result. There are more parameters which can affect in the same function, including k.anchor, k.filter and k.score, although they may not be the first parameters that you want to start with. Similarly, in thefunction IntegrateData used at the next step there is also the dims parameter, that you may want to change as well.

It is worth to mention that Seurat also provides another strategy for integrative analysis, which is data transfer. It is used when there is an existed annotated reference data, and one wants to use the reference data to assist cell type/state annotation of a new query data. The major differences between data integration and data transfer include:

Instead of generating a joint space using CCA when doing data integration, data transfer by default applies the same PCA transformation in the reference data to the query data set to identify anchors
No expression value is corrected, and therefore no joint embedding of the two data sets is created; instead, one can project cells in the query data to the reference embedding. Besides the embedding, cell labels can also be projected so that one can 'transfer' labels in the reference atlas to the query data set for annotation.
This tutorial won't cover this part as it doesn't match with the data set we have in hand. For people would like to try, it won't be difficult to follow the respective Seurat tutorial.

***`Step 2-2. Data integration using Harmony`***
Besides Seurat, there are more data integration methods available now. Harmony, developed by Soumya Raychaudhurils lab, is one of them. It is also the most highlighted integration method in the first benchmark on scRNA-seq batch effect correction tools. In brief, Harmony uses fuzzy clustering to assign every cell to multiple clusters. For each cluster, it then calculates a correction factor for each data set to move the centroid of the cluster of this data set towards the global centroid of the cluster. Since every cell is represented as a combination of multiple clusters, a cell-specific correction factor is calculated by averaging the correction factors of clusters that the cell belongs to while weighting by the cluster assignment ratio. This process will be iterated until convergence happens or reaching the iteration limits. To get more details of the method, please refer to the paper.

Harmony provides a simple API for Seurat object, which is a function called RunHarmony, so it is very easy to use. It takes the merged Seurat object (the one generated at Step 1) as the input and one needs to tell the function which metadata feature to use as the batch identity. It returns a Seurat object, with a more reduction called harmony added. It is like the corrected PCA so one should then explicitly tell Seurat to use the harmony reduction for following analysis including making UMAP embedding and identifying cell clusters.
```r
seurat <- merge(seurat_DS1, seurat_DS2) %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
library(harmony)
seurat <- RunHarmony(seurat, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.6) 

# You may also want to save the object
saveRDS(seurat, file="integrated_harmony.rds")
```


P.S. The dims.use parameter determines which dimensions (by default, of PCA) to be used for the fuzzy clustering and to be corrected. By default it uses all the calculated dimensions. The max.iter.harmony controls the maximum number of iterations to be done. By default it is 10 but since Harmony is pretty fast, it is completely fine to increase the limit so that convergence can be ensured.

We can then visualize the integration results similar to before
```r
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```



Not bad. Cells of the two samples are quite nicely mixed, and we can see some nice trajectories. Question marks may need to put at some of the mixed groups, particularly those of non-dorsal-telencephalic cells, whether or not they are indeed cells of the same cell type that should be mixed.

As you may have noticed, Harmony by default takes the PCA result as the input and iterations of correction are done to the PCs of each cell. Therefore, parameters affecting original PCA, including nfeatures in FindVariableFeatures to identify highly variable genes, should have effect on the integration. Next, when there is not specified parameter provided, the RunHarmony function takes all the available dimensions in the provided input (PCA by default). One can specify which dimensions to use by setting the dims.use parameter (this parameter is similar to the dims parameters in many Seurat functions).

***`Step 2-3. Data integration using LIGER`***
Together with Harmony and Seurat, LIGAR, developed by Evan Macosko's lab, is another data integration tool that was highlighted by the benchmark paper. It adapts integrative non-negative matrix factorization to identifying shared and dataset-specific factors for joint analysis. The detailed mathematics of the method can be found in the paper. It is implemented as the liger package in R, and it provides a wrapper for Seurat object, which relies also on the additional package SeuratWrappers in R.
```r
library(liger)
library(SeuratWrappers)

seurat <- merge(seurat_DS1, seurat_DS2) %>%
    FindVariableFeatures(nfeatures = 3000)
seurat <- ScaleData(seurat, split.by = "orig.ident", do.center = FALSE)
seurat <- RunOptimizeALS(seurat, k = 20, lambda = 5, split.by = "orig.ident")
seurat <- RunQuantileAlignSNF(seurat, split.by = "orig.ident")
seurat <- RunUMAP(seurat, dims = 1:ncol(seurat[["iNMF"]]), reduction = "iNMF")
seurat <- FindNeighbors(seurat, reduction = "iNMF", dims = 1:ncol(Embeddings(seurat, "iNMF"))) %>%
    FindClusters(resolution = 0.6)

# You may also want to save the object
saveRDS(seurat, file="integrated_liger.rds")
```


P.S. To install LIGER, do devtools::install_github('MacoskoLab/liger'). If you have a Mac machine and there is any error happened, there are some suggestions on its page. To install SeuratWrappers, do devtools::install_github('satijalab/seurat-wrappers')

Similar to above, we next visualize the integration results with the UMAP showing data sets, clusters and also some feature plots.
```r
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```



The result doesn't seem to be very easy to understand.

In case you want to improve the LIGER integration, besides the nfeatures parameter in the FindVariableFeatures function just like all the other methods, parameters in the RunOptimizeALS function also matters, such as k and lambda. LIGER has two functions called suggestK and suggestLambda to help to set these two parameters. Unfortunately these two parameters don't have their corresponding Seurat wrapper functions, or one would have to use the standalone liger package with its LIGER data type in order to use these two functions, and they are actually pretty slow. One can also change by guess with some principles, such as a larger kwould be needed when there are more sub-structure of the data; a larger lambda penalizes dataset-specific effects more strongly, so should better mixing cells from different data sets but potentially at the cost of over-integration (e.g. mixing cells with different expression signatures).

***`Step 2-4. Data integration using MNN`***
MNN, developed by John Marioni's lab in EMBL-EBI, is one of the first algorithms developed for scRNA-seq data integration or batch correction. It estimates a cell-specific correction vector based on the mutual nearest neighbors between cells from two different samples/batches to introduce correction to the dimension reduction (e.g. PCA) of the query cells. It also introduces an ordering mechanism so that it also supports integration of more than two samples/batches. Although not being the most highlighted methods in the benchmarking paper mentioned above, it is one of the best methods according to other benchmark effort (e.g. Luecken et al.). To get more details of the method, please refer to the paper. In R, the MNN algorithm is implemented in the batchelor package, and the wrapper function for a Seurat object is included in the SeuratWrappers package (RunFastMNN function).

The RunFastMNN function uses a list of Seurat objects, each of which is for one sample/batch, as the input. One can use the SplitObject function in the Seurat package to split a Seurat object given a metadata column.
```r
library(SeuratWrappers)

seurat_samples <- SplitObject(seurat, "orig.ident")
seurat_mnn <- RunFastMNN(seurat_samples)
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat),], key="MNN_")
seurat <- RunUMAP(seurat, dims = 1:20, reduction = "mnn")
seurat <- FindNeighbors(seurat, reduction = "mnn", dims = 1:20) %>%
    FindClusters(resolution = 0.6)

# You may also want to save the object
saveRDS(seurat, file="integrated_mnn.rds")
```

P.S. To install batchelor, do BiocManager::install("batchelor"). The batchelor package is required for the RunFastMNN function to work.

We can next check the the integration method via its UMAP embedding.
```r
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```



The integration looks pretty promising. In most of the time MNN performs pretty well with default parameters. Still, one can easily introduce some tuning by e.g. changing the number of features or providing a fully customized feature set for the integration. This can be done by setting up the features parameter in the RunFastMNN wrapper function. There are also more parameters that one can pass to the original function (fastMNN in the batchelor package, e.g. number of PCs to calculate).

***`Step 2-5. Data integration using RSS to BrainSpan`***
Seurat, Harmony, LIGER and MNN are probably the most commonly used methods designed for generic scRNA-seq data integration, but there are also more methods and concepts available which can be applied to data integration. One of the concept is, if there is a reference data set with multiple sample, where differences among those samples contain information of the cell type heterogeneity in the samples, representing each cell by its transcriptome similarities to those reference samples rather than its transcriptome profile itself may efficiently clean up technical noise while preserving the essential information. The method derived from this concept is called reference component analysis (RCA) or reference similarity spectrum.

To do this analysis, one firstly needs a good reference. For cerebral organoid samples, the BrainSpan bulk RNA-seq data set of human brains from early fetal development to adult by Allen Brain Atlas is a very good one.

```ref_brainspan <- readRDS("data/ext/brainspan_fetal.rds")```
Next we need to calculate similarity, or normalized Pearson's correlation between every cell and samples in the reference. There is a wrapper function for this step in the simspec package. The resulted representation is stored as one dimension reduction in the Seurat object (called rss by default). One can then use this dimension reduction for analysis including tSNE/UMAP and clustering.
```r
library(simspec)
seurat <- merge(seurat_DS1, seurat_DS2)
seurat <- ref_sim_spectrum(seurat, ref)
seurat <- RunUMAP(seurat, reduction="rss", dims = 1:ncol(Embeddings(seurat, "rss")))
seurat <- FindNeighbors(seurat, reduction = "rss", dims = 1:ncol(Embeddings(seurat, "rss"))) %>%
    FindClusters(resolution = 0.6)

plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```

P.S. If you don't have simspec package, install it via devtools::install_github("quadbiolab/simspec")



We got nice trajectories and cells from the two samples seem to mix in a reasonable way. Still, you may have realized problems when comparing the clustering results and for instance LHX9 expression.

Even if you like this result very much, there is a very obvious limitation of RCA/RSS, that there has to be a nice reference data set available so that one can calculate the similarities without lossing too much information. If your data set happened to have some interesting signals which are unavailable at all in the reference data, you would very likely miss it. As RSS represents the data purely by similarities to the reference data, if there is no change applied to the reference data, there is no much space for improving its result. The only effective parameter in the function which could be beneficial to change is the method parameter in the ref_sim_spectrum which defines the type of correlation to calculate. By default it is Pearson correlation (method = "pearson") but using Spearman correlation is also possible (method = "spearman").

***`Step 2-6. Data integration using CSS`***
At the end we would try the last data integration method in this tutorial, which is the extended version of RCA/RSS, which is cluster similarity spectrum (CSS) developed by our group. Instead of using external reference data set to represent cells in the data by similarities, it firstly does cell clustering to scRNA-seq data of each sample to be integrated, and uses the average expression profiles of the resulted clusters as the reference to calculate these similarities. More detailed description of the method can be seen in this paper.
```r
library(simspec)
seurat <- merge(seurat_DS1, seurat_DS2) %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
seurat <- cluster_sim_spectrum(seurat, label_tag = "orig.ident", cluster_resolution = 0.3)
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat, "css")))
seurat <- FindNeighbors(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css"))) %>%
    FindClusters(resolution = 0.6)

plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
```



The result doesn't seem to be worse than the others, but the trajectories look a bit odds.

Since CSS does clustering on each data set using the PCA defined when data sets were merged, nfeatures in the FindVaraibleFeatures, as well as the dims parameter in the cluster_sim_spectrum both affect the used PCs. In addition, CSS applies clustering to each data set separately, with the cluster resolution defined in the cluster_resolution parameter in the cluster_sim_spectrum function (by default cluster_resolution = 0.6). A higher resolution considers finer structure of the data which may enhance the capacity of retaining data structure but potentially at the cost of keeping more data-set-specific differences.
