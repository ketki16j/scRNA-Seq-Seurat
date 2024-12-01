***`Now starts Part 3: when you have an annotated reference data set and want it to facilitate the analysis of a new data`***

As more and more scRNA-seq data being generated all around the world, especially thanks to the effort by the Human Cell Atlas (HCA), there are more and more nicely annotated atlas-level scRNA-seq data set publicly available. It would be therefore a waste if, when analysing new but related data set, not to take the advantage of them to assist the analysis, at least to help with the annotation. This is a different scenario from the previous one, where we integrate multiple data sets with equal status. Here there is one reference data set which has been well annotated, and one query data set to be analyzed. Of course, it is still possible to use the same methods mentioned in the previous part integrating the reference and query data sets, followed by either the reference-query-joint or query-centric analysis. However, the more straightforward way is to do data transfer or projection analysis, where we fix the reference data set, and try to anchor cells or cell populations in the query data set to their counterpart in the reference.

In this part of the tutorial, we will introduce two or three strategies of how to use the reference data to help annotating the query data. We will use the DS1 which has been described above as the query data, and an annotated data subset from the same paper as the reference data. The reference data can be retrieved from this link. Among the metadata frame of the reference data, there is one column called "celltype" which shows the annotated cell type of each cell.

***`Step 0. Load data`***
We start with importing Seurat and the data we need, including the saved Seurat object of DS1 and the newly downloaded Seurat object of the reference data set.
```r
library(Seurat)
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
seurat_ref <- readRDS("ref_seurat_obj.rds")
```
Let's look at how the reference data set looks.
```r
library(patchwork)
plot1 <- UMAPPlot(seurat_ref, group.by="branch")
plot2 <- UMAPPlot(seurat_ref, group.by="celltype")
plot3 <- FeaturePlot(seurat_ref, c("SOX2","DCX","FOXG1","EMX1","DLX2","LHX9"),
                     ncol=3, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,3))
```
![image](https://github.com/user-attachments/assets/201a0fb7-d4de-4677-b031-eaafa07364e2)



We can see that the reference data set has been properly annotated, and it contains cell types representing different brain regions and neuronal subtypes.

***`Method 1-1. Transcriptome similarity on cell cluster level`***

The first strategy is very simple. We can compare the transcriptome profile of each cell population in the query data set, to the transcriptome profiles of different cell types in the reference data set. The query cell cluster can be then referred to the cell type in the reference data set which shows the highest similarity of transcriptome profiles. To do that, we need to firstly decide two things:

Based on which genes to calculate the transcriptome similarity.
How to define the similarity between two transcriptome profiles.
There are different options one can use. For the first issue, a very straightforward option is to use the highly variable genes of the reference data set. We can also intersect this gene set with the variable genes of the query data. Alternatively, we can firstly identify marker genes for each cell type in the reference data set and union them as the gene list used to represent the transcriptomic signatures.

For the second issue, one commonly used option is the correlation coefficient across genes. There are different types of correlations. The most commonly used ones include Pearson correlation and Spearman correlation. Pearson correlation focuses on the linear relationship between two vectors, while Spearman correlation is equivalent to Pearson correlation of the two ranked vector. This ranking operation allows Spearman correlation to assess also the non-linear monotonic relationship between the two vectors. Spearman correlation usually provides more robust estimate for its robustness to small amount of outliers. On the other hand, the typical ranking operation is usually time and resource consuming, especially for high-throughput data, which makes the calculation of Spearman correlation significantly slower than Pearson correlation.

In the following example, we will use the intersect of highly variable genes in the reference and query data set to calculate the Spearman correlation to represent transcriptome similarities.

First, we need to calculate the average transcriptome profiles for every annotated cell type in the reference data set and every cell cluster in the query data set.
```r
avg_expr_ref <- sapply(sort(unique(seurat_ref$celltype)), function(ct) rowMeans(seurat_ref@assays$RNA@data[,which(seurat_ref$celltype == ct)] ))
avg_expr_ds1 <- sapply(levels(seurat_DS1@active.ident), function(ct) rowMeans(seurat_ds1@assays$RNA@data[,which(seurat_ds1@active.ident == ct)]))
```
Next, we get the genes to represent transcriptome and calculate pairwise Spearman correlation across those genes' average expression between reference cell types and query clusters.
```r
genes2cor <- intersect(VariableFeatures(seurat_ref), rownames(seurat_ds1))
corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="spearman")
```

P.S. In the output matrix of the ```cor``` function, different entries (columns) of the first input matrix are represented by rows, and those of the second input matrix are represented by columns. In this case, every row in the correlation matrix is one cluster in the query data set, and every column is one cell type in the reference data set.

Now we can use a heatmap to visualize the correlation matrix.
```r
library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
```
P.S. The ```colorRampPalette``` function is used to generate a customized color palette function, which can be then used to generate a list of colors. In the example script, this function is used to generate a blue-white-red color palette and then create 30 colors along the color axis.
![image](https://github.com/user-attachments/assets/3c5ef716-f104-43e3-b095-4ce6bab98e0d)




From the heatmap we can already judge, based on the transcriptome similarity to the annotated cell types in the reference data, whether the annotation we made previously for this data set makes sense. Many of the clusters in the query data set, e.g. the query cell cluster that we annotated as "Dorsal telen. IP", indeed shows the highest similarity to the "Cortical IP" cell type in the reference data set.

***`Method 1-2. Transcriptome similarity on cell level`***

The first strategy tries to link clusters or cell types in the two data sets. While being simple, such a method also has an obvious limitation, that the clusters or cell types in the two data sets may not be defined with comparable resolution, and thus may not be comparable. This is particularly important for dynamic systems, e.g. those represent development or regeneration, where continuous cell states exist and the clustering analysis may break the continuums differently for different data sets. In that scenario, one alternative solution is thus to calculate also the transcriptome similarities to different reference cell types, but instead of doing for each query cell cluster, do it for each query cell.

Similarly, we use the intersect of highly variable genes of the two data sets as transcriptome signatures. In terms of the type of correlation to use, because of the high sparseness of the measured single-cell transcriptome profile, Spearman correlation is usually better in performance. However, as mentioned above, calculating Spearman correlation requires a ranking step. Ranking the expression of thousands of genes for thousand of cells is not only time-consuming, it also results in a huge dense matrix which needs a lot of memory, and sometimes it may be even out of the R environment capacity when the cell number is tremendous. Therefore, we shouldn't rely on the basic R function cor to do the calculation. We need a more elegant way to do the ranking while keeping the matrix sparseness, and then calculate Pearson correlation given the ranked sparse matrix using specific package designed to calculate Pearson correlation for sparse matrix.

Firstly, for the cell type average expression profiles of the reference data set, as it is not sparse and the number of entries (here the cell types) shouldn't be huge, we can use the basic rank function directly.
```r
ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
```
Next, we will introduce two ways of fast ranking of sparse matrix. The first way is to use the rank_matrix function implemented in the presto package. Yes, the package introduced before for marker identification. This function implements a fast ranking algorithm in C and is therefore very fast.

```library(presto)
ranked_expr_ds1 <- rank_matrix(seurat_DS1@assays$RNA@data[genes2cor,])$X_ranked
```
Alternatively, we can implement the sparse ranking algorithm by ourselves. You don't have to get into all the details of the function, but just copy and paste. Explaining its details is out of the scope of this tutorial. But of course, if you are interested in, it is also not too complicated to understand.
```r
rank_matrix <- function (mat) 
{
    if (is.matrix(mat) | is.data.frame(mat)) {
        ranked_mat <- apply(mat, 2, rank)
    }
    else {
        df_mat <- Matrix::summary(mat)
        dfs_mat <- split(df_mat, df_mat$j)
        df_mat_ranked <- do.call(rbind, lapply(dfs_mat, function(df) {
            num_zeros <- nrow(mat) - nrow(df)
            ranks_nonzero <- rank(df$x)
            df$x <- ranks_nonzero + num_zeros - (1 + num_zeros)/2
            return(df)
        }))
        ranked_mat <- sparseMatrix(i = df_mat_ranked$i, j = df_mat_ranked$j, 
            x = df_mat_ranked$x, dims = dim(mat), dimnames = dimnames(mat))
    }
    return(ranked_mat)
}
ranked_expr_ds1 <- rank_matrix(seurat_DS1@assays$RNA@data[genes2cor,])
```
Finally, to quickly calculate Pearson correlation between two sparse matrix or one sparse matrix and one dense matrix, the ```corSparse``` function in the ```qlcMatrix``` package is highly recommended. Afterwards, a cell in the query data set can be assigned to a reference cell type if its transcriptome similarity is the highest to that cell type.
```r
library(qlcMatrix)
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
seurat_DS1$celltype_maxcor <- ct_maxcor
```
Let's compare the annotation we did before and this projected annotation.
```r
plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(seurat_DS1, group.by="celltype_maxcor", label=T)
plot1 | plot2
```
![image](https://github.com/user-attachments/assets/d6cd2711-0ed6-4803-b7fe-b117b483fc5f)


It is not perfect, but it doesn't look bad.

We can also summarize the cell-level similarities to the query cell clusters by averaging the scaled similarities of cells in the cluster to different reference cell types, and then visualize as a heatmap similar to above.
```r
corr2ref_scaled <- scale(t(corr2ref_cell))
corr2ref_sum2cl <- t(sapply(levels(seurat_DS1@active.ident), function(cl)
  rowMeans(corr2ref_scaled[,which(seurat_DS1@active.ident == cl)]) ))
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))

```
![image](https://github.com/user-attachments/assets/b3b0f8ae-0577-42c4-8b08-55e0344d0f2c)

***`Method 2. Seurat-based label transfer`***

The above two methods are principally simple and straightforward. However, such simplicity also limits its performance. While every gene in the signature list is considered equally in those methods, those genes could have completely different power in distinguishing cell types. The projection performance would be then difficult to reach optimal if such differences on feature importances are not taken into account. Therefore, you may want to try some other more sophisticated approaches that do more than just calculating the correlations.

Here we will introduce the currently most commonly used label transfer method, which is the anchor-based label transfer implemented in Seurat. Its principle is not too complicated. It firstly applies the same dimension reduction transformation used in the reference data set (e.g. PCA) to the query data. Next, it tries to identify so-call anchors between the two data sets. Each anchor is a pair of cells, one from the reference and one from the query, which are mutual nearest neighbors with each other when calculating distances based on the transformed data. Those anchors get further filtered, requiring certain similarity on the original expression space between the two cells. Afterwards, a weights matrix is constructed, defining associations between each query cell and each anchor. This weights matrix is then used to summarize labels/values of the anchor reference cell to each query cell using the weights.

All these steps have been implemented into two functions in Seurat: ```FindTransferAnchors``` and ```TransferData```. Following are the example scripts applying to our reference and query cerebral organoid data sets.
```r
anchors <- FindTransferAnchors(reference = seurat_ref, query = seurat_DS1, dims = 1:30, npcs = 30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$celltype, dims = 1:30)
seurat_DS1$celltype_transfer <- predictions$predicted.id

plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(seurat_DS1, group.by="celltype_transfer", label=T)
plot1 | plot2
```
P.S. There are several parameters in this pipeline that could be critical to the result. By default, the FindTransferAnchors function reruns PCA on the reference data, and then applies the same transformation to the query data, using highly variable genes in the reference data set which is also detected in the query data. Besides PCA, there are other possible transformation, with CCA as the one which worths a mentioning. When doing data integration, it is CCA being used instead of PCA to maximize similarity between data sets. This is not necessary when doing data transfer, while it might distort the data too much and create artifact that may affect the performance of data transfer. For scRNA-seq data, reciprocal PCA is the third option, which not only projects the reference PCA to query but also the query PCA to reference. One can try when the default PCA doesn't work well enough. What's more, the ```weight.reduction``` parameter in the TransferData function is also important, as it defines the dimension reduction representation used to calculate weights between every query cells and the anchors. By default it uses the reference-projected PCA to define the distances, but one can also change it to other dimension reduction representation of the query data for it. Last but not least, which dimensions to use (the dims parameters) also matters.

P.S.2. If you use SCT instead of the typical logNormalize method for normalization, do set ```normalization.method="SCT"``` in the ```FindTransferAnchors``` function, and better to also set ```recompute.residuals = T```


![image](https://github.com/user-attachments/assets/3c4e3db8-4f72-4c0a-a89f-575cad0ea5dd)


The output of the ```TransferData``` function contains more information than just the predicted labels. For the categorical information to be transferred (like the cell type labels here), the output data fram also contains prediction scores of each cell to different reference cell types, which can be seen as the estimated probabilities that every query cell being different cell types. We can thus summarize those scores to the query cell clusters to compare.
```r
pred_scores_sum2cl <- t(sapply(levels(seurat_DS1@active.ident), function(cl)
  colMeans(predictions[which(seurat_DS1@active.ident == cl),-c(1,ncol(predictions))]) ))

heatmap.2(pred_scores_sum2cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = unique(seurat_ref$celltype), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))

```
![image](https://github.com/user-attachments/assets/50b665aa-29ee-4813-a2b6-9c41657e172f)


Other methods, and more to say
Above we only introduce two (or three) methods, but there are of course more. For instance, some integration methods, such as CSS and Harmony mentioned above, supports query data projection to the reference data (natively supported by CSS, and Symphony for Harmony-integrated reference). The limitation, though, is that the reference data would have to be processed using those methods. There are also deep-learning-based models which can be used to represent a given reference data set, and then being applied to other data sets for query. Examples include Sfaira developed by the Theis lab. Obviously, this also requires the reference data being processed with the framework. All those methods could work better than the introduced methods in some scenarios, while the two/three being introduced here don't have much limitation on analysis being done with the reference data set, and therefore would be the ones we usually try first.

One issue that we should always to keep in mind is that the above analysis actually assume that the reference data set is comprehensive, i.e. it contains all the cell types/states that the query data set contains. This, of course, is not always correct. Therefore, we shouldn't blindly rely on the comparison to the reference without further checking marker gene expression, and/or to use some quantitative metrics to assess how similar the projected cell types/states are to the cell populations we have in the query data.

