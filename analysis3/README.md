Now starts Part 3: when you have an annotated reference data set and want it to facilitate the analysis of a new data
As more and more scRNA-seq data being generated all around the world, especially thanks to the effort by the Human Cell Atlas (HCA), there are more and more nicely annotated atlas-level scRNA-seq data set publicly available. It would be therefore a waste if, when analysing new but related data set, not to take the advantage of them to assist the analysis, at least to help with the annotation. This is a different scenario from the previous one, where we integrate multiple data sets with equal status. Here there is one reference data set which has been well annotated, and one query data set to be analyzed. Of course, it is still possible to use the same methods mentioned in the previous part integrating the reference and query data sets, followed by either the reference-query-joint or query-centric analysis. However, the more straightforward way is to do data transfer or projection analysis, where we fix the reference data set, and try to anchor cells or cell populations in the query data set to their counterpart in the reference.

In this part of the tutorial, we will introduce two or three strategies of how to use the reference data to help annotating the query data. We will use the DS1 which has been described above as the query data, and an annotated data subset from the same paper as the reference data. The reference data can be retrieved from this link. Among the metadata frame of the reference data, there is one column called "celltype" which shows the annotated cell type of each cell.

Step 0. Load data
We start with importing Seurat and the data we need, including the saved Seurat object of DS1 and the newly downloaded Seurat object of the reference data set.

library(Seurat)
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
seurat_ref <- readRDS("ref_seurat_obj.rds")
Let's look at how the reference data set looks.

library(patchwork)
plot1 <- UMAPPlot(seurat_ref, group.by="branch")
plot2 <- UMAPPlot(seurat_ref, group.by="celltype")
plot3 <- FeaturePlot(seurat_ref, c("SOX2","DCX","FOXG1","EMX1","DLX2","LHX9"),
                     ncol=3, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,3))



We can see that the reference data set has been properly annotated, and it contains cell types representing different brain regions and neuronal subtypes.

Method 1-1. Transcriptome similarity on cell cluster level
The first strategy is very simple. We can compare the transcriptome profile of each cell population in the query data set, to the transcriptome profiles of different cell types in the reference data set. The query cell cluster can be then referred to the cell type in the reference data set which shows the highest similarity of transcriptome profiles. To do that, we need to firstly decide two things:

Based on which genes to calculate the transcriptome similarity.
How to define the similarity between two transcriptome profiles.
There are different options one can use. For the first issue, a very straightforward option is to use the highly variable genes of the reference data set. We can also intersect this gene set with the variable genes of the query data. Alternatively, we can firstly identify marker genes for each cell type in the reference data set and union them as the gene list used to represent the transcriptomic signatures.

For the second issue, one commonly used option is the correlation coefficient across genes. There are different types of correlations. The most commonly used ones include Pearson correlation and Spearman correlation. Pearson correlation focuses on the linear relationship between two vectors, while Spearman correlation is equivalent to Pearson correlation of the two ranked vector. This ranking operation allows Spearman correlation to assess also the non-linear monotonic relationship between the two vectors. Spearman correlation usually provides more robust estimate for its robustness to small amount of outliers. On the other hand, the typical ranking operation is usually time and resource consuming, especially for high-throughput data, which makes the calculation of Spearman correlation significantly slower than Pearson correlation.

In the following example, we will use the intersect of highly variable genes in the reference and query data set to calculate the Spearman correlation to represent transcriptome similarities.

First, we need to calculate the average transcriptome profiles for every annotated cell type in the reference data set and every cell cluster in the query data set.

avg_expr_ref <- sapply(sort(unique(seurat_ref$celltype)), function(ct) rowMeans(seurat_ref@assays$RNA@data[,which(seurat_ref$celltype == ct)] ))
avg_expr_ds1 <- sapply(levels(seurat_DS1@active.ident), function(ct) rowMeans(seurat_ds1@assays$RNA@data[,which(seurat_ds1@active.ident == ct)]))
Next, we get the genes to represent transcriptome and calculate pairwise Spearman correlation across those genes' average expression between reference cell types and query clusters.

genes2cor <- intersect(VariableFeatures(seurat_ref), rownames(seurat_ds1))
corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="spearman")
P.S. In the output matrix of the cor function, different entries (columns) of the first input matrix are represented by rows, and those of the second input matrix are represented by columns. In this case, every row in the correlation matrix is one cluster in the query data set, and every column is one cell type in the reference data set.

Now we can use a heatmap to visualize the correlation matrix.

library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
P.S. The colorRampPalette function is used to generate a customized color palette function, which can be then used to generate a list of colors. In the example script, this function is used to generate a blue-white-red color palette and then create 30 colors along the color axis.




From the heatmap we can already judge, based on the transcriptome similarity to the annotated cell types in the reference data, whether the annotation we made previously for this data set makes sense. Many of the clusters in the query data set, e.g. the query cell cluster that we annotated as "Dorsal telen. IP", indeed shows the highest similarity to the "Cortical IP" cell type in the reference data set.

Method 1-2. Transcriptome similarity on cell level
The first strategy tries to link clusters or cell types in the two data sets. While being simple, such a method also has an obvious limitation, that the clusters or cell types in the two data sets may not be defined with comparable resolution, and thus may not be comparable. This is particularly important for dynamic systems, e.g. those represent development or regeneration, where continuous cell states exist and the clustering analysis may break the continuums differently for different data sets. In that scenario, one alternative solution is thus to calculate also the transcriptome similarities to different reference cell types, but instead of doing for each query cell cluster, do it for each query cell.

Similarly, we use the intersect of highly variable genes of the two data sets as transcriptome signatures. In terms of the type of correlation to use, because of the high sparseness of the measured single-cell transcriptome profile, Spearman correlation is usually better in performance. However, as mentioned above, calculating Spearman correlation requires a ranking step. Ranking the expression of thousands of genes for thousand of cells is not only time-consuming, it also results in a huge dense matrix which needs a lot of memory, and sometimes it may be even out of the R environment capacity when the cell number is tremendous. Therefore, we shouldn't rely on the basic R function cor to do the calculation. We need a more elegant way to do the ranking while keeping the matrix sparseness, and then calculate Pearson correlation given the ranked sparse matrix using specific package designed to calculate Pearson correlation for sparse matrix.

Firstly, for the cell type average expression profiles of the reference data set, as it is not sparse and the number of entries (here the cell types) shouldn't be huge, we can use the basic rank function directly.

ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
Next, we will introduce two ways of fast ranking of sparse matrix. The first way is to use the rank_matrix function implemented in the presto package. Yes, the package introduced before for marker identification. This function implements a fast ranking algorithm in C and is therefore very fast.

library(presto)
ranked_expr_ds1 <- rank_matrix(seurat_DS1@assays$RNA@data[genes2cor,])$X_ranked
Alternatively, we can implement the sparse ranking algorithm by ourselves. You don't have to get into all the details of the function, but just copy and paste. Explaining its details is out of the scope of this tutorial. But of course, if you are interested in, it is also not too complicated to understand.

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
Finally, to quickly calculate Pearson correlation between two sparse matrix or one sparse matrix and one dense matrix, the corSparse function in the qlcMatrix package is highly recommended. Afterwards, a cell in the query data set can be assigned to a reference cell type if its transcriptome similarity is the highest to that cell type.

library(qlcMatrix)
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
seurat_DS1$celltype_maxcor <- ct_maxcor
Let's compare the annotation we did before and this projected annotation.

plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(seurat_DS1, group.by="celltype_maxcor", label=T)
plot1 | plot2



It is not perfect, but it doesn't look bad.

We can also summarize the cell-level similarities to the query cell clusters by averaging the scaled similarities of cells in the cluster to different reference cell types, and then visualize as a heatmap similar to above.

corr2ref_scaled <- scale(t(corr2ref_cell))
corr2ref_sum2cl <- t(sapply(levels(seurat_DS1@active.ident), function(cl)
  rowMeans(corr2ref_scaled[,which(seurat_DS1@active.ident == cl)]) ))
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))



Method 2. Seurat-based label transfer
The above two methods are principally simple and straightforward. However, such simplicity also limits its performance. While every gene in the signature list is considered equally in those methods, those genes could have completely different power in distinguishing cell types. The projection performance would be then difficult to reach optimal if such differences on feature importances are not taken into account. Therefore, you may want to try some other more sophisticated approaches that do more than just calculating the correlations.

Here we will introduce the currently most commonly used label transfer method, which is the anchor-based label transfer implemented in Seurat. Its principle is not too complicated. It firstly applies the same dimension reduction transformation used in the reference data set (e.g. PCA) to the query data. Next, it tries to identify so-call anchors between the two data sets. Each anchor is a pair of cells, one from the reference and one from the query, which are mutual nearest neighbors with each other when calculating distances based on the transformed data. Those anchors get further filtered, requiring certain similarity on the original expression space between the two cells. Afterwards, a weights matrix is constructed, defining associations between each query cell and each anchor. This weights matrix is then used to summarize labels/values of the anchor reference cell to each query cell using the weights.

All these steps have been implemented into two functions in Seurat: FindTransferAnchors and TransferData. Following are the example scripts applying to our reference and query cerebral organoid data sets.

anchors <- FindTransferAnchors(reference = seurat_ref, query = seurat_DS1, dims = 1:30, npcs = 30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$celltype, dims = 1:30)
seurat_DS1$celltype_transfer <- predictions$predicted.id

plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(seurat_DS1, group.by="celltype_transfer", label=T)
plot1 | plot2
P.S. There are several parameters in this pipeline that could be critical to the result. By default, the FindTransferAnchors function reruns PCA on the reference data, and then applies the same transformation to the query data, using highly variable genes in the reference data set which is also detected in the query data. Besides PCA, there are other possible transformation, with CCA as the one which worths a mentioning. When doing data integration, it is CCA being used instead of PCA to maximize similarity between data sets. This is not necessary when doing data transfer, while it might distort the data too much and create artifact that may affect the performance of data transfer. For scRNA-seq data, reciprocal PCA is the third option, which not only projects the reference PCA to query but also the query PCA to reference. One can try when the default PCA doesn't work well enough. What's more, the weight.reduction parameter in the TransferData function is also important, as it defines the dimension reduction representation used to calculate weights between every query cells and the anchors. By default it uses the reference-projected PCA to define the distances, but one can also change it to other dimension reduction representation of the query data for it. Last but not least, which dimensions to use (the dims parameters) also matters.

P.S.2. If you use SCT instead of the typical logNormalize method for normalization, do set normalization.method="SCT" in the FindTransferAnchors function, and better to also set recompute.residuals = T




The output of the TransferData function contains more information than just the predicted labels. For the categorical information to be transferred (like the cell type labels here), the output data fram also contains prediction scores of each cell to different reference cell types, which can be seen as the estimated probabilities that every query cell being different cell types. We can thus summarize those scores to the query cell clusters to compare.

pred_scores_sum2cl <- t(sapply(levels(seurat_DS1@active.ident), function(cl)
  colMeans(predictions[which(seurat_DS1@active.ident == cl),-c(1,ncol(predictions))]) ))

heatmap.2(pred_scores_sum2cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = unique(seurat_ref$celltype), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))



Other methods, and more to say
Above we only introduce two (or three) methods, but there are of course more. For instance, some integration methods, such as CSS and Harmony mentioned above, supports query data projection to the reference data (natively supported by CSS, and Symphony for Harmony-integrated reference). The limitation, though, is that the reference data would have to be processed using those methods. There are also deep-learning-based models which can be used to represent a given reference data set, and then being applied to other data sets for query. Examples include Sfaira developed by the Theis lab. Obviously, this also requires the reference data being processed with the framework. All those methods could work better than the introduced methods in some scenarios, while the two/three being introduced here don't have much limitation on analysis being done with the reference data set, and therefore would be the ones we usually try first.

One issue that we should always to keep in mind is that the above analysis actually assume that the reference data set is comprehensive, i.e. it contains all the cell types/states that the query data set contains. This, of course, is not always correct. Therefore, we shouldn't blindly rely on the comparison to the reference without further checking marker gene expression, and/or to use some quantitative metrics to assess how similar the projected cell types/states are to the cell populations we have in the query data.

Now starts Part 4: more optional advanced analysis for scRNA-seq data
The analysis mentioned above are mostly about scRNA-seq data preprocessing (e.g. normalization, dimension reduction and data integration) as well as the most basic analysis (e.g. clustering and marker identification). Depending on the systems that the scRNA-seq data represents, more analysis can be potentially applied to investigate the relevant biological insight. These analysis include but not limit to pseudotime analysis (which has been mentioned above), differential expression analysis between conditions, RNA velocity analysis, branching point analysis, cell-cell communication analysis with ligand-receptor pairing, and gene regulatory network inferences. In the following section, we will briefly introduce some of those advanced analysis on scRNA-seq data.

Part 4-1. Cluster connectivity analysis with PAGA
The clustering analysis described above is the most commonly used way to summarize the cell type/state heterogeneity in the data. However, the basic clustering analysis does not provide the information about how each cluster of cells may have connected with each other. This is very important for understanding dynamical systems related to development and regeneration, for instance. There are different ways to complement this issue. One option is to do branching analysis, which instead of defining clusters, describes the cell type/state landscape as a tree, with the very initial cell state as the root and every time when two different cell types/states are derived from the same ancestor state, a branching point is defined. There are quite some tools available for this branching analysis, with the most commonly used ones including monocle/monocle2/monocle3, URD, and Slingshot. Actually, the branch analysis is usually coupled with pseudotime reconstruction, both as parts of the trajectory analysis. Therefore, those tools for branching analysis usually contain the function to estimate pseudotimes; and many tools developed to generate pseudotimes, e.g. diffusion map as described above, also include the function to identify branching points.

Besides the branching analysis, another strategy is to rely on the clustering results, and apply certain statistcs to evaluate the strength of connectivity between every two clusters. In most of the time, this connectivity is defined by how likely cells in one cluster being one of the nearest neighbors of a cell in another cluster. The most commonly used tool for this analysis is PAGA developed by Fabian Theis lab in Helmholtz Center Munich, which is included as a part of the scanpy toolkit in Python.

Next, we will use DS1 as the example to show how to run PAGA using the scanpy package in R, with the help of the Python interface provided by the reticulate package.

First, we should load the DS1 Seurat object with the above preprocessing being done.

library(Seurat)
library(Matrix)
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
As scanpy is a python package and doesn't support a Seurat object as the input, we need to store its information in a format that scanpy supports. Possible format include h5ad and loom. Following is the example we use the loomR package to create a loom file with the information needed (e.g. PCA, UMAP, and cell type annotation).

library(loomR)
cell_attrs <- list(pca = Embeddings(seurat_DS1,"pca")[,1:20],
                   umap = Embeddings(seurat_DS1,"umap"),
                   celltype = seurat_DS1@active.ident)
loom <- loomR::create("DS1/loom_obj.loom",
                      data = seurat_DS1[['RNA']]@data,
                      layers = list(counts = seurat[['RNA']]@counts),
                      cell.attrs = cell_attrs)
loom$close_all()
P.S. To install loomR, one can use devtools as following devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

These operations generate a new file (DS1/loom_obj.loom) which can be then used as the input to run scanpy.

Alternatively, one can use the anndata package to create a h5ad file with similar information. Now this is actually more recommended, as anndata is more accepted as the default file format for many scRNA-seq analysis tools in Python, such as scanpy, scvelo, and cellrank.

library(anndata)
adata <- AnnData(X = t(seurat_DS1[['RNA']]@data),
                 obs = data.frame(celltype = seurat_DS1@active.ident, row.names = colnames(seurat_DS1)),
                 var = seurat_DS1[['RNA']]@meta.features,
                 layers = list(counts = t(seurat_DS1[['RNA']]@counts)),
                 obsm = list(pca = Embeddings(seurat_DS1,"pca")[,1:20],
                             umap = Embeddings(seurat_DS1,"umap"))
                )
adata$write_h5ad("DS1/anndata_obj.h5ad")
P.S. Installing anndata can be easily done via install.packages("anndata")

These generate the new file DS1/anndata_obj.h5ad which can also be used as the input to run scanpy.

Next, one has the option to then swtich to Python (>=3.6) with scanpy installed to do the next steps. If you want to stay in R, you need the R package called reticulate. This package, developed by RStudio, provides the R interface to Python, so that one can easily run Python scripts in the R environment. In principle, it should have been installed when the Seurat package was installed. One can make sure by explicitly installing it with install.packages(reticulate). One can also install the develop branch of reticulate with remotes::install_github("rstudio/reticulate").

Please be aware that reticulate only provides the interface to Python, not the Python itself. Therefore, one still need to install a Python (>=3.6) which can be called by reticulate. Next, we import the reticulate package, install the scanpy package, and import it to the R environment.

library(reticulate)
py_install("scanpy", pip=T)
sc <- import("scanpy")
If you don't see any error, you have scanpy successfully installed in your Python and also have it imported. As mentioned above, scanpy is a comprehensive toolkit for scRNA-seq data analysis, similar to Seurat but implemented in Python. It therefore contains a lot more functionalities than just PAGA. Here, we just use the PAGA analysis in the package, and rely on the Seurat analysis result to provide the information needed.

adata_DS1 <- sc$read_loom("DS1/loom_obj.loom") # option1: load the loom file
adata_DS1 <- sc$read("DS1/anndata_obj.h5ad") # option2: load the h5ad file

sc$pp$neighbors(adata_DS1, n_neighbors=20L, use_rep='pca')
sc$tl$paga(adata_DS1, groups='celltype')
adata_DS1$write_h5ad("DS1/anndata_obj.h5ad")
The rationale of PAGA is to firstly construct a k-nearest neighbor network of cells using the provided dimension reduction matrix (in the example it is 'pca'), and then count the number of edges connecting two cell clusters. This number will be then compared with the expected edge number in a random network to estimate the PAGA connectivity between the two clusters. The resulted connectivity is a value between 0 and 1, with 0 meaning no connectivity and 1 meaning full connectivity. It is therefore straightforward to see that critical parameters that could have affected the PAGA results, besides the clustering result, include 1) the dimension reduction used to identify cell neighbors (the use_rep parameter in the neighbors function) and 2) the number of neighbors per cell (the n_neighbors parameter). It is also worth mentioning that there are additional parameters in the neighbors function, e.g. the method parameter to determine the method to quantify cell connectivity, and the metric parameter to determine the type of distance being calculated between cells, which could lead to changes of the neighboring network, and therefore the PAGA result.

Next, we can plot the PAGA estimated cluster connectivity.

plt <- import("matplotlib")
plt$use("Agg", force = TRUE)
sc$pl$paga(adata_DS1,
           color='celltype',
           fontsize=7,
           frameon=FALSE,
           save="DS1_paga.png")
This generates a PNG image in the figures subfolder (figures/DS1_paga.png) 

It is not perfect, but one can see that the summarized graph on the cluster level somehow recapitulates the general heterogeneity of the data set. Many of those thickest edges with the strongest connectivity between two clusters, indeed represent the differentiation process, e.g. the strong connectivity between dorsal telen. IP and dorsal telen. neurons.

P.S. When running with reticulate, it is possible to encounter errors like ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found. That is due to the C version conflict between the one in the system that being used in R, and the one used in the conda environment that reticulate thought to be used. In that case, running PAGA in Python directly would be the easiest solution.

Part 4-2. Pseudotime reconstruction without subseting into an unbranched trajectory
In a previous section we have shown how to reconstruct pseudotimes for cells at an unbranched trajectory using diffusion map. However, when the system is complex and includes lots of different terminal, as well as initial/intermediate states which are shared by trajectories of multiple terminal states, finding cells representing a single unbranched trajectory becomes a very difficult task, and by itself becomes a challenge that we would have to develop or apply specific tools to solve. Interestingly, many of those tools rely on random walk across the neighborhood graph of cells and strongly recommend to include pseudotimes as one input. Obviously, this makes the pseudotime analysis and trajectory analysis into a loop. There are multiple approaches to solve this issue. For example, we can use other analysis than the transcriptome similarity based pseudotimes to provide clues about cell state progression, and one of the most pronounced methods in this category is the RNA velocity analysis which will be introduced in the next section. Alternatively, we can try to perform pseudotime analysis without subseting into a single unbranched trajectory.

For example, we can directly apply diffusion map and diffusion pseudotime (DPT), which have been introduced earlier, to the complete DS1 data rather than just the cortical portion.

library(Seurat)
library(destiny)

seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
dm <- DiffusionMap(Embeddings(seurat_DS1, "pca")[,1:20])
dpt <- DPT(dm)
seurat_DS1$dpt <- rank(dpt$dpt)

FeaturePlot(seurat_DS1, c("dpt","SOX2","NHLH1","DCX"), ncol=4)


It doesn't look bad! Although one obvious issue is that the reconstructed pseudotimes seem to be flipped, with the SOX2+ progenitors having higher pseudotime than the DCX+ neurons. The easy way to fix it of course is to manually flip it.

seurat_DS1$dpt <- max(seurat_DS1$dpt) - seurat_DS1$dpt
FeaturePlot(seurat_DS1, c("dpt","SOX2","NHLH1","DCX"), ncol=4)


Alternatively, we can specify the 'tip' cell for the random walk to start. Let's randomly pick three progenitor cells as the tips, and then run DPT again.

set.seed(12345)
idx <- sample(which(seurat@active.ident %in% c('Dorsal telen. NPC',
                                               'G2M dorsal telen. NPC',
                                               'Dien. and midbrain NPC',
                                               'G2M Dien. and midbrain NPC')),3)
dpt2 <- DPT(dm, tips=idx)
seurat_DS1$dpt2 <- rank(dpt2$dpt)

FeaturePlot(seurat_DS1, c("dpt","dpt2"), ncol=2)


The problem of this approach is that as the NPC clusters actually represent a range of cell states along the progenitor to neuron differentiation, the randomly picked tip cells might not be at the actual beginning of the trajectory but somewhere in the middle. To further solve this issue, we can try to firstly use the random_root function in the destiny package explicitly, which firstly randomly pick a cell and then identify the cell with the largest DPT from the selected cell, to identify tip candidates (this is also how the default DPT function works). Next, we can subset the candidates with the cells annotated as progenitors, and use the intersect as the tip cells for DPT.

tips_cand <- sapply(1:100, function(i){ random_root(dm) })
idx_NPC <- which(seurat@active.ident %in% c('Dorsal telen. NPC',
                                            'G2M dorsal telen. NPC',
                                            'Dien. and midbrain NPC',
                                            'G2M Dien. and midbrain NPC'))
tips_cand <- as.numeric(names(which.max(table(tips_cand[tips_cand %in% idx_NPC]))))
dpt3 <- DPT(dm, tips=tips_cand)
seurat_DS1$dpt3 <- rank(dpt3$dpt)

FeaturePlot(seurat_DS1, c("dpt","dpt2", "dpt3"), ncol=3)



And of course, there are a lot of other approaches and algorithms developed in the past years to reconstruct pseudotimes. In R, the famous options besides the diffusion pseudotime used in the above examples include Monocle developed by the Trapnell lab, Slingshot by the Dudoit lab, and URD by the Regev lab and Schier lab. There are also more options available in Python, which can also been used in R by using reticulate, similar to how we do the PAGA analysis. One of those Python options is Palantir developed by the Pe'er lab. Besides, diffusion pseudotime (implemented in Scanpy) and Slingshot (pyslingshot) are also available in Python. Many of those methods not only provide pseudotime reconstruction, but also trajectory analysis inspired by the calculated pseudotime.

Part 4-3. RNA velocity analysis
RNA velocity analysis was firstly proposed by La Manno et al. in Sten Linnarsson lab in Karolinska institute and Peter Kharchenko lab in Harvard Medical School in 2018. It is an analysis based on a simple model of transcriptional dynamics. In this model, the transcript number of a certain gene that one can measure in a cell is determined by several processes: transcription, splicing and degradation. Considering that the current mature RNAs represent the current state of the cell, such a state may stay steady if the dynamics of RNA degradation and transcription+splicing reach equilibrium, or it may change over time. Assuming the cell state is not steady, the time lag between transcription and RNA processing (e.g. splicing) make it possible to infer how the cell state is going to change, if the degradation rates of different transcripts are known. Based on this concept, La Manno et al. developed the first algorithm, to use the exonic reads as the proxy of mature RNA transcripts, and the intronic reads as the proxy of the immature RNA transcripts to be spliced. The details of the method can be found in the published paper. The R implementation of this method is available in the velocity.R package.

Based on their work, Bergen et al. in Alexander Wolf lab and Fabian Theis lab in Helmholtz Center Munich further generalized the transcriptional dynamics model estimation procedure, so that it no longer relies on the assumption of steady cell states. The description of the method can be found in the paper. They also developed the python package scvelo, which is not only methodologically more general, but also computationally more efficient.

Next we will use DS1 as the example to show how to apply RNA velocity analysis using the scvelo package in R, similar to above with PAGA.

As RNA velocity analysis requires the exonic and intronic count matrices separately for the cells, these two matrices need to be generated. Unfortunately, for 10x Genomics scRNA-seq platform which is the most commonly used one right now, its routine preprocessing pipeline CellRanger only counts the transcript number per cell for those overlapping with exonic regions. Therefore, one needs to do extra counting to generate the matrices in need. There are in general two strategies:

Use the CellRanger mapping result (in BAM format) and the gene annotation, generate count matrices with transcript counting software (e.g. dropEst).
Use pseudomapping to generate exonic and intronic count matrices, with tools like kallisto.
Here, the example matrices of DS1 was generated with dropEst.

Firstly, let's load the DS1 Seurat object with the above preprocessing being done, as well as the exonic and intronic count matrices in R. The dropEst-generated count matrices include all detected cellular barcodes, so we shall subset only cells included in the Seurat object.

library(Seurat)
library(Matrix)

seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
mats <- readRDS("DS1/mats_dropest.rds")
mats <- lapply(mats, function(mat) mat[,colnames(seurat_DS1)])
Next, we need to create the loom or h5ad file which contains the information we need for running scvelo. This is very similar to above for PAGA, but the spliced and unspliced data matrix would have to be used.

The following is how to create the loom file:

library(loomR)
cell_attrs <- list(pca = Embeddings(seurat_DS1,"pca")[,1:20],
                   umap = Embeddings(seurat_DS1,"umap"),
                   celltype = seurat_DS1@active.ident)
shared_genes <- intersect(rownames(mats$exon),rownames(mats$intron))
loom <- loomR::create("DS1/loom_obj_scvelo.loom",
                      data = mats$exon[shared_genes,],
                      layers = list(spliced = mats$exon[shared_genes,],
                                    unspliced = mats$intron[shared_genes,]),
                      cell.attrs = cell_attrs)
loom$close_all()
Following is the way to create the h5ad file:

library(anndata)
shared_genes <- intersect(rownames(mats$exon),rownames(mats$intron))
adata <- AnnData(X = t(mats$exon[shared_genes,]),
                 obs = data.frame(seurat_DS1@meta.data, celltype=seurat_DS1@active.ident),
                 var = NULL,
                 layers = list(spliced = t(mats$exon[shared_genes,]),
                               unspliced = t(mats$intron[shared_genes,])),
                 obsm = list(pca = Embeddings(seurat_DS1,"pca")[,1:20],
                             umap = Embeddings(seurat_DS1,"umap"))
                )
adata$write_h5ad("DS1/anndata_obj_scvelo.h5ad")
And then, similar to when doing PAGA analysis, one can swtich to Python (>=3.6) with scvelo installed to do the next steps, or to use the reticulate R package.

Next, we import the reticulate package, install the scvelo package, and import it to the R environment.

library(reticulate)
py_install("scvelo", pip=T)
scvelo <- import("scvelo")
If you don't see any error, you have scvelo successfully installed in your Python and also have it imported. Next, it's time to run the RNA velocity.

adata_DS1 <- scvelo$read_loom("DS1/loom_obj_scvelo.loom") # option1: load the loom file
adata_DS1 <- scvelo$read("DS1/anndata_obj_scvelo.h5ad") # option2: load the h5ad file

scvelo$pp$filter_and_normalize(adata_DS1,
                               min_shared_counts=as.integer(10),
                               n_top_genes=as.integer(3000))
scvelo$pp$moments(adata_DS1,
                  n_neighbors = as.integer(30),
                  use_rep = "pca")
scvelo$tl$velocity(adata_DS1)
scvelo$tl$velocity_graph(adata_DS1)
P.S. You may have noticed the as.integer function. This is used because Python is type-sensitive and the called Python function expects an integer rather than a float number. However, numbers are considered as double float type by default in R, so it would return errors if the values are not converted to integer explicitly with as.integer.

There are parameters that one can change and tune in scvelo. Also with velocity estimated, one can do more analysis with scvelo, e.g. velocity pseudotime estimation; but these are out of the scope of this tutorial. To get more details, please visit the scvelo manual page (https://scvelo.readthedocs.io/index.html). At the end of this part, let's visualize the velocity estimates.

plt <- import("matplotlib")
plt$use("Agg", force = TRUE)
scvelo$pl$velocity_embedding_stream(adata_DS1,
                                    basis="umap",
                                    color="celltype",
                                    dpi=120,
                                    figsize = c(8,8),
                                    save="DS1_scvelo_stream.png")
This generates a PNG image in the figures subfolder (figures/DS1_scvelo_stream.png) 

See the nice arrows! They point out the estimated cell state transitions and one can clearly see the transition from NPC to neurons, which indicates the neuronal differentiation and maturation processes.

Besides, scvelo also allows to generate pseudotimes based on the estimated RNA velocity. There are two possible velocity pseudotimes supported by scvelo. One is called "velocity pseudotime" which is estimated using a similar way as the diffusion pseudotime introduced above. The difference is that when doing random walk, it is the symmetric transition matrix that is estimated based on diffusion map embedding being used in the diffusion pseudotime estimation, while it is the directional transition matrix estimated by RNA velocity being used in the velocity pseudotime estimation. The other is called "latent pseudotime" which is fully based on the velocity dynamics.

scvelo$tl$velocity_pseudotime(adata_DS1)
scvelo$tl$recover_dynamics(adata_DS1)
scvelo$tl$recover_latent_time(adata_DS1)
P.S. It is required to have the full splicing dynamics being recovered to estimate the latent time (the recover_dynamics function). However, this step is very slow and has to go through every gene. One has to prepare for that.

Afterwards, we can extract the estimated pseudotimes, and project them on the UMAP as feature plots for visualization and comparison.

seurat_DS1$velocity_pseudotime <- adata_DS1[['obs']]$velocity_pseudotime
seurat_DS1$latent_time <- adata_DS1[['obs']]$latent_time

FeaturePlot(seurat_DS1,
            c("velocity_pseudotime", "latent_time")) & NoAxes()



Here it is quite obvious that the velocity pseudotime can represent the NPC-to-neuron differentiation trajectory much better than the latent time. This often happens but also not always. There are also data sets where the latent time shows better the trajectory. One needs to check and compare.

P.S. You may have expected some interactions between scvelo and scanpy, as they are developed by the same group. And you are right! For instance, in the PAGA function there is one parameter called use_rna_velocity, which is False by default. However, if your anndata object contains the RNA velocity results by scvelo, you can set this parameter to True, and PAGA will then run based on the directed cell graph estimated by the RNA velocity analysis.

Last, let's again save the processed AnnData and the Seurat object with the velocity-based pseudotime information

adata_DS1$write_h5ad('DS1/anndata_obj_scvelo.h5ad')
saveRDS(seurat_DS1, file='DS1/seurat_obj_all.rds')
Part 4-4. Trajectory analysis with CellRank
The pseudotime analysis and RNA velocity analysis give us very rich information about cell state dynamics on the single cell level. And when we have a complex and dynamic system like developing primary tissues or stem cell models like organoids, we would definitely want to use those information to characterize how cells move from the initial stem cell state through different intermediate cell states and eventrally reach different terminal cell types. PAGA, which has been introduced earlier, is somehow on this line, but it is performed on the cluster level and somehow lack of flexibility of which information to use (either similarity-based or velocity-based connectivities). Other trajectory analysis tools such as Monocle, Slingshot and URD can estimate how the cell state trajectories branch, and assign every cell to one partition of the branched trajectory. However, those methods mostly rely on only the similarity-based pseudotime and it is not very easy to hack them in order to use also the RNA velocity information. Considering limitations of those approaches, people in the field are therefore eagerly looking for a method which can perform analysis on the single-cell level and have the flexibility to consider different information, including transcriptome similarity, pseudotime, RNA velocity, as well as other meta information (e.g. time points and metabolic labeling) and analysis (e.g. CytoTRACE) which can inform the system dynamics, in a customizable manner. Further considering the growing amount of single-cell data to analyze, the scalability of a computational method is also getting more and more important. CellRank, a Python-based framework to study cellular dynamics based on Markov state modeling of single-cell data, provides a great balance and functionalities on those topics.

Similar to other Python-based tools we have introduced above, it is possible to import CellRank via reticulate in R. The standard procedure of CellRank is also not very complicated. There are comprehensive tutorials in the CellRank webpage that guide you through the whole process.

Here in our example, we can start with the H5AD file we saved after the scVelo analysis, after making sure that CellRank is successfully installed (with py_install('cellrank', pip=TRUE) in R).

library(reticulate)
sc <- import('scanpy')
cr <- import('cellrank')
adata_DS1 <- sc$read_h5ad('DS1_scvelo.h5ad')
Next, we can initialize different kernels based on pseudotime, transcriptomic similarity (PCA), as well as scVelo results. In CellRank, kernels are used to compute cell-cell transition matrices, which are later analyzed. Here with what we have we can construct three different kernels: the pseudotime kernel which can calculate transition probabilities between cells based on the provided pseudotime, for example, the diffusion pseudotime we reconstructed earlier; the connectivity kernel which do the calculation based on transcriptomic similarities between cells (usually estimated as Euclidean distance in the dimensionality reduction space like PCA or any integrated latent representations); and the velocity kernel which calculates the transition probabilities based on the RNA velocity analysis we did above.

pk <- cr$kernels$PseudotimeKernel(adata_DS1, time_key="dpt3")$compute_transition_matrix()
ck <- cr$kernels$ConnectivityKernel(adata_DS1)$compute_transition_matrix()
vk <- cr$kernels$VelocityKernel(adata_DS1)$compute_transition_matrix()
We can further create combined kernels based those kernels with different customized weights. Unfortunately from this step it is a bit problematic to run the code in R. Luckily, reticulate provides the option (repl_python) to open a Python interactive environment directly in R, from where we can still access objects we make in the R environment (via r.object). Also, any object we create in that Python environment can also be accessed later when we are back to the R (via py$object). Let's first open the Python interactive session.

repl_python()
Now we can continue with coding in the opened Python session.

import cellrank as cr
import numpy as np

combined_kernel = 0.5 * r.vk + 0.3 * r.pk + 0.2 * r.ck
Now we have the combined kernel created. It is a combination of the three kernels, with the velocity kernel contributing to 50% of the final kernel, and the pseudotime kernel and connectivity kernel contributing 30% and 20% respectively. Next we can try to infer terminal states based on the combined kernel.

g = cr.estimators.GPCCA(combined_kernel)
g.fit(n_states=15, cluster_key="celltype")
g.predict_terminal_states(method="top_n", n_states=6)
g.plot_macrostates(which="terminal")


So actually it doesn't look too bad! Different types of neurons are successfully inferred as terminal states. On the other hand, the result is also not perfect, for example, there is a group of NPCs considered as potential terminal states. Of course, one should in principle look into that population and see whether it actually represents any unknown terminal cell states in the system, which could lead to exciting findings. But let's assume this is an artifact, then we can also choose to manually set the terminal states as the combined randomly subsets (30 cells in the following example) of different annotated neuronal cell types.

g = cr.estimators.GPCCA(combined_kernel)
g.set_terminal_states({"Midbrain-hindbrain boundary neuron": r.adata_DS1[r.adata_DS1.obs['celltype'] == "Midbrain-hindbrain boundary neuron"].obs_names[:30],
                       "Ventral telen. neuron": r.adata_DS1[r.adata_DS1.obs['celltype'] == "Ventral telen. neuron"].obs_names[:30],
                       "Dorsal telen. neuron": r.adata_DS1[r.adata_DS1.obs['celltype'] == "Dorsal telen. neuron"].obs_names[:30],
                       "Dien. and midbrain excitatory neuron": r.adata_DS1[r.adata_DS1.obs['celltype'] == "Dien. and midbrain excitatory neuron"].obs_names[:30],
                       "Dien. and midbrain inhibitory neuron": r.adata_DS1[r.adata_DS1.obs['celltype'] == "Dien. and midbrain inhibitory neuron"].obs_names[:30]})
g.plot_macrostates(which="terminal")


P.S. It is also possible to do even more flexible terminal state specification. For example, one can firstly run the data-driven terminal state prediction, and then extract the inferred states via g.terminal_states and manipulate the returned Series object (for example, to exclude some labels that are not supposed to be terminal states), and then assign the manipulated result back to the estimator using the g.set_terminal_states function.

Now we can start to calculate for every cell the fate probabilities of how likely it will in the end turn into each of the terminal states.

g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="right", basis='umap', same_plot=False)


We can then extract the fate probability estimates, save them into a data frame.

import pandas as pd
prob = pd.DataFrame(g.fate_probabilities).set_index(g.terminal_states.index).set_axis(g.terminal_states.cat.categories, axis=1)
And now, we can quit the interactive Python session with exit to go back to R, and save the fate probability data frame into the Seurat object

exit
library(Seurat)
seurat_DS1 <- readRDS(file='DS1/seurat_obj_all.rds')
seurat_DS1@meta.data[,colnames(py$prob)] <- py$prob[colnames(seurat_DS1),]

FeaturePlot(seurat_DS1, colnames(py$prob), ncol=5)


With the fate probabilities estimated on the cell level, one can potentially do more analysis, for example, to identify potential driver genes which correlate to changes of probability into each cell fate. One can try to place the cells into different part of the branched differentiation trajectory based on the estimated probability. For example, in our previous papers studying cell fate specification in brain organoids (Fleck, et cl. 2023) and retinal organoids (Wahle, et al. 2023), we summarized the single-cell level fate probabilities into meta-clusters and visualize how different cell types in the system are gradually specficied.

Part 4-5. Cell communication analysis
The above analysis focus mostly on cell states of single cells. In biology, what can be equally or even more important is communications among different cells. Unfortunately, such communications cannot be directly measured by scRNA-seq data. However, as the communications usually rely on receptor proteins and ligand molecules that specifically bind to its corresponding receptor, given a list of the paired ligand and receptor, it is possible to infer the existence of such communications, assuming that cells/cell types that communicate with each other co-express the interacting ligands and receptors. This is the principle of most of the existed tools to investigate putative cell-cell communications. Among those tools, CellPhoneDB, developed by Sarah Teichmann's lab in Wellcome Sanger Institute, is one of the most commonly used one. More details of the method can also been found in the paper.

In this part of the tutorial, we will use DS4 as the example to show how to infer communications between cell types using scRNA-seq data and cellphonedb. Just to mention, DS4 is not about cerebral organoid, but a subset of developing human duodenum scRNA-seq data presented in this paper. It is suggested that the interactions between epithelial and mesenchymal populations are critical for gut development. Therefore, it would be interesting to see whether we can infer such interaction from the scRNA-seq data and identify ligand-receptor pairs contributing to it.

First of all, we need to have CellPhoneDB installed. It is a Python package and one can install it following the tutorial in its github page. If you are a conda user, you can also use conda to manage the environment. For example,

conda create -n env_cellphonedb python=3.7
conda activate env_cellphonedb
pip install cellphonedb
Next, we need to generate the input files for CellPhoneDB, which include:

A TSV table of gene expression in single cells (rows as genes, columns as cells)
A TSV table with two columns, the first column being cell IDs (the same as the column names in the expression table), the second column being cell type annotation
Let's load the data into Seurat and take a look at it. Please note that the cell annotation is already included in the provided metadata.

library(Seurat)
library(Matrix)
library(dplyr)

counts <- readMM("DS4/matrix.mtx.gz")
metadata <- read.table("DS4/metadata.tsv.gz")
features <- read.table("DS4/features.tsv.gz")[,1]
dimnames(counts) <- list(features, rownames(metadata))
seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)

seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims=1:20)

UMAPPlot(seurat_int, group.by="Cell_type", label=T) & NoAxes() & NoLegend()



Now it's time to generate the tables that CellPhoneDB requires

expr_mat <- as.matrix(seurat@assays$RNA@data)
write.table(expr_mat, file="DS4/expr.txt", quote=F, sep="\t")
write.table(data.frame(Cell_bc = colnames(seurat), Cell_type = seurat$Cell_type),
            file="DS4/metadata.txt", row.names=F, quote=F, sep="\t")
P.S. CellPhoneDB requires the densed expression table being stored in a text file, and this becomes impractical when the scRNA-seq data is big. In that case, one needs to decrease the number of cells for CellPhoneDB to read by e.g. subseting the data. This is not necessary in this example because the data has already been subset.

Next, we move back to the shell to run CellPhoneDB

# go to the folder with the generated expr.txt and metadata.txt files
cd DS4

# make sure to activate the cellphonedb conda environment, if conda is used
conda activate env_cellphonedb

# run cellphonedb
cellphonedb method statistical_analysis --counts-data gene_name metadata.txt expr.txt
P.S. if you are using the newest numpy package, running CellPhoneDB may return the error of AttributeError: type object 'object' has no attribute 'dtype'. This is because the pandas version which cellphonedb depends on does not work with the latest numpy version (see this page Teichlab/cellphonedb#266). As mentioned in the same page, there are at least two solutions: 1. if the cellphonedb is installed in a conda environment, installing tensorflow with conda install tensorflow can make sure that a compatible numpy version is installed; 2. downgrade the numpy version pip install --force-reinstall numpy==1.19.5

CellPhoneDB estimates significance of ligand-receptor interactions between cell types using permutation of cell type labels (by default 1000 times). Therefore, it takes quite some time to get the final results. One needs to be patient here. When it is done, a folder called out by default is created with the cellphonedb output inside. It should include four files:

deconvoluted.txt
means.txt
pvalues.txt
significant_means.txt
CellPhoneDB has two plotting functions implemented (dotplot and heatmap) that one can try. Alternatively, one can use R to do the plotting. The following example is to plot the number of inferred interacting pairs between every two cell types. Columns represent cell types secreting ligands, while rows represent cell types receiving the signals.

library(gplots)

p <- read.csv("DS4_cellphonedb/out/pvalues.txt", header=T, sep="\t", check.names=F)
num_pairs <- colSums(p[,-(1:11)] < 0.01)
num_pairs <- data.frame(partner1 = sapply(strsplit(names(num_pairs),"\\|"),"[",1),
                        partner2 = sapply(strsplit(names(num_pairs),"\\|"),"[",2),
                        num = num_pairs)
mat_num_pairs <- sapply(sort(unique(num_pairs$partner1)), function(x)
  sapply(sort(unique(num_pairs$partner2)), function(y)
    num_pairs$num[which(num_pairs$partner1 == x & num_pairs$partner2 == y)]))

bluered_colscheme <- colorRampPalette(c("#4575B4","#9DCEE3","#EFE9C3","#FA9D59","#D73027"))
heatmap.2(mat_num_pairs + t(mat_num_pairs) - diag(diag(mat_num_pairs)),
          trace="none", scale="none", col = bluered_colscheme(30), key=F, keysize=0.8, margins = c(9,9))



One can also separate source and target based on the information in the CellPhoneDB output

p <- p[p$secreted=="True" &
         ((p$receptor_a == "True" & p$receptor_b == "False") |
          (p$receptor_a == "False" & p$receptor_b == "True")),]

idx <- which(p$receptor_a == "False")
num_pairs <- colSums(p[idx,-(1:11)] < 0.05)
num_pairs <- data.frame(from = sapply(strsplit(names(num_pairs),"\\|"),"[",1),
                        to = sapply(strsplit(names(num_pairs),"\\|"),"[",2),
                        num = num_pairs)
idx <- which(p$receptor_a == "True")
num_pairs_2 <- colSums(p[idx,-(1:11)] < 0.05)
num_pairs_2 <- data.frame(from = sapply(strsplit(names(num_pairs_2),"\\|"),"[",2),
                          to = sapply(strsplit(names(num_pairs_2),"\\|"),"[",1),
                          num = num_pairs_2)
num_pairs$num <- num_pairs$num + num_pairs_2$num
mat_num_pairs <- sapply(sort(unique(num_pairs$from)), function(x)
  sapply(sort(unique(num_pairs$to)), function(y)
    num_pairs$num[which(num_pairs$from == x & num_pairs$to == y)]))

bluered_colscheme <- colorRampPalette(c("#4575B4","#9DCEE3","#EFE9C3","#FA9D59","#D73027"))
heatmap.2(mat_num_pairs,
          trace="none", scale="none", col = bluered_colscheme(30), key=F, keysize=0.8, margins = c(9,9),
          xlab="FROM", ylab="TO")



More can be done with the CellPhoneDB output. Here we are going to try another tool called COMUNET, which is a R package developed by Antonio Scialdone's lab in Helmholtz Zentrum Munich. It provides additional statistics (e.g. clustering of communications) given the CellPhoneDB output. More details can be found in its paper, and more tutorials are available in its github page.

Now let's go back to R.

# install COMUNET
devtools::install_github("ScialdoneLab/COMUNET/COMUNET")
# You may fail to install COMUNET with the error of SMDTools being missing and cannot be installed.
# This is because SMDTools has been removed from the CRAN repository.
# In that case, one can install the package from its archive, as following
install.packages("R.utils")
install.packages("https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz")
devtools::install_github("ScialdoneLab/COMUNET/COMUNET")

library(COMUNET)
# read CellPhoneDB complex and gene info
complex_input <- read.csv("https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/complex_input.csv")
complex_input$complex_name <- gsub("_", " ", complex_input$complex_name)
gene_input <- read.csv("https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv")

# read CellPhoneDB output
CellPhoneDB_output <- read.csv("out/significant_means.txt", sep = "\t", check.names = F)
CellPhoneDB_output <- CellPhoneDB_output[!duplicated(CellPhoneDB_output$interacting_pair),]
rownames(CellPhoneDB_output) <- CellPhoneDB_output$interacting_pair
CellPhoneDB_output$receptor_a <- CellPhoneDB_output$receptor_a == "True"
CellPhoneDB_output$receptor_b <- CellPhoneDB_output$receptor_b == "True"

# Convert to COMUNET format
interactions <- convert_CellPhoneDB_output(CellPhoneDB_output = CellPhoneDB_output,
                                           complex_input = complex_input,
                                           gene_input = gene_input)

# Run communication clusters analysis
lrp_clusters <- lrp_clustering(weight_array = interactions$weight_array,
                               ligand_receptor_pair_df = interactions$ligand_receptor_pair_df,
                               nodes = interactions$nodes)

# Do heatmap
plot_cluster_heatmap(dissim_matrix = lrp_clusters$dissim_matrix,
                    lrp_clusters = lrp_clusters$clusters)

# Plot the communication mode of cluster 10 pairs as an example
cluster <- "cluster 10"
plot_communication_graph(LRP = cluster,
                         weight_array = lrp_clusters$weight_array_by_cluster[,,cluster],
                         ligand_receptor_pair_df = interactions$ligand_receptor_pair_df,
                         nodes = interactions$nodes,
                         is_cluster = T)
 
