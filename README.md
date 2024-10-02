# RescueSC
Single cell RNA sequencing (scRNA-seq) is a revolutionary technology which allows to detect gene expression signatures specific to cell types. However, it yields highly complex data and there is no gold standard in how to analyze it. Thus, we developed RescueSC tool, which tackles 3 specific problematic areas in scRNA-seq data secondary analysis:
1. Demultiplexing
2. Quality control (QC)
3. Automating cluster annotation with cell types

The pipeline consists of a number of functions and is generally divided into three parts: RescueTag for demultiplexing; RescueCluster for keeping more biological information during clustering; and AutoClusterType for annotation of clusters with cell types. Here you will find a tutorial on how to analyze your scRNA-seq data with RescueSC.

## Dependencies
In this repository you can find an R file called "Dependencies_for_RescueSC.R". Inside this file all the necessary libraries which should be installed before using RescueSC and their versions at the moment of the tool's development are listed. It is also written how to install the packages. The version of R which was used for the RescueSC development is written in the file as well. 

After installation, load all the libraries into your workspace.
```
library(devtools)
library(remotes)
library(BiocManager)
library(dplyr)
library(Seurat) 
library(patchwork)
library(ggplot2)
library(SingleR) 
library(celldex) 
library(RColorBrewer)
library(SingleCellExperiment) 
library(pheatmap)
library(shinythemes)
library(scMCA) 
library(tidyverse) 
library(harmony) 
library(scales)
library(ggridges)
library(mixtools)
library(lessR)
```
## Tutorial
First of all, you need to download the "Final_RescueSC.R" file and load all the functions from it:
```
source("Final_RescueSC.R")
```
RescueSC works with a Seurat object, which contains read counts (nCount_RNA) and number of genes detected (nFeature_RNA) for every cell.

NB! RescueSC works with the raw data after aligment. It is important to make sure that no preliminary quality control filtering was performed on the data after the alignment! 

Therefore, you need to create a variable with the Seurat object.
```
seurat_obj<-readRDS(file.choose()) #choose the proper RDS file from your file system
```
### RescueTag
Now, if your experiment did not have multiplexing, skip this part and go directly to RescueCluster.

The first step that RescueTag does is calculating the most abundant sample tag for every cell (Best) and the difference between the most abundant and the second most abundant tag (Delta). Then the ratio Delta/Best is calculated and based on this parameter the putative tags are assigned to every cell. 
For this purpose a ScoringTags function is used. As an input it takes:
1. A path to a table with sample tag read counts for every cell. The table should be in csv format.
2. Numbers of the first and the last column in the table, which represent the sample tags used for the experiment (for example, if a table has 9 columns, and only columns from 4th to 7th have reads of the sample tags, used in the experiment, then the numbers will be 4 and 7). (x, y)
3. Numbers of first and the last tags used in the experiment (e.g. if in the experiment were used SampleTag1-SampleTag4, you should write 1,4) (first_tag_number, last_tag_number)
4. Seurat object (seurat_obj)

NB! Make sure that the table does not have any additional text or information apart from the rows and columns with the sample tag read counts! 
```
seurat_obj <- ScoringTags(path_to_table, x, y, first_tag_number, last_tag_number, seurat_obj)
```
After assigning putative tags to the cells the first quality control step can be performed. We filter out cells which contain less than 200 genes detected, and those that had less than 3 reads. 
```
seurat_obj <- PreQCFilter(seurat_obj)
```
Then we normalize the First Best parameter to the whole sequencing depth (i.e. gene reads + sample tag reads).
```
seurat_obj <- NormToDepth(seurat_obj)
```
Before performing another quality control filtering step, it is necessary to look at the distribution of the First Best parameter normalized to the whole sequencing depth by running the following lines of code.
```
df_for_plot<-data.frame('First_best'=seurat_obj$NTD_first_best)
df_for_plot %>% ggplot(aes(x=First_best))+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+expand_limits(x=c(0,1))
```
The filtering of low-quality cells (i.e. low sequencing depth cells) is performed by the FilterLowSeqDepth function and depends on the obtained distribution. 

FilterLowSeqDepth takes only two parameters as input:
1. Seurat object (seurat_obj)
2. Distribution. The default value is 0.5.

So if the distribution is bimodal and the local minimum between peaks falls around 0.5, you should go with the default value of the "distribution" parameter.
If the distribution is normal, the function uses 95% confidence interval (CI) for filtering, and in case of bimodal distribution it implement expectation-maximization (EM) algorithm for identifying 2 separate peaks in the distribution (as a reference for implementing the EM algorithm [this webpage](https://rpubs.com/H_Zhu/246450)  was taken). 

After filtering First Best and Delta are normalized to the sample tag sequencing depth for bringing them to 0-1 scale. Also, the Ratio parameter is calculated.
```
seurat_obj <- FilterLowSeqDepth(seurat_obj, distribution) # for the distribution write either "normal" or "bimodal"
seurat_obj <- NormTagQCParams(seurat_obj)
```
Now a threshold for tagging can be set and the tags can be assigned to cells. For this, the user needs to run the FinalTagging function, which takes the Seurat object and the Ratio threshold as an input. The default value is 0.5, meaning that we require that the difference between Best and Second Best is at least half of the Best. The user can set up custom thresholds as well. NB! The custom threshold cannot be equal to or lower than 0.1!

```
seurat_obj<-FinalTagging(seurat_obj, threshold_ratio) # for the thresholds insert the variables which you obtained in the previous step
```

### RescueCluster
The RescueCluster strategy is based on a regular [Seurat clustering](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) with two major differences:
1. No filtering of cells with mitochondrial gene percent per cell (MGPC)
2. Utilization of fuzzy clustering by using [Harmony](https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06c_integration_harmony.md)

If you did not use RescueTag, you have to perform first the primary filtering step with a PreQCFilter function to filter out cells which contain less than 200 genes detected, and those that had less than 3 reads.
```
seurat_obj <- PreQCFilter(seurat_obj)
```
Now, we need to calculate MGPC and visualize its distribution.
```
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^mt") #calculate percent mt
Idents(seurat_obj) <- "orig.ident"
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=1)
```
Then we can perform all the steps from the abovementioned Harmony tutorial for normalization, dimensionality reduction and clustering.

NB! If you don't have multiplexing in your experiment, use the RunHarmony function only if you integrate different samples into one Seurat object. In such case, put the Sample ID variable into group.by.vars argument. However, if the Seurat object contains only one sample, then SCTransform normalization will reduce technical batch effects, and RunHarmony should not be run. 
```
merged_seurat <- seurat_obj %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% SCTransform(vars.to.regress = c("percent.mt"))
merged_seurat <- RunPCA(merged_seurat, assay = "SCT")
ElbowPlot(merged_seurat)
unfiltered_seurat <- RunHarmony(merged_seurat, group.by.vars = c("Final_tags"), reduction = "pca", assay.use = "SCT", reduction.save = "harmony", dims.use=1:n) # instead of n insert a number of PCs to use for clustering based on ElbowPlot
unfiltered_seurat <- RunUMAP(unfiltered_seurat, reduction = "harmony", assay = "SCT", dims = 1:n) 
unfiltered_seurat <- FindNeighbors(object = unfiltered_seurat, reduction = "harmony", dims = 1:n) 
unfiltered_seurat <- FindClusters(unfiltered_seurat, resolution = 0.5) # you can change the resolution parameter; the bigger the resolution - the more clusters you will obtain
```
After these procedures you will obtain a Seurat object with clustering information inside it. Now, you can visualize the results of clustering, coloring cells with high MGPC in grey. 

NB! This function can be applied to the data with maximum 52 clusters!
```
RefinedDimPlot(unfiltered_seurat, percent.mt_threshold) # instead of percent.mt_threshold insert a number which will represent the maximum MGPC after which the cells will be considered as having high MGPC level (e.g. 20, 30, or 35, etc.)
```
You can as well get some information about quality of every cluster by displaying the following plots.
```
VlnPlot(unfiltered_seurat, features = 'percent.mt', group.by = 'seurat_clusters') # displays MGPC distribution in every cluster
VlnPlot(unfiltered_seurat, features = 'nFeature_RNA', group.by = 'seurat_clusters') # displays number of genes per cell distribution in every cluster
ClusterQC(unfiltered_seurat, "seurat_clusters") # displays the same information in a scatter plot manner
```
### AutoClusterType
After clustering cell type annotation must be performed. We use [scMCA](https://github.com/ggjlab/scMCA) mouse cell atlas for this purpose. 

At first the following lines of code should be ran to obtain the annotations.

NB! scMCA function takes a lot of memory, so make sure that you have at least 32 GB of RAM in your PC or laptop!
```
querry <- GetAssayData(unfiltered_seurat, slot = "data")
querryMatrix <- as.matrix(querry)
mca_result <- scMCA(scdata = querryMatrix, numbers_plot = 3)
unfiltered_seurat[["MCA"]] <- mca_result$scMCA
```
Then we can utilize AutoClusterType function to robustly annotate every cluster with a cell type. This function takes the following as an input:
1. A path to a reference file (you can find it in this repository by the name of "Reference_annotations_scMCA.csv")
2. mca_result (a variable generated by previous lines of code)
3. unfiltered_seurat

```
unfiltered_seurat <- AutoClusterType(reference_file_path, mca_result, unfiltered_seurat)
```
Now we can visualize the clusters again, but from now on instead of the cluster numbers the cluster cell types will be displayed.
```
Idents(unfiltered_seurat) <- "cell.ID"
RefinedDimPlot(unfiltered_seurat, percent.mt_threshold) # instead of percent.mt_threshold insert a number which will represent the maximum MGPC after which the cells will be considered as having high MGPC level (e.g. 20, 30, or 35, etc.)
ClusterQC(unfiltered_seurat, "cell.ID")
```
