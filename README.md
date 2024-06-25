# RescueSC
Single cell RNA sequencing (scRNA-seq) is a revolutionary technology which allows to detect gene expression signatures specific to cell types. However, it yields highly complex data and there is no gold standard in how to analyze it. Thus, we developed RescueSC tool, which tackles 3 specific problematic areas in scRNA-seq data secondary analysis:
1. Demultiplexing
2. Quality control (QC)
3. Automating cluster annotation with cell types

The pipeline consists of a number of functions and is generally divided into three parts: RescueTag for demultiplexing; RescueCluster for keeping more biological information during clustering; and AutoClusterType for annotation of clusters with cell types. Here you will find a tutorial on how to analyze your scRNA-seq data with RescueSC.

## Dependencies
In this repository you can find an R file called "Dependencies_for_RescueSC.R". Inside this file there are all the necessary libraries which should be installed before using RescueSC. It is also written how to install the packages, which are not in CRAN. Versions of all the packages and the version of R which was used for the RescueSC development are written in the file as well. 

## Tutorial
First of all, you need to download the "Final_RescueSC.R" file and load all the functions from it:
```
source("Final_RescueSC.R")
```
RescueSC works with a Seurat object, which contains read counts (nCount_RNA) and number of genes detected (nFeature_RNA) for every cell.

NB! RescueSC works with the raw data after aligment. It is important to make sure that no preliminary quality control filtering was performed on the data after the alignment! 

Therefore, the you need to create a variable with the Seurat object.
```
scAD<-readRDS(file.choose())
```
### RescueTag
Now, if your experiment did not have multiplexing, skip this part and go directly to RescueCluster.

The first step that RescueTag does is calculating the most abundant sample tag for every cell (First Best) and the difference between the most abundant and the second most abundant tag (Delta). Based on these parameters the putative tags are assigned to every cell. 
For this purpose a ScoringTags function is used. As an input it takes:
1. A path to a table with sample tag read counts for every cell. The table should be in csv format.
2. Numbers of the first and the last column in the table, which represent the sample tags used for the experiment (x, y)
3. Numbers of first and the last tags used in the experiment (e.g. if in the experiment were used SampleTag1-SampleTag4, you should write 1,4; and if, for instance, SampleTag6-SampleTag10 were used from the same kit, then the respective numbers will be 6,10) (first_tag_number, last_tag_number)
4. Seurat object (scAD)

```
scAD <- ScoringTags(path_to_table, x, y, first_tag_number, last_tag_number, scAD)
```
