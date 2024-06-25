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

Therefore, you need to create a variable with the Seurat object.
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
After assigning putative tags to the cells the first quality control step can be performed. We filter out cells which contain less than 200 genes detected, and those that had less than 3 reads. 
```
scAD <- PreQCFilter(scAD)
```
Then we normalize the First Best parameter to the whole sequencing depth (i.e. gene reads + sample tag reads).
```
scAD <- NormToDepth(scAD)
```
Before performing another quality control filtering step, it is necessary to look at the distribution of the First Best parameter normalized to the whole sequencing depth by running the following lines of code.
```
df_for_plot<-data.frame('First_best'=scAD$NTD_first_best)
df_for_plot %>% ggplot(aes(x=First_best))+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+expand_limits(x=c(0,1))
```
The filtering of low-quality cells (i.e. low sequencing depth cells) is performed by the FilterLowSeqDepth function and depends on the obtained distribution. 

FilterLowSeqDepth takes only two parameters as input:
1. Seurat object (scAD)
2. Distribution. It can be either "normal" or "bimodal". The default value is 0.5.

So if the distribution is bimodal and the local minimum between peaks falls around 0.5, you should choose the default value of the "distribution" parameter.
If the distribution is normal, the function uses 95% confidence interval (CI) for filtering, and in case of bimodal distribution it implement expectation-maximization (EM) algorithm for identifying 2 separate peaks in the distribution (the implementation of the EM algorithm was taken from [https://rpubs.com/H_Zhu/246450](https://rpubs.com/H_Zhu/246450)
```
scAD <- FilterLowSeqDepth(scAD, distribution)
```
