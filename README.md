# RescueSC
Single cell RNA sequencing (scRNA-seq) is a revolutionary technology which allows to detect gene expression signatures specific to cell types. However, it yields highly complex data and there is no gold standard in how to analyze it. Thus, we developed RescueSC tool, which tackles 3 specific problematic areas in scRNA-seq data secondary analysis:
1. Demultiplexing
2. Quality control (QC)
3. Automating cluster annotation with cell types
The pipeline consists of a number of functions and is generally divided into three parts: RescueTag for demultiplexing; RescueCluster for keeping more biological information during clustering; and AutoClusterType for annotation of clusters with cell types.
