# R version 4.3.1

# devtools version 2.4.5
# remotes version 2.5.0
# BiocManager) version 1.30.23
# dplyr version 1.1.4
# Seurat version 5.1.0
# patchwork version 1.2.0
# ggplot2 version 3.5.1
# SingleR version 2.2.0
# celldex version 1.10.1
# RColorBrewer version 1.1-3
# SingleCellExperiment version 1.22.0
# pheatmap version 1.0.12
# scMCA version 0.2.0
# DoubletDecon version 1.1.6
# tidyverse version 2.0.0
# harmony version 1.2.0
# scales version 1.3.0
# ggridges version 0.5.6
# mixtools version 2.0.0
# lessR version 4.3.0

install.packages("devtools")
install.packages("remotes") 
install.packages("BiocManager")
install.packages("dplyr")
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
install.packages("patchwork")
install.packages("ggplot2")
BiocManager::install("SingleR")
BiocManager::install("celldex")
install.packages("RColorBrewer")
BiocManager::install("SingleCellExperiment")
install.packages("pheatmap")
devtools::install_github("ggjlab/scMCA")
devtools::install_github('EDePasquale/DoubletDecon')
install.packages("tidyverse")
BiocManager::install("harmony", version = "3.8")
install.packages("scales")
install.packages("ggridges")
install.packages("mixtools")
install.packages("lessR")
