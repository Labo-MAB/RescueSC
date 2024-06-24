# R version 4.3.1

library(devtools) #version 2.4.5
library(remotes) #version 2.5.0
library(BiocManager)#version 1.30.23
library(dplyr)#version 1.1.4
library(Seurat) # this one should be installed with remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE); version 5.1.0
library(patchwork)#version 1.2.0
library(ggplot2)#version 3.5.1
library(SingleR) # this one should be installed with BiocManager::install("SingleR"); version 2.2.0
library(celldex) # this one should be installed with BiocManager::install("celldex"); version 1.10.1
library(RColorBrewer)#version 1.1-3
library(SingleCellExperiment) # this one should be installed with BiocManager::install("SingleCellExperiment") version 1.22.0
library(pheatmap)#version 1.0.12
library(scMCA) # this one should be installed with devtools::install_github("ggjlab/scMCA"); version 0.2.0
library(DoubletDecon) # this one should be installed with devtools::install_github('EDePasquale/DoubletDecon'); version 1.1.6
library(tidyverse)#version 2.0.0
library(harmony) # this one should be installed with BiocManager::install("harmony", version = "3.8"); version 1.2.0
library(scales)#version 1.3.0
library(ggridges)#version 0.5.6
library(mixtools)#version 2.0.0
library(lessR)#version 4.3.0