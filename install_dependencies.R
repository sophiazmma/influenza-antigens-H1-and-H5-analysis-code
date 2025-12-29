# 00. install.packages --------------------------------------------------------
#0.1 安装系统依赖库
#sudo apt-get update  # 更新软件源
#sudo apt-get install -y build-essential pkg-config libcurl4-openssl-dev libssl-dev  # 关键依赖[1][2][9]

#0.2 installing
install.packages("curl")      # 基础网络库[1][6][9]
install.packages("openssl")   # 加密支持[1][6][9]
install.packages("httr")      # HTTP请求工具[1][6][8]
install.packages("plotly")    # 可视化依赖[6][8][9]
install.packages("Seurat")  
install.packages("clustree") 
install.packages("eulerr")
install.packages("ggforce")
#解决clustree安装问题：
# sudo apt-get update
# sudo apt-get install aptitude
# sudo aptitude install libmount-dev libblkid-dev libselinux1-dev uuid-dev libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev
# no
## cran下载
cran_packages <- c(
  "stringr", "survival", "glmnet", "survminer", "timeROC",
  "data.table", "ggpubr", "dplyr", "patchwork", "readr",
  "tibble", "ggplot2", "tidyverse", "future", "pheatmap",
  "RColorBrewer", "hdf5r", "devtools", "Seurat"
)

# 批量安装（自动跳过已安装的包）
install.packages(cran_packages, dependencies = TRUE)
install.packages("harmony")
  
# Bioconductor包安装
  if (!require("BiocManager")) install.packages("BiocManager")

bioc_packages <- c(
  "msigdbr", "clusterProfiler", "GSVA", "AUCell", "monocle",
  "cellchat", "SingleR", "celldex"
)
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)
BiocManager::install("cellchat", dependencies = TRUE)
BiocManager::install("AUCell", dependencies = TRUE, force = TRUE)
BiocManager::install("Biostrings", dependencies = TRUE, force = TRUE)
BiocManager::install("alakazam", dependencies = TRUE, force = TRUE)
BiocManager::install("shazam", dependencies = TRUE, force = TRUE)
BiocManager::install("tigger", dependencies = TRUE, force = TRUE)
BiocManager::install("stringr", dependencies = TRUE, force = TRUE)
BiocManager::install("tidyverse", dependencies = TRUE, force = TRUE)
BiocManager::install("ggpubr", dependencies = TRUE, force = TRUE)

##Github特殊包安装
# 若 DoubletFinder/copykat 等未通过上述方式安装，手动从 GitHub 安装
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
remotes::install_github("navinlabcode/copykat")
remotes::install_github("aertslab/SCENIC") 

# 从 CRAN 安装依赖
install.packages(c("ggplot2", "dplyr", "tidyr", "viridis", "ggrepel", "igraph", "Rcpp", "RcppEigen"))
# 从 Bioconductor 安装依赖
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "limma", "BiocGenerics"))

install.packages("devtools")
pkgbuild::check_build_tools(debug = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages("scRepertoire")
install.packages("vegan")
install.packages("optparse")
install.packages("circlize")
install.packages("ragg")
install.packages("tidyverse")
install.packages("treemap")
install.packages("viridis")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2) 安装 dowser（以及常用依赖）
BiocManager::install(c("dowser", "ggtree", "treeio"))

# 3) 载入测试
library(dowser)
packageVersion("dowser")
