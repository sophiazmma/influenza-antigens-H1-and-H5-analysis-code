
# 01. loading library ------------------------------------------------------------
# 定义需要加载的包名向量
pkgs <- c(
  "stringr", "survival", "glmnet", "survminer", "data.table",
  "ggpubr", "dplyr", "patchwork", "matrixStats", "readr", "tibble", "ggplot2",
  "future", "pheatmap", "msigdbr", 
  "Seurat",  "harmony", 
   "RColorBrewer"
)

# 检查包是否已安装，没安装的先安装，然后加载
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
library(Seurat) 
library(dplyr)
library(patchwork)
library(Matrix)
library(vegan) 
library(future)
library(Biostrings)
library(GenomicAlignments)
library(alakazam)
library(shazam)
library(tigger)
library(dplyr)
library(stringr)
library(readr)
library(clustree)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(scales)
library(purrr) 
library(optparse)
library(circlize)
library(RColorBrewer)
library(treemapify)
library(treemap)
library(readr)
library(ggpubr)
library(viridis)


# 02. project1-- different HA（H1 and H5） from influenza --------------------------------------------------
## 02.0 scbcr annotation --------------------------------------------------

### Igblast --------------------------------------------------
#igblast download: https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz
#tar -xvzf ncbi-igblast-1.22.0-x64-linux.tar.gz
#export PATH=$PATH:/usr/local/igblast/ncbi-igblast-1.22.0/bin

### convert imgt to IgBLAST --------------------------------------------------
# IGHV/IGKV/IGLV均出现问题
# sed -i '/^>/! s/[^ATCGNatcgn]//g' IGHV.fasta
# sed -i '/^>/! s/[^ATCGNatcgn]//g' IGKV.fasta
# sed -i '/^>/! s/[^ATCGNatcgn]//g' IGLV.fasta

# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in IGHV.fasta -out mouse_HV
# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in IGHJ.fasta -out mouse_HJ
# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in IGHD.fasta -out mouse_HD

# /root/miniconda3/envs/bcr_analysis/bin/igblastn   
# -query ./cell_reports/VDJ/G19-013_SampleID_2_13mar19/filtered_contig.fasta 
# -germline_db_V /Q/10_imgt/mouse/IGHV.fasta  
# -germline_db_D /Q/10_imgt/mouse/IGHD.fasta   
# -germline_db_J /Q/10_imgt/mouse/IGHJ.fasta   
# -auxiliary_data /root/miniconda3/envs/bcr_analysis/share/igblast/optional_file/mouse_gl.aux  
# -organism mouse 
# -ig_seqtype Ig 
# -outfmt 19
# -num_threads 8   
# -out Heavy_bcr_igblast.out

# /root/miniconda3/envs/bcr_analysis/bin/MakeDb.py igblast 
# -i Heavy_bcr_igblast.out 
# -s ./cell_reports/VDJ/G19-013_SampleID_2_13mar19/filtered_contig.fasta 
# -r /Q/10_imgt/mouse/IGHV.fasta  
#    /Q/10_imgt/mouse/IGHD.fasta   
#    /Q/10_imgt/mouse/IGHJ.fasta 
# --extended（执行错误，尚未解决问题）

# cat IGHV.fasta IGKV.fasta IGLV.fasta > mouse_gl_V
# cat IGHJ.fasta IGKJ.fasta IGLJ.fasta > mouse_gl_J
# cat IGHD.fasta >mouse_gl_D

# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_V -out mouse_gl_V
# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_D -out mouse_gl_D
# /Q/software/ncbi-igblast-1.22.0/bin/makeblastdb -parse_seqids -dbtype nucl -in mouse_gl_J -out mouse_gl_J

# 终端执行
# flu_H1_mouse1:  /root/miniconda3/envs/bcr_analysis/bin/igblastn -query ./cell_reports/VDJ/G19-013_SampleID_2_13mar19/filtered_contig.fasta  -germline_db_V /Q/10_imgt/mouse/mouse_gl_V -germline_db_D /Q/10_imgt/mouse/mouse_gl_D  -germline_db_J /Q/10_imgt/mouse/mouse_gl_J   -auxiliary_data /root/miniconda3/envs/bcr_analysis/share/igblast/optional_file/mouse_gl.aux  -organism mouse  -ig_seqtype Ig -outfmt 19 -num_threads 8  -out mouse_H1_1.out
# flu_H1_mouse2: /root/miniconda3/envs/bcr_analysis/bin/igblastn -query ./cell_reports/VDJ/G19-013_SampleID_1_9apr19/filtered_contig.fasta  -germline_db_V /Q/10_imgt/mouse/mouse_gl_V -germline_db_D /Q/10_imgt/mouse/mouse_gl_D  -germline_db_J /Q/10_imgt/mouse/mouse_gl_J   -auxiliary_data /root/miniconda3/envs/bcr_analysis/share/igblast/optional_file/mouse_gl.aux  -organism mouse  -ig_seqtype Ig -outfmt 19 -num_threads 8  -out mouse_H1_2.out
# naive_mouse1: /root/miniconda3/envs/bcr_analysis/bin/igblastn -query ./cell_reports/VDJ/G19-013_SampleID_1_16apr19/filtered_contig.fasta  -germline_db_V /Q/10_imgt/mouse/mouse_gl_V -germline_db_D /Q/10_imgt/mouse/mouse_gl_D  -germline_db_J /Q/10_imgt/mouse/mouse_gl_J   -auxiliary_data /root/miniconda3/envs/bcr_analysis/share/igblast/optional_file/mouse_gl.aux  -organism mouse  -ig_seqtype Ig -outfmt 19 -num_threads 8  -out naive_mouse_1.out
# naive_mouse2: /root/miniconda3/envs/bcr_analysis/bin/igblastn -query ./cell_reports/VDJ/G19-035_SampleID_1_27may19/filtered_contig.fasta  -germline_db_V /Q/10_imgt/mouse/mouse_gl_V -germline_db_D /Q/10_imgt/mouse/mouse_gl_D  -germline_db_J /Q/10_imgt/mouse/mouse_gl_J   -auxiliary_data /root/miniconda3/envs/bcr_analysis/share/igblast/optional_file/mouse_gl.aux  -organism mouse  -ig_seqtype Ig -outfmt 19 -num_threads 8  -out naive_mouse_2.out

### Annotation --------------------------------------------------
# merge_cellranger_igblast
# # merge_cellranger_igblast <- function(cellranger_csv, igblast_tsv) {
#   # 读取 CellRanger 输出
#   cellranger <- read_csv(cellranger_csv, show_col_types = FALSE)
#   
#   # 读取 IgBLAST AIRR 格式输出
#   igblast <- read_tsv(igblast_tsv, show_col_types = FALSE)
#   
#   # 合并，contig_id 对应 sequence_id
#   merged <- cellranger %>%
#     left_join(igblast, by = c("contig_id" = "sequence_id"))
#   
#   # 避免字段冲突，重命名
#   merged <- merged %>%
#     rename(
#       productive_cellranger = productive.x,
#       productive_igblast = productive.y,
#       v_gene_cellranger = v_gene,
#       v_call_igblast = v_call,
#       j_gene_cellranger = j_gene,
#       j_call_igblast = j_call
#     )
#   
#   return(merged)
# }
# merge_cellranger_igblast函数不稳定，每次运行结果不同
#### flu_H1_mouse1 --------------------------------------------------
#flu_H1_mouse1_merged <- merge_cellranger_igblast(
#  cellranger_csv = "./cell_reports/VDJ/G19-013_SampleID_2_13mar19/filtered_contig_annotations.csv",
#  igblast_tsv   = "./cell_reports/VDJ/annotation_results/mouse_H1_1.out"
#)
#flu_H1_mouse1_merged.clean <- process_and_format_bcr(flu_H1_mouse1_merged)
#write_csv(flu_H1_mouse1_merged.clean, "./cell_reports/VDJ/annotation_results/flu_H1_mouse1_merged_bcr_annotations.csv")

#### flu_H1_mouse2 --------------------------------------------------
# flu_H1_mouse2_merged <- merge_cellranger_igblast(
#   cellranger_csv = "./cell_reports/VDJ/G19-013_SampleID_1_9apr19/filtered_contig_annotations.csv",
#   igblast_tsv   = "./cell_reports/VDJ/annotation_results/mouse_H1_2.out"
# )
# flu_H1_mouse2_merged.clean <- process_and_format_bcr(flu_H1_mouse2_merged)
# write_csv(flu_H1_mouse2_merged.clean, "./cell_reports/VDJ/annotation_results/flu_H1_mouse2_merged_bcr_annotations.csv")

#### naive_mouse1 --------------------------------------------------
# naive_mouse1_merged <- merge_cellranger_igblast(
#   cellranger_csv = "./cell_reports/VDJ/G19-013_SampleID_1_16apr19/filtered_contig_annotations.csv",
#   igblast_tsv   = "./cell_reports/VDJ/annotation_results/naive_mouse_1.out"
# )
# naive_mouse1_merged.clean <- process_and_format_bcr(naive_mouse1_merged)
# write_csv(naive_mouse1_merged.clean, "./cell_reports/VDJ/annotation_results/naive_mouse1_merged_bcr_annotations.csv")

#### naive_mouse2 --------------------------------------------------
# naive_mouse2_merged <- merge_cellranger_igblast(
#   cellranger_csv = "./cell_reports/VDJ/G19-035_SampleID_1_27may19/filtered_contig_annotations.csv",
#   igblast_tsv   = "./cell_reports/VDJ/annotation_results/naive_mouse_2.out"
# )
# naive_mouse2_merged.clean <- process_and_format_bcr(naive_mouse2_merged)
# write_csv(naive_mouse2_merged.clean, "./cell_reports/VDJ/annotation_results/naive_mouse2_merged_bcr_annotations.csv")

#2. Calculating nearest neighbor distances based on heavy chains

library(shazam)#专门用于体细胞高频突变和克隆分析的包
library(ggplot2)
data(ExampleDb, package="alakazam")
native_mouse2_heavy_parse <- read.table("./cell_reports/VDJ/annotation_results/heavy_native_mouse2_parse-select.tsv",sep = "\t",header = T,stringsAsFactors = F)
native_mouse2_dist_ham <- distToNearest(native_mouse2_heavy_parse, sequenceColumn="junction", #计算每条序列到其“最近邻居”的距离 (使用汉明距离模型)
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=1)
# 这是最核心的函数调用。
# distToNearest: shazam包的函数，它的作用是：对于数据集中的每一条序列，
#                找到另一条与它共享相同V基因和J基因、且序列差异最小的序列（即“最近邻居”），
#                并计算它们之间的距离。
# sequenceColumn="junction": 指定使用 junction 区域 (包含了最具多样性的CDR3区) 的序列来进行距离计算。
# vCallColumn="v_call", jCallColumn="j_call": 这是一个重要的约束条件，
#                它要求只有在V基因和J基因完全相同的情况下，两条序列才会被认为是“邻居”并计算距离。
#                这是定义克隆的生物学基础。
# model="ham": 指定距离计算模型为“汉明距离 (Hamming distance)”，即两个等长序列中，
#              对应位置上不同字符的数量。这是最简单、最常用的模型。
# normalize="len": 将计算出的汉明距离除以序列的长度，进行归一化。这使得不同长度的CDR3可以被公平比较。
# nproc=1: 使用1个CPU核心进行计算。
native_mouse2_output_ham <- findThreshold(native_mouse2_dist_ham$dist_nearest, method="density")
# findThreshold: shazam包的另一个核心函数。
# dist_ham$dist_nearest: 输入是上一步计算出的所有“最近邻居距离”的集合。
# 这个距离的分布通常是双峰的：
#   - 第一个峰 (靠近0)：代表了“克隆内”的距离，因为同一克隆的成员彼此非常相似。
#   - 第二个峰 (距离较远)：代表了“克隆间”的距离，因为不同克隆的序列差异很大。
# method="density": 告诉函数使用基于核密度估计的方法来自动识别这两个峰之间的“谷底”，
#                   这个“谷底”就是区分两种距离的最佳分界线，即我们需要的阈值。
# dist_s5f <- distToNearest(native_mouse2_heavy_parse, sequenceColumn="junction",
#                           vCallColumn="v_call", jCallColumn="j_call",
#                           model="hh_s5f", normalize="none", nproc=1)
# 下面两行代码只是为了演示另一种名为 "hh_s5f" 的模型，它考虑了不同核苷酸的突变倾向，
# 在您的实际分析中，通常先从简单的 "ham" 模型开始。
# output_s5f <- findThreshold(dist_s5f$dist_nearest, method="density")
native_mouse2_output_ham@threshold
# output_s5f@threshold# show the threshold
# 打印出两种模型计算出的阈值。例如，output_ham@threshold 的值可能是 0.16。
# 这个 0.16 就是您下一步在 `DefineClones.py` 脚本中需要通过 `--dist` 参数提供的值。


file_paths <- c(
  "fluH1_mouse1" = "./cell_reports/VDJ/annotation_results/heavy_fluH1_mouse1_parse-select.tsv",
  "fluH1_mouse2" = "./cell_reports/VDJ/annotation_results/heavy_fluH1_mouse2_parse-select.tsv",
  "fluH5_mouse"  = "./cell_reports/VDJ/annotation_results/heavy_fluH5_mouse_parse-select.tsv",
  "naive_mouse1" = "./cell_reports/VDJ/annotation_results/heavy_native_mouse1_parse-select.tsv",
  "naive_mouse2" = "./cell_reports/VDJ/annotation_results/heavy_native_mouse2_parse-select.tsv"
)
print(file.exists(file_paths))
calculate_clonal_threshold_from_path <- function(file_path, sample_name) {
  
  print(paste("Calculating threshold for sample:", sample_name))
  
  # 直接使用传入的完整路径
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NA)
  }
  
  heavy_data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  dist_ham <- distToNearest(heavy_data, 
                            sequenceColumn = "junction",
                            vCallColumn = "v_call", 
                            jCallColumn = "j_call",
                            model = "ham", 
                            normalize = "len", 
                            nproc = 10)
  
  output_ham <- findThreshold(dist_ham$dist_nearest, method = "density")
  
  return(output_ham@threshold)
}
all_thresholds <- mapply(
  FUN = calculate_clonal_threshold_from_path,
  file_path = file_paths,       # 第一个参数，对应函数中的 file_path
  sample_name = names(file_paths) # 第二个参数，对应函数中的 sample_name
)
print(all_thresholds)

###combine Heavy and Light chain------

sample_prefixes <- c(
  "fluH1_mouse1", 
  "fluH1_mouse2", 
  "fluH5_mouse", 
  "naive_mouse1", 
  "naive_mouse2"
)

input_directory <- "./cell_reports/VDJ/annotation_results"

output_directory <- "./cell_reports/VDJ/annotation_results"

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

for (prefix in sample_prefixes) {
  
  print(paste("--- Processing sample:", prefix, "---"))
  
  heavy_chain_file <- file.path(input_directory, paste0("heavy_", prefix, "_10X_clone-pass_germ-pass.tsv"))
  light_chain_file <- file.path(input_directory, paste0("light_", prefix, "_parse-select.tsv"))
  
  if (!file.exists(heavy_chain_file) || !file.exists(light_chain_file)) {
    warning(paste("Input file(s) missing for sample:", prefix, ". Skipping."))
    next 
  }
  
  heavy_df <- read_tsv(heavy_chain_file, col_types = cols(.default = "c"))
  light_df <- read_tsv(light_chain_file, col_types = cols(.default = "c"))
  
  light_info_to_merge <- light_df %>%
    select(
      cell_id,
      light_sequence_id = sequence_id,
      light_fwr1 = fwr1,
      light_fwr2 = fwr2,
      light_fwr3 = fwr3,
      light_fwr4 = fwr4,
      light_cdr1 = cdr1,
      light_cdr2 = cdr2,
      light_cdr3 = cdr3,
      light_consensus_count = consensus_count,
      light_umi_count = umi_count,
      light_v_call_10x = v_call_10x,
      light_j_call_10x = j_call_10x,
      light_junction_10x = junction_10x,
      light_junction_10x_aa = junction_10x_aa
    ) %>%
    group_by(cell_id) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # 2.5. 合并轻重链数据
  merged_df <- left_join(heavy_df, light_info_to_merge, by = "cell_id")
  
  # 2.6. (关键步骤) 构建输出文件名并保存
  #    文件名格式为: [prefix].merge.bcr.tsv
  output_filename <- paste0(prefix, ".merge.bcr.tsv")
  output_filepath <- file.path(output_directory, output_filename)
  
  # 使用 write_tsv() 保存文件
  write_tsv(merged_df, output_filepath)
  
  print(paste("Successfully merged and saved to:", output_filepath))
}

print("--- All samples processed and saved. ---")

list.files(path = output_directory, pattern = "\\.merge\\.bcr\\.tsv$")


###translate DNA TO aa------
sample_prefixes <- c(
  "fluH1_mouse1", 
  "fluH1_mouse2", 
  "fluH5_mouse", 
  "naive_mouse1", 
  "naive_mouse2"
)

# 数据目录
data_directory <- "./cell_reports/VDJ/annotation_results"

# 要翻译的列
heavy_chain_fields <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")
light_chain_fields <- c("light_fwr1", "light_cdr1", "light_fwr2", "light_cdr2", "light_fwr3", "light_cdr3")

all_bcr_data_translated <- lapply(sample_prefixes, function(prefix) {
  
  message(paste0("--- Processing and translating sample: ", prefix, " ---"))
  
  # 输入文件路径
  file_path <- file.path(data_directory, paste0(prefix, ".merge.bcr.tsv"))
  if (!file.exists(file_path)) {
    warning(paste("File not found for sample:", prefix, ". Skipping."))
    return(NULL)
  }
  
  # 读取
  bcr_df <- read_tsv(file_path, col_types = cols(.default = "c"))
  
  # 对每个需要翻译的字段逐列翻译
  for (field in c(heavy_chain_fields, light_chain_fields)) {
    if (field %in% colnames(bcr_df)) {
      aa_field <- paste0(field, "_aa")
      bcr_df[[aa_field]] <- sapply(bcr_df[[field]], function(seq_nt) {
        if (is.na(seq_nt) || nchar(seq_nt) < 3) {
          return(NA)
        } else {
          # 捕获翻译异常
          tryCatch({
            alakazam::translateDNA(seq_nt, trim = TRUE)
          }, error = function(e) {
            NA
          })
        }
      })
    }
  }
  
  message(paste0("Translation complete for sample: ", prefix))
  
  # --- 保存结果 ---
  out_file <- file.path(data_directory, paste0(prefix, ".merge.final.bcr.tsv"))
  write_tsv(bcr_df, out_file)
  message(paste0("Saved translated file: ", out_file))
  
  return(bcr_df)
})


###shm calculation for Heavy Chain------
# --- 步骤 0: 加载所有必要的库 ---
# 确保这些包已经安装: install.packages(c("alakazam", "shazam", "dplyr", "readr"))
library(alakazam)
library(shazam)
library(dplyr)
library(readr)

# --- 步骤 1: (关键) 设置您的工作目录 ---
# 请将这里的路径替换为您存放 .merge.final.bcr.tsv 文件的实际目录
work_dir <- "./cell_reports/VDJ/annotation_results"

# --- 步骤 2: 查找目录中所有需要处理的文件 ---
cat("Finding files to process...\n")

# 查找所有以 ".merge.final.bcr.tsv" 结尾的文件
bcr_files <- list.files(path = work_dir, 
                        pattern = "\\.merge\\.final\\.bcr\\.tsv$", 
                        full.names = TRUE)

# 检查是否找到了文件
if (length(bcr_files) == 0) {
  stop("No '*.merge.final.bcr.tsv' files found in the specified directory. Please check your 'work_dir' path.")
} else {
  cat("Found the following files to process:\n")
  print(basename(bcr_files))
}

# --- 步骤 3: 使用 for 循环，逐一处理每个文件 ---
cat("\nStarting mutation quantification for each file...\n")

for (file_path in bcr_files) {
  
  cat('--------------------------------------------------\n')
  cat('Processing file:', basename(file_path), '\n')
  
  # 1. 读取数据
  db <- tryCatch({
    read_tsv(file_path, col_types = cols(.default = "c"))
  }, error = function(e) {
    warning(paste("Could not read file:", basename(file_path), ". Skipping. Error:", e$message))
    return(NULL)
  })
  
  # 如果文件读取失败，则跳到下一个文件
  if (is.null(db)) {
    next
  }
  
  # 2. 检查必需的列是否存在
  required_cols <- c("sequence_alignment", "germline_alignment_d_mask", "locus")
  if (!all(required_cols %in% colnames(db))) {
    warning(paste("File", basename(file_path), "is missing required columns. Skipping."))
    next
  }
  
  # 3. 筛选出重链进行计算
  heavy_chains_db <- db %>% filter(locus == "IGH")
  
  if (nrow(heavy_chains_db) > 0) {
    # **计算总突变数量 (COUNT)**
    count_db <- observedMutations(heavy_chains_db, 
                                  sequenceColumn = "sequence_alignment", 
                                  germlineColumn = "germline_alignment_d_mask",
                                  regionDefinition=NULL,
                                  frequency = FALSE)
    
    # **计算总突变频率 (FREQUENCY)**
    freq_db <- observedMutations(heavy_chains_db, 
                                 sequenceColumn = "sequence_alignment", 
                                 germlineColumn = "germline_alignment_d_mask",
                                 regionDefinition=NULL,
                                 frequency = TRUE, # 计算频率
                                 combine = TRUE)
    # 4. 准备用于合并的突变信息
    mutation_counts <- count_db %>% 
      select(sequence_id, MU_COUNT_HEAVY_seq_r = mu_count_seq_r,MU_COUNT_HEAVY_seq_s = mu_count_seq_s)
    
    mutation_freqs <- freq_db %>%
      select(sequence_id, MU_FREQ_HEAVY_TOTAL = mu_freq)
    
    # 5. 将突变数量和频率信息合并回原始的完整数据框
    db_with_mutations <- db %>%
      left_join(mutation_counts, by = "sequence_id") %>%
      left_join(mutation_freqs, by = "sequence_id")
    
  } else {
    warning(paste("No heavy chains (locus == 'IGH') found in file:", basename(file_path)))
    # 即使没有重链，也创建空的列以保持所有输出文件格式一致
    db_with_mutations <- db %>%
      mutate(
        MU_COUNT_HEAVY_TOTAL = NA_integer_, 
        MU_FREQ_HEAVY_TOTAL = NA_real_
      )
  }
  
  # 6. 构建新的输出文件名并保存
  #    将 ".tsv" 替换为 ".shm.tsv"
  output_filepath <- gsub("\\.tsv$", ".shm.tsv", file_path)
  
  # 使用 write_tsv() 保存文件
  write_tsv(db_with_mutations, output_filepath)
  
  cat('Successfully calculated mutations and saved to:', basename(output_filepath), '\n')
}

cat('--------------------------------------------------\n')
cat('All samples processed and saved successfully!\n')


###shm calculation for Light chain------

# --- 步骤 0: 安装和加载所有必要的R包 ---

# 加载所有包
library(dplyr)
library(readr)
library(shazam)

# --- 步骤 1: 加载、合并并准备所有数据 ---

# 1.1 定义文件路径和名称
base_path <- "./cell_reports/VDJ/annotation_results" # 表示当前工作目录
file_paths <- c(
  fluH1_mouse1 = file.path(base_path, "light_fluH1_mouse1_parse-select.tsv"),
  fluH1_mouse2 = file.path(base_path, "light_fluH1_mouse2_parse-select.tsv"),
  fluH5_mouse = file.path(base_path, "light_fluH5_mouse_parse-select.tsv"),
  naive_mouse1 = file.path(base_path, "light_naive_mouse1_parse-select.tsv"),
  naive_mouse2 = file.path(base_path, "light_naive_mouse2_parse-select.tsv")
)

# 1.2 逐一读取文件
list_of_dfs <- lapply(file_paths, function(path) {
  if (file.exists(path)) {
    read_tsv(path, col_types = cols(.default = "c"))
  } else {
    warning("找不到文件: ", path, "。已跳过。")
    NULL
  }
})

list_of_dfs <- list_of_dfs[!sapply(list_of_dfs, is.null)]

if (length(list_of_dfs) == 0) {
  stop("错误：所有指定的文件都未能读取。请检查文件名是否正确。")
}

# 1.3 创建2个合并文件
list_of_dfs$fluH1_mouse_combined <- bind_rows(list_of_dfs$fluH1_mouse1, list_of_dfs$fluH1_mouse2)
list_of_dfs$naive_mouse_combined <- bind_rows(list_of_dfs$naive_mouse1, list_of_dfs$naive_mouse2)

message("所有7个数据样本已成功加载并准备就绪。")
print(names(list_of_dfs))

# --- 步骤 2: 设置输出目录 ---
results_directory <- "light_shm_results"
if (!dir.exists(results_directory)) {
  dir.create(results_directory)
}

# --- 步骤 3: 定义一个可以为单个样本计算SHM的函数 (封装了您的核心逻辑) ---
calculate_shm_for_light_chain <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- 开始为样本处理轻链SHM: ", sample_name, " ---")
  
  # 检查必需的列
  required_cols <- c("sequence_id", "sequence_alignment", "germline_alignment")
  if (!all(required_cols %in% colnames(sample_data))) {
    warning("样本 '", sample_name, "' 缺少必需的列，跳过处理。")
    return(NULL)
  }
  
  db <- sample_data
  
  if (nrow(db) > 0) {
    message("找到 ", nrow(db), " 条轻链序列进行计算。")
    
    # --- 分步计算体细胞高频突变 (SHM) ---
    
    # 计算总突变数量 (COUNT)
    message("  正在计算突变数量 (COUNT)...")
    count_db <- observedMutations(db, 
                                  sequenceColumn = "sequence_alignment", 
                                  germlineColumn = "germline_alignment",
                                  regionDefinition = NULL,
                                  frequency = FALSE,
                                  nproc = 1)
    
    # 计算总突变频率 (FREQUENCY)
    message("  正在计算突变频率 (FREQUENCY)...")
    freq_db <- observedMutations(db, 
                                 sequenceColumn = "sequence_alignment", 
                                 germlineColumn = "germline_alignment",
                                 regionDefinition = NULL,
                                 frequency = TRUE,
                                 nproc = 1)
    
    # 分别整理数量和频率的结果
    mutation_counts <- count_db %>%
      rowwise() %>%
      mutate(MU_COUNT_LIGHT_TOTAL = sum(c_across(starts_with("mu_count")), na.rm = TRUE)) %>%
      ungroup() %>%
      select(sequence_id, MU_COUNT_LIGHT_TOTAL)
    
    mutation_freqs <- freq_db %>%
      rowwise() %>%
      mutate(MU_FREQ_LIGHT_TOTAL = sum(c_across(starts_with("mu_freq")), na.rm = TRUE)) %>%
      ungroup() %>%
      select(sequence_id, MU_FREQ_LIGHT_TOTAL)
    
    # 将数量和频率信息依次合并回原始数据框
    db_with_mutations <- db %>%
      left_join(mutation_counts, by = "sequence_id") %>%
      left_join(mutation_freqs, by = "sequence_id")
    
  } else {
    warning("在样本 '", sample_name, "' 中没有找到任何序列。")
    db_with_mutations <- db %>%
      mutate(
        MU_COUNT_LIGHT_TOTAL = NA_integer_, 
        MU_FREQ_LIGHT_TOTAL = NA_real_
      )
  }
  
  # 构建新的输出文件名并保存
  output_filename <- paste0(sample_name, "_light_with_shm.tsv")
  output_filepath <- file.path(output_dir, output_filename)
  
  write_tsv(db_with_mutations, output_filepath)
  
  message("成功计算轻链SHM并保存到: ", output_filename)
}

# --- 步骤 4: 循环处理列表中的所有数据框 ---
for (name in names(list_of_dfs)) {
  current_data <- list_of_dfs[[name]]
  tryCatch({
    # 调用我们封装好的函数
    calculate_shm_for_light_chain(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}

message("\n\n所有样本处理完毕！请检查 '", results_directory, "' 文件夹。")


## 02.1 samples --------------------------------------------------
### flu_H5_mouse ----------------------------------------------------------

# reading data 
flu_H5_mouse.counts <- Read10X(data.dir = "./AgB/outs/per_sample_outs/AgB/count/sample_filtered_feature_bc_matrix/")
flu_H5_mouse.obj <- CreateSeuratObject(counts = flu_H5_mouse.counts, project = "flu_H5_mouse",min.cells = 3, min.features = 200)

#QC
flu_H5_mouse.obj[["percent.mt"]] <- PercentageFeatureSet(flu_H5_mouse.obj, pattern = "^mt-")
flu_H5_mouse.obj <- subset(flu_H5_mouse.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) 

#reading BCR and filtering 
flu_H5_mouse.bcr <- read.delim("./cell_reports/VDJ/annotation_results/fluH5_mouse.merge.final.bcr.shm.tsv", sep = "\t", row.names = NULL, fileEncoding = "UTF-8")
 process_and_format_bcr <- function(raw_bcr_df) {
  
  # --- 步骤 1: 初始全局质控 ---
  # 使用 as.character() 将列强制转为文本，然后再转大写
  # 这可以同时处理因子(factor)和各种大小写的字符串("true", "True")
  clean_bcr <- raw_bcr_df %>%
    filter(toupper(as.character(productive)) == "T", 
           toupper(as.character(high_confidence)) == "TRUE",
           !is.na(umis))
  
  # 如果在最开始的筛选后就没有数据了，提前终止并警告
  if(nrow(clean_bcr) == 0) {
    warning("在初始筛选 (productive/high_confidence) 后没有数据剩下。这极不寻常，请手动检查 flu_H5.bcr 中这两列的内容！")
    return(tibble()) # 返回一个空的tibble
  }
  
  # --- 步骤 2: 清洗重链 (IGH) ---
  heavy_chains_clean <- clean_bcr %>%
    filter(chain == "IGH") %>%
    group_by(barcode) %>%
    slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # --- 步骤 3: 清洗轻链 (IGK & IGL) ---
  light_chains_clean <- clean_bcr %>%
    filter(chain %in% c("IGK", "IGL")) %>%
    group_by(barcode) %>%
    slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  final_bcr_table <- full_join(heavy_chains_clean, light_chains_clean, by = "barcode")
  
  return(final_bcr_table)
}
# flu_H5.bcr.clean <- process_and_format_bcr(flu_H5.bcr)
# print(flu_H5.bcr.clean)

#add bcr to rna
# head(colnames(flu_H5_mouse.obj))
# [1] "AAACCAAAGAACACCC-1" "AAACCAAAGATTACGT-1" "AAACCAAAGCGAAGCA-1" "AAACCAAAGCTCGAGT-1"
# [5] "AAACCAAAGCTGACGT-1" "AAACCAAAGGCATCGT-1"
# head(flu_H5.bcr.clean$barcode)
# [1] "AAACCAAAGAACACCC-1" "AAACCAAAGCTCGAGT-1" "AAACCAAAGCTGACGT-1" "AAACCAAAGGGCATAG-1"
# [5] "AAACCAGCAACGTGGG-1" "AAACCAGCACCCTTTA-1"

add_bcr_to_seurat <- function(seurat_obj, cleaned_bcr_df) {
  
  # 1. 获取 Seurat 对象当前的元数据
  seurat_metadata <- seurat_obj@meta.data
  
  # 2. 将 Seurat 元数据的行名（即 barcode）变成一列
  seurat_metadata$barcode_seurat <- rownames(seurat_metadata)
  
  # 3. 处理 BCR 数据：删除 sequence_id 列中后缀 "_contig_X" 并重命名为 barcode
  cleaned_bcr_df$barcode <- sub("_contig_\\d+", "", cleaned_bcr_df$sequence_id)
  
  # 4. 将 barcode 列移到第一列
  cleaned_bcr_df <- cleaned_bcr_df[, c("barcode", setdiff(names(cleaned_bcr_df), "barcode"))]
  
  # 5. 使用 left_join 进行合并
  merged_metadata <- left_join(seurat_metadata, cleaned_bcr_df, by = c("barcode_seurat" = "barcode"))
  
  # 6. 检查是否有 NA（未匹配的条形码），并发出警告
  if (sum(is.na(merged_metadata$barcode_seurat)) > 0) {
    warning(paste(sum(is.na(merged_metadata$barcode_seurat)), "cells in Seurat object did not match any BCR data."))
  }
  
  # 7. 保证列名唯一（如果合并后有重复列名）
  colnames(merged_metadata) <- make.unique(colnames(merged_metadata))
  
  # 8. 重新设置行名
  rownames(merged_metadata) <- merged_metadata$barcode_seurat
  
  # 9. 将合并后的新元数据重新赋给 Seurat 对象
  seurat_obj@meta.data <- merged_metadata
  
  # 10. 输出消息
  message("BCR数据已成功合并到Seurat对象的元数据中！")
  
  # 11. 返回更新后的 Seurat 对象
  return(seurat_obj)
}

flu_H5_mouse.obj_bcr <- add_bcr_to_seurat(flu_H5_mouse.obj, flu_H5_mouse.bcr)

#Normalization and Data scaling
flu_H5_mouse.obj_bcr <- NormalizeData(flu_H5_mouse.obj_bcr)
flu_H5_mouse.obj_bcr <- FindVariableFeatures(flu_H5_mouse.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
#dir.create("./processed_seurat_objects")
saveRDS(flu_H5_mouse.obj_bcr, file = "./processed_seurat_objects/flu_H5_mouse.obj_bcr.processed.rds")


### flu_H1_mouse1 --------------------------------------------------------

# reading data 
flu_H1_mouse1.counts <- Read10X(data.dir = "./cell_reports/GEX/G19-013_SampleID_3_27feb19/")
flu_H1_mouse1.obj <- CreateSeuratObject(counts = flu_H1_mouse1.counts, project = "flu_H1_mouse1",min.cells = 3, min.features = 200)

# QC
flu_H1_mouse1.obj[["percent.mt"]] <- PercentageFeatureSet(flu_H1_mouse1.obj, pattern = "^mt-")
flu_H1_mouse1.obj <- subset(flu_H1_mouse1.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) 

#reading BCR and filtering
flu_H1_mouse1.bcr <- read.delim("./cell_reports/VDJ/annotation_results/fluH1_mouse1.merge.final.bcr.shm.tsv",, sep = "\t", row.names = NULL, fileEncoding = "UTF-8")

#flu_H1_mouse1.bcr.clean <- process_and_format_bcr(flu_H1_mouse1.bcr)

#add bcr to rna
flu_H1_mouse1.obj_bcr <- add_bcr_to_seurat(flu_H1_mouse1.obj, flu_H1_mouse1.bcr)

#Normalization and Data scaling
flu_H1_mouse1.obj_bcr <- NormalizeData(flu_H1_mouse1.obj_bcr)
flu_H1_mouse1.obj_bcr <- FindVariableFeatures(flu_H1_mouse1.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
saveRDS(flu_H1_mouse1.obj_bcr, file = "./processed_seurat_objects/flu_H1_mouse1.obj_bcr.processed.rds")

### flu_H1_mouse2 --------------------------------------------------------

# reading data 
flu_H1_mouse2.counts <- Read10X(data.dir = "./cell_reports/GEX/G19-013_SampleID_4_1apr19/")
flu_H1_mouse2.obj <- CreateSeuratObject(counts = flu_H1_mouse2.counts, project = "flu_H1_mouse2",min.cells = 3, min.features = 200)

# QC
flu_H1_mouse2.obj[["percent.mt"]] <- PercentageFeatureSet(flu_H1_mouse2.obj, pattern = "^mt-")
flu_H1_mouse2.obj <- subset(flu_H1_mouse2.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) 

# reading BCR and filtering 
flu_H1_mouse2.bcr <- read.delim("./cell_reports/VDJ/annotation_results/fluH1_mouse2.merge.final.bcr.shm.tsv", sep = "\t", row.names = NULL, fileEncoding = "UTF-8")

#flu_H1_mouse2.bcr.clean <- process_and_format_bcr(flu_H1_mouse2.bcr)

#add bcr to rna
flu_H1_mouse2.obj_bcr <- add_bcr_to_seurat(flu_H1_mouse2.obj, flu_H1_mouse2.bcr)

#Normalization and Data scaling
flu_H1_mouse2.obj_bcr <- NormalizeData(flu_H1_mouse2.obj_bcr)
flu_H1_mouse2.obj_bcr <- FindVariableFeatures(flu_H1_mouse2.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
saveRDS(flu_H1_mouse2.obj_bcr, file = "./processed_seurat_objects/flu_H1_mouse2.obj_bcr.processed.rds")

### naive_mouse1 ------------------------------------------------------------
# reading data 
naive_mouse1.counts <- Read10X(data.dir = "./cell_reports/GEX/G19-013_SampleID_1_27feb19/")
naive_mouse1.obj <- CreateSeuratObject(counts = naive_mouse1.counts, project = "naive_mouse1",min.cells = 3, min.features = 200)

# QC
naive_mouse1.obj[["percent.mt"]] <- PercentageFeatureSet(naive_mouse1.obj, pattern = "^mt-")
naive_mouse1.obj <- subset(naive_mouse1.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) 

# reading BCR and filtering 
naive_mouse1.bcr <- read.delim("./cell_reports/VDJ/annotation_results/naive_mouse1.merge.final.bcr.shm.tsv", sep = "\t", row.names = NULL, fileEncoding = "UTF-8")

#naive_mouse1.bcr.clean <- process_and_format_bcr(naive_mouse1.bcr)

#add bcr to rna
naive_mouse1.obj_bcr <- add_bcr_to_seurat(naive_mouse1.obj, naive_mouse1.bcr)

#Normalization and Data scaling
naive_mouse1.obj_bcr <- NormalizeData(naive_mouse1.obj_bcr)
naive_mouse1.obj_bcr <- FindVariableFeatures(naive_mouse1.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
saveRDS(naive_mouse1.obj_bcr, file = "./processed_seurat_objects/naive_mouse1.obj_bcr.processed.rds")


### naive_mouse2 ------------------------------------------------------------

# reading data 
naive_mouse2.counts <- Read10X(data.dir = "./cell_reports/GEX/G19-013_SampleID_2_29apr19/")
naive_mouse2.obj <- CreateSeuratObject(counts = naive_mouse2.counts, project = "naive_mouse2",min.cells = 3, min.features = 200)

# QC
naive_mouse2.obj[["percent.mt"]] <- PercentageFeatureSet(naive_mouse2.obj, pattern = "^mt-")
naive_mouse2.obj <- subset(naive_mouse2.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) 

# reading BCR and filtering 
naive_mouse2.bcr <- read.delim("./cell_reports/VDJ/annotation_results/naive_mouse2.merge.final.bcr.shm.tsv", sep = "\t", row.names = NULL, fileEncoding = "UTF-8")

#naive_mouse2.bcr.clean <- process_and_format_bcr(naive_mouse2.bcr)

# add bcr to rna
naive_mouse2.obj_bcr <- add_bcr_to_seurat(naive_mouse2.obj, naive_mouse2.bcr)

# Normalization and Data scaling
naive_mouse2.obj_bcr <- NormalizeData(naive_mouse2.obj_bcr)
naive_mouse2.obj_bcr <- FindVariableFeatures(naive_mouse2.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
saveRDS(naive_mouse2.obj_bcr, file = "./processed_seurat_objects/naive_mouse2.obj_bcr.processed.rds")

# naive_mouse1.obj <- NormalizeData(naive_mouse1.obj)
# naive_mouse1.obj <- FindVariableFeatures(naive_mouse1.obj, selection.method = "vst", nfeatures = 2000)
# saveRDS(naive_mouse1.obj, file = "./processed_seurat_objects/test/naive_mouse1.obj.rds")
# 
# naive_mouse2.obj <- NormalizeData(naive_mouse2.obj)
# naive_mouse2.obj <- FindVariableFeatures(naive_mouse2.obj, selection.method = "vst", nfeatures = 2000)
# saveRDS(naive_mouse2.obj, file = "./processed_seurat_objects/test/naive_mouse2.obj.rds")
# 
# flu_H1_mouse1.obj <- NormalizeData(flu_H1_mouse1.obj)
# flu_H1_mouse1.obj <- FindVariableFeatures(flu_H1_mouse1.obj, selection.method = "vst", nfeatures = 2000)
# saveRDS(flu_H1_mouse1.obj, file = "./processed_seurat_objects/test/flu_H1_mouse1.obj.rds")
# 
# flu_H1_mouse2.obj_bcr <- NormalizeData(flu_H1_mouse2.obj_bcr)
# flu_H1_mouse2.obj_bcr <- FindVariableFeatures(flu_H1_mouse2.obj_bcr, selection.method = "vst", nfeatures = 2000)
# saveRDS(flu_H1_mouse2.obj, file = "./processed_seurat_objects/test/flu_H1_mouse2.obj.rds")
# 
# flu_H5_mouse.obj <- NormalizeData(flu_H5_mouse.obj)
# flu_H5_mouse.obj <- FindVariableFeatures(flu_H5_mouse.obj, selection.method = "vst", nfeatures = 2000)
# saveRDS(flu_H5_mouse.obj, file = "./processed_seurat_objects/test/flu_H5_mouse.obj.rds")

## 02.2 Quality control, QC ------------------------------------------------

### load and merge seurat ------------------------------------------------------------
# flu_mouse_object.list <- merge(flu_H1_mouse2.obj_bcr, y = c(flu_H1_mouse1.obj.bcr, flu_H5_mouse.obj_bcr,
#                                                        naive_mouse1.obj_bcr, naive_mouse2.obj_bcr))
processed_dir <- "./processed_seurat_objects/"
flu_mouse_files <- list.files(path = processed_dir, pattern = "\\.rds$", full.names = TRUE)
flu_mouse_object.list <- lapply(flu_mouse_files, readRDS)
names(flu_mouse_object.list) <- gsub("\\.rds$", "", basename(flu_mouse_files))

### Integrate  ----------------------------------------------------

flu_mouse_object.list <- mapply(
  function(obj, nm) RenameCells(obj, add.cell.id = nm),
  flu_mouse_object.list, names(flu_mouse_object.list),
  SIMPLIFY = FALSE
)

#SCTransform
for (i in 1:length(flu_mouse_object.list)) {
  flu_mouse_object.list[[i]][["percent.mt"]] <- PercentageFeatureSet(flu_mouse_object.list[[i]], pattern = "^mt-")
  flu_mouse_object.list[[i]] <- SCTransform(flu_mouse_object.list[[i]], 
                                            verbose = TRUE,
                                            vars.to.regress = c("percent.mt"))
}

#FindIntegrationAnchors
features <- SelectIntegrationFeatures(object.list = flu_mouse_object.list, nfeatures = 3000)
features <- setdiff(features, grep("(?i)^(ighm|ighg|igha|ighe|ighd)$", features, value = TRUE))

#PrepSCTIntegration
flu_mouse_object.list <- PrepSCTIntegration(object.list = flu_mouse_object.list, 
                                            anchor.features = features)

# ★ 每个对象先做 PCA（用于 rpca 邻域）
flu_mouse_object.list <- lapply(flu_mouse_object.list, function(x) {
  DefaultAssay(x) <- "SCT"
  RunPCA(x, features = features, npcs = 50, verbose = FALSE)
})

#reference index
ref.samples <- c("flu_H1_mouse1.obj_bcr.processed","flu_H1_mouse2.obj_bcr.processed","naive_mouse1.obj_bcr.processed","naive_mouse2.obj_bcr.processed")
ref.idx <- which(names(flu_mouse_object.list) %in% ref.samples)

#
flu_mouse.anchors <- FindIntegrationAnchors(
  object.list = flu_mouse_object.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30,
  k.anchor = 10,
  k.score = 30
)

#
flu_mouse.integrated <- IntegrateData(
  anchorset = flu_mouse.anchors,
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 60
)

### interation before and after view ----------------------------------------------

merged.obj <- merge(flu_mouse_object.list[[1]], y = flu_mouse_object.list[2:5])
DefaultAssay(merged.obj) <- "RNA"
merged.obj <- NormalizeData(merged.obj)
merged.obj <- FindVariableFeatures(merged.obj,selection.method = "vst",nfeatures = 3000)
merged.obj <- ScaleData(merged.obj, vars.to.regress = "percent.mt")
merged.obj <- RunPCA(merged.obj, npcs = 50, verbose = FALSE)
merged.obj <- FindNeighbors(merged.obj, dims = 1:30)
merged.obj <- RunUMAP(merged.obj, dims = 1:30)

merged.obj <- FindClusters(merged.obj, resolution = 0.6)
plotS0a <- DimPlot(merged.obj, group.by = "orig.ident") + 
  ggtitle("Before integration")+
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15)   
  )
plotS0d <- DimPlot(merged.obj, group.by = "orig.ident", split.by = "orig.ident") + 
  ggtitle("Before integration")+
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15)   
  )

DefaultAssay(flu_mouse.integrated) <- "integrated"
flu_mouse.integrated <- RunPCA(flu_mouse.integrated, npcs = 50, verbose = FALSE)
flu_mouse.integrated <- FindNeighbors(flu_mouse.integrated, dims = 1:30)
flu_mouse.integrated <- RunUMAP(flu_mouse.integrated, dims = 1:30)
flu_mouse.integrated <- FindClusters(flu_mouse.integrated, resolution = 0.6)
DimPlot(flu_mouse.integrated, group.by = "orig.ident") + ggtitle("After integration")
plotS0b <- DimPlot(flu_mouse.integrated, group.by = "orig.ident") + 
  ggtitle("After integration")+
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15)   
  )

plotS0c <- DimPlot(flu_mouse.integrated, group.by = "orig.ident", split.by = "orig.ident")+ 
  ggtitle("After integration")+
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15)   
  )

plotS0 <- (plotS0a + plotS0b) /plotS0d/plotS0c
print(plotS0)
ggsave("./results/FigureS0.pdf", plot = plotS0, width = 12, height = 12)
ggsave("./results/FigureS0.png", plot = plotS0, width = 12, height = 12)

### Dimensionality reduction & clustering  ----------------------------------------------------

####Reduction ----------------------------------------------------
DefaultAssay(flu_mouse.integrated) <- "integrated"
flu_mouse.combined <- RunPCA(flu_mouse.integrated, verbose = FALSE)
plotS1a<- ElbowPlot(flu_mouse.combined, ndims=50, reduction="pca") 
plotS1b<- DimPlot(flu_mouse.combined, reduction = "pca", group.by="orig.ident")
plotS1 <- plotS1a+plotS1b
print(plotS1)
ggsave("./results/FigureS1.pdf", plot = plotS1, width = 12, height = 6)
ggsave("./results/FigureS1.png", plot = plotS1, width = 12, height = 6)

#clustree
flu_mouse.combined <- FindNeighbors(flu_mouse.combined, dims = 1:40)
flu_mouse.combined <- FindClusters(flu_mouse.combined, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
plotS2 <- clustree(flu_mouse.combined)
print(plotS2)
ggsave("./results/FigureS2.pdf", plot = plotS2, width = 18, height = 9)
ggsave("./results/FigureS2.png", plot = plotS2, width = 18, height = 9)

#Dimplot
plotS3a <- DimPlot(flu_mouse.combined, group.by = "integrated_snn_res.0.4", label = TRUE) + ggtitle("Resolution 0.4")
plotS3b <- DimPlot(flu_mouse.combined, group.by = "integrated_snn_res.0.6", label = TRUE) + ggtitle("Resolution 0.6")
plotS3c <- DimPlot(flu_mouse.combined, group.by = "integrated_snn_res.0.8", label = TRUE) + ggtitle("Resolution 0.8")
plotS3 <- plotS3a + plotS3b + plotS3c
print(plotS3)
ggsave("./results/FigureS3.pdf", plot = plotS3, width = 18, height = 6)
ggsave("./results/FigureS3.png", plot = plotS3, width = 18, height = 6)

####Cluster ----------------------------------------------------
flu_mouse.combined <- FindNeighbors(flu_mouse.combined, dims = 1:40)
flu_mouse.combined <- FindClusters(flu_mouse.combined, resolution = 0.8) 
flu_mouse.combined <- RunUMAP(flu_mouse.combined, dims = 1:40)
plotS4a <- DimPlot(flu_mouse.combined, reduction = "umap", label = TRUE)
plotS4b <- DimPlot(flu_mouse.combined, reduction = "umap", group.by='orig.ident',cols = c("#5471ab","#5471ab","#72ab77","#d1895d","#d1895d"))
plotS4 <- plotS4a+plotS4b
print(plotS4)
ggsave("./results/FigureS4.pdf", plot = plotS4, width = 18, height = 9)
ggsave("./results/FigureS4.png", plot = plotS4, width = 18, height = 9)

table(flu_mouse.combined@meta.data$seurat_clusters)
metadata <- flu_mouse.combined@meta.data
ncol(flu_mouse.combined)#统计项目中细胞数目
#[1] 28707
flu_mouse_cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(flu_mouse_cell_cluster,'./results/TableS1.csv',row.names = F)

## 02.3 Results ------------------------------------------------

### Figure 1: ------------------------------------------------
# Antigenic challenge with H5 and H1 influenza hemagglutinins dramatically reshapes the B cell transcriptional landscape
#### FigureS5------------------------------------------------
####Figure
features_B<- c("Cd19","Cd79a","Ebf1","Iglc3","Fcer2a")
features_T<- c("Cd3e", "Cd3g")
features_Macrophage<- c("C1qa", "Cd68","Trf")
features_Monocyte<- c("Fn1","Ace")
features_Neutrophil<- c("S100a9","S100a8","Retnlg")
features_NKC<- c("Ncr1","Klrk1")
#features_Dendriticcell<- c("Ccl22","Il4i1")
feature_Erythroid <-c("Hba-a1", "Hba-a2","Hbq1b")
features_Granulocyte <- c("Ccl4","Csf1")
genes_alls = c(features_B, features_T, features_Macrophage, features_Neutrophil, features_NKC,
               feature_Erythroid,features_Granulocyte)
plotS5d1 <- DotPlot(flu_mouse.combined, features = genes_alls, assay = "RNA")+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
print(plotS5d1)
ggsave("./results/FigureS5d1.pdf", plot = plotS5d, width = 18, height = 9)
ggsave("./results/FigureS5d1.png", plot = plotS5d, width = 18, height = 9)
#查看某些特定cluster高表达基因确定cluster所属于的细胞亚群
#cluster13
#cluster13_markers <- FindMarkers(flu_mouse.combined, ident.1 = 13)
#head(cluster13_markers)
#cluster13_up_markers <- subset(cluster13_markers, avg_log2FC > 0 & p_val_adj < 0.05)

# B_cell:1,2,3,4,5,6,10,11,12,13,14,16,17,18,21,26,28,29,30,32,33,34,36;
# T_cell:24,25;
# Macrophage:19,20,23,27;
# Neutrophil:15;
# NKC:7;
# Erythroid:35,37;
# Granulocyte:0,8,9,22,31

#####FigureS5a umap -----
flu_mouse.combined@meta.data$spleen_cell_subpopulations <- plyr::mapvalues(from = c (1,2,3,4,5,6,10,11,12,13,14,16,17,18,21,26,28,29,30,32,33,34,36,
                                                                                     24,25,
                                                                                     19,20,23,27,
                                                                                     15,
                                                                                     7,
                                                                                     35,37,
                                                                                     0,8,9,22,31),
                                                                         to = c(rep("B Cell",23), 
                                                                                rep("T Cell",2),
                                                                                rep("Macrophage",4),
                                                                                rep("Neutrophil",1),
                                                                                rep("NKC",1),
                                                                                rep("Erythroid",2),
                                                                                rep("Granulocyte",5)),
                                                                         x = flu_mouse.combined@meta.data$seurat_clusters)

#提取scRNAi_QC@meta.data的数据
flu_mouse.annotated_clusters <- flu_mouse.combined@meta.data
write.csv(flu_mouse.annotated_clusters,'./results/TableS2.csv',row.names = T)


umap_data <- data.frame(
  UMAP_1 = flu_mouse.combined@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = flu_mouse.combined@reductions$umap@cell.embeddings[, 2],
  cell_type = flu_mouse.combined$spleen_cell_subpopulations # 假设这是你的注释列
)

cell_type_medians <- umap_data %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

allcolour <- c("#4169E1","#c53b33","#228B22","#FF8C00", "#875b51", "#8d69b8","#d57dbe")

axis_origin_x <- min(umap_data$UMAP_1) - 1
axis_origin_y <- min(umap_data$UMAP_2) - 1
arrow_length <- 4 

plotS5a <- ggplot() +
  geom_point(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type), size = 1) +
  geom_label_repel(
    data = cell_type_medians, 
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    fontface = "bold",
    size = 7,
    box.padding = 0.5,
    point.padding = 0.5,
    label.size = 0.5,
    segment.color = 'grey50',
    family = "Arial" 
  ) +
  
  geom_segment(
    aes(x = axis_origin_x, y = axis_origin_y, xend = axis_origin_x + arrow_length, yend = axis_origin_y),
    arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1
  ) +
  geom_segment(
    aes(x = axis_origin_x, y = axis_origin_y, xend = axis_origin_x, yend = axis_origin_y + arrow_length),
    arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1
  ) +
  
  annotate("text", x = axis_origin_x + arrow_length / 2, y = axis_origin_y - 0.5, label = "UMAP 1", size = 6, fontface = "bold") +
  annotate("text", x = axis_origin_x - 0.5, y = axis_origin_y + arrow_length / 2, label = "UMAP 2", size = 6, fontface = "bold", angle = 90) +
  scale_color_manual(values = allcolour) +
  theme_void() +
  theme(
    legend.position = "none", # 因为有标签，所以不需要图例
    plot.background = element_rect(fill = "white", color = NA) # 确保背景是白色无边框
  )
print(plotS5a)
ggsave("./results/FigureS5a.pdf", plot = plotS5a, width = 9, height = 9)
ggsave("./results/FigureS5a.png", plot = plotS5a, width = 9, height = 9, dpi = 300)

#####FigureS5b cell counts------
table(flu_mouse.combined$orig.ident)
table(flu_mouse.combined@meta.data$spleen_cell_subpopulations)
meta_df <- flu_mouse.combined@meta.data
meta_df$spleen_cell_subpopulations <- unlist(meta_df$spleen_cell_subpopulations)
plot_data <- meta_df %>%
  dplyr::count(orig.ident, spleen_cell_subpopulations)
total_counts_data <- plot_data %>%
  group_by(orig.ident) %>%
  summarise(total_n = sum(n))
cell_type_colors <- c(
   "B Cell" = "#c53b33", 
   "Granulocyte" = "#4169e1", 
   "NKC" = "#228b22", 
   "T Cell" = "#8d69b8", 
   "Neutrophil" = "#ff8c00", 
   "Macrophage" = "#875b51", 
   "Erythroid" = "#d57dbe")
plotS5b <- ggplot(plot_data, aes(x = orig.ident, y = n, fill = spleen_cell_subpopulations)) +
  geom_col(position = "fill") +
  geom_text(
    data = total_counts_data,
    aes(x = orig.ident, y = 1, label = paste("n =", total_n)), 
    hjust = -0.1, 
    size = 5.5, 
    fontface = "bold",
    inherit.aes = FALSE 
  ) +
  coord_flip() +
  scale_fill_manual(values = cell_type_colors) +
  
  # --- 主要修改在这里 ---
  # 将 0.15 增加到 0.2 或 0.25，为右侧提供更多空间
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, .35))) + 
  
  labs(
    title = "Cell Proportions per Sample by Cell Type",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.title = element_text(size = 19, face = "bold"),
    legend.text = element_text(size = 18)
  )

print(plotS5b)
ggsave("./results/FigureS5b.pdf", plot = plotS5b, width = 9, height = 6)
ggsave("./results/FigureS5b.png", plot = plotS5b, width = 9, height = 6, dpi = 300)

#小插曲！统计每个样本中B细胞的细胞数量
b_cell_seurat_object <- subset(
  flu_mouse.combined, 
  idents = c(1,2,3,4,5,6,10,11,12,13,14,16,17,18,21,26,28,29,30,32,33,34,36)
)
b_cell_counts_from_subset <- table(b_cell_seurat_object$orig.ident)
print(b_cell_counts_from_subset)

#####FigureS5c heatmap ------
table(Idents(flu_mouse.combined))
flu_mouse.combined <- RenameIdents(object = flu_mouse.combined,
                                 "1" = "B Cell",
                                 "2" = "B Cell",
                                 "3" = "B Cell",
                                 "4" = "B Cell",
                                 "5" = "B Cell",
                                 "6" = "B Cell",
                                 "10" = "B Cell",
                                 "11" = "B Cell",
                                 "12" = "B Cell",
                                 "13" = "B Cell",
                                 "14" = "B Cell",
                                 "16" = "B Cell",
                                 "17" = "B Cell",
                                 "18" = "B Cell",
                                 "21" = "B Cell",
                                 "26" = "B Cell",
                                 "28" = "B Cell",
                                 "29" = "B Cell",
                                 "30" = "B Cell",
                                 "32" = "B Cell",
                                 "33" = "B Cell",
                                 "34" = "B Cell",
                                 "36" = "B Cell",
                                 "24" = "T Cell",
                                 "25" = "T Cell",
                                 "19" = "Macrophage",
                                 "20" = "Macrophage",
                                 "23" = "Macrophage",
                                 "27" = "Macrophage",
                                 "15" = "Neutrophil",
                                 "7" = "NKC",
                                 "35" = "Erythroid",
                                 "37" = "Erythroid",
                                 "0" = "Granulocyte",
                                 "8" = "Granulocyte",
                                 "9" = "Granulocyte",
                                 "22" = "Granulocyte",
                                 "31" = "Granulocyte")
table(Idents(flu_mouse.combined))
DefaultAssay(flu_mouse.combined) <- "RNA"
flu_mouse.combined <- JoinLayers(flu_mouse.combined, assay = "RNA")
dim(flu_mouse.combined[["RNA"]]$data) 
flu_mouse.combined_markers.gene <- FindAllMarkers(flu_mouse.combined, logfc.threshold = 0.25, test.use = "wilcox", assay = "RNA",
                                    min.pct = 0.1, only.pos =  TRUE)
top5_c <- flu_mouse.combined_markers.gene %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top8_c <- flu_mouse.combined_markers.gene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

flu_mouse.combined <- ScaleData(flu_mouse.combined, features = top8_c$gene)
DoHeatmap(flu_mouse.combined,top8_c$gene,size=8)

# 1) 确保 levels 与颜色名一一对应
Idents(flu_mouse.combined) <- factor(
  Idents(flu_mouse.combined),
  levels = c("B Cell","T Cell","Macrophage","Neutrophil","NKC","Erythroid","Granulocyte")
)
grp_cols <- c(
  "B Cell"      = "#f8766d",
  "T Cell"      = "#c49a00",
  "Macrophage"  = "#53b400",
  "Neutrophil"  = "#00c094",
  "NKC"         = "#00b6eb",
  "Erythroid"   = "#a58aff",
  "Granulocyte" = "#fb61d7"
)
grp_cols <- grp_cols[levels(Idents(flu_mouse.combined))]  # 对齐顺序

genes_use <- intersect(top8_c$gene, rownames(flu_mouse.combined))
# 2) 画图：连续色带用 gradientn（表达量），群组图例用 color_manual 且 drop=FALSE
plotS5c <- DoHeatmap(
  object       = flu_mouse.combined,
  features     = genes_use,
  size         = 5,
  draw.lines   = FALSE,
  group.by     = NULL,          # 用 Idents
  group.colors = grp_cols
) +
  scale_fill_gradientn(colors = c("#421863", "#437F8C", "#F4E755")) +  # 表达量色条
  scale_color_manual(                                             # ← 关键：强制包含全部群
    values = grp_cols,
    breaks = names(grp_cols),
    drop   = FALSE,
    name   = "Identity"
  ) +
  theme(
    axis.text.x      = element_text(family = "Arial", size = 20),
    axis.title.x     = element_text(family = "Arial", size = 20),
    axis.text.y      = element_text(family = "Arial", size = 20),
    legend.title     = element_text(family = "Arial", size = 18),
    legend.text      = element_text(family = "Arial", size = 18),
    legend.key.height = unit(0.8, "cm")
  )
print(plotS5c)
ggsave(filename="./results/FigureS5c.pdf", plot = plotS5c, width = 9, height = 18)
ggsave(filename="./results/FigureS5c.png", plot = plotS5c, width = 9, height = 18)

#####FigureS5d Dotplot------

plotS5d <- DotPlot(flu_mouse.combined, features = genes_alls) +
  theme(
    axis.text.x = element_text(family = "Arial", angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(family = "Arial", size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(family = "Arial"),
    legend.text = element_text(family = "Arial", size = 20),
    legend.title = element_text(family = "Arial", size = 20, face = "bold")
  ) +
  guides(
    color = guide_colorbar(
      title = "Expression", 
      title.theme = element_text(family = "Arial", size = 20),
      label.theme = element_text(family = "Arial", size = 20)
    ),
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(family = "Arial", size = 20),
      label.theme = element_text(family = "Arial", size = 20)
    )
  )
print(plotS5d)
ggsave("./results/FigureS5d.pdf", plot = plotS5d, width = 18, height = 6)
ggsave("./results/FigureS5d.png", plot = plotS5d, width = 18, height = 6)

####FigureS6 umap for 5 samples------
umap_data <- data.frame(
  UMAP_1 = flu_mouse.combined@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = flu_mouse.combined@reductions$umap@cell.embeddings[, 2],
  cell_type = flu_mouse.combined$spleen_cell_subpopulations,
  orig.ident = flu_mouse.combined$orig.ident
)
allcolour <- c("#c53b33","#4169E1","#228B22","#FF8C00", "#875b51", "#8d69b8","#d57dbe")
sample_list <- unique(umap_data$orig.ident)
plot_list <- list()
for (sample_id in sample_list) {
  print(paste("Generating plot for sample:", sample_id))
  current_sample_data <- filter(umap_data, orig.ident == sample_id)
  background_data <- filter(umap_data, orig.ident != sample_id)
  p_sample <- ggplot() +
    geom_point(
      data = background_data, 
      aes(x = UMAP_1, y = UMAP_2), 
      color = "grey85",
      size = 0.5 # 背景点可以更小
    ) +
    geom_point(
      data = current_sample_data, 
      aes(x = UMAP_1, y = UMAP_2, color = cell_type), 
      size = 1 
    ) +
    scale_color_manual(values = allcolour) +
    labs(title = sample_id) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none" 
    )
  plot_list[[sample_id]] <- p_sample
}
if (length(plot_list) == 5) {
  plotS6 <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = 'collect') 
    plot_annotation( 
      title = 'UMAP Projections per Sample',
      theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
    )
  
} else {
  plotS6 <- wrap_plots(plot_list, ncol = 3)
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = 'UMAP Projections per Sample',
      theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
    )
}
print(plotS6)
ggsave("./results/FigureS6.pdf", plot = plotS6, width = 18, height = 12)
ggsave("./results/FigureS6.png", plot = plotS6, width = 18, height = 12, dpi = 300)

####FigureS7------
B_cell_subset_flu_mouse.combined <- subset(flu_mouse.combined, idents = "B Cell")
#####B cell bach effect
B_cell_subset_flu_mouse.list <- SplitObject(B_cell_subset_flu_mouse.combined, split.by = "orig.ident")
for (i in seq_along(B_cell_subset_flu_mouse.list)) {
  B_cell_subset_flu_mouse.list[[i]][["percent.mt"]] <- PercentageFeatureSet(B_cell_subset_flu_mouse.list[[i]], pattern = "^mt-")
  B_cell_subset_flu_mouse.list[[i]] <- SCTransform(
    B_cell_subset_flu_mouse.list[[i]],
    vars.to.regress = "percent.mt", 
    verbose = FALSE
  )
}
features <- SelectIntegrationFeatures(object.list = B_cell_subset_flu_mouse.list, nfeatures = 3000)
B_cell_subset_flu_mouse.list <- PrepSCTIntegration(object.list = B_cell_subset_flu_mouse.list, anchor.features = features)
B_cell_subset_flu_mouse.anchors <- FindIntegrationAnchors(
  object.list = B_cell_subset_flu_mouse.list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:40,
  k.anchor = 10
)
B_cell_subset_flu_mouse.integrated <- IntegrateData(
  anchorset = B_cell_subset_flu_mouse.anchors,
  normalization.method = "SCT",
  dims = 1:40
)

B_cell_subset_flu_mouse.integrated <- RunPCA(B_cell_subset_flu_mouse.integrated, verbose = FALSE)
B_cell_subset_flu_mouse.integrated <- RunUMAP(B_cell_subset_flu_mouse.integrated, dims = 1:30)
B_cell_subset_flu_mouse.integrated <- FindNeighbors(B_cell_subset_flu_mouse.integrated, dims = 1:30)
B_cell_subset_flu_mouse.integrated <- FindClusters(B_cell_subset_flu_mouse.integrated, resolution = 0.6)

plotS7a<- DimPlot(B_cell_subset_flu_mouse.integrated, group.by = "orig.ident", reduction = "umap") +
  ggtitle("After B cell subset integration")
plotS7b<- DimPlot(B_cell_subset_flu_mouse.integrated, group.by = "seurat_clusters", reduction = "umap", label = TRUE) +
  ggtitle("After B cell subset integration")

####Before B cell interation

B_cell_subset_flu_mouse.combined <- FindNeighbors(B_cell_subset_flu_mouse.combined, dims = 1:30)
B_cell_subset_flu_mouse.combined <- FindClusters(
  B_cell_subset_flu_mouse.combined, 
  resolution = 0.6,
  graph.name = "integrated_snn"
)
B_cell_subset_flu_mouse.combined <- RunUMAP(B_cell_subset_flu_mouse.combined, dims = 1:30)
plotS7c <- DimPlot(B_cell_subset_flu_mouse.combined, group.by = "orig.ident", reduction = "umap") +
  ggtitle("Before B cell subset integration")
plotS7d <- DimPlot(B_cell_subset_flu_mouse.combined, group.by = "seurat_clusters", reduction = "umap", label = TRUE) +
  ggtitle("Before B cell subset integration")
plotS7 <- (plotS7a+plotS7b)/(plotS7c+plotS7d)
print(plotS7)
ggsave("./results/FigureS7.pdf", plot = plotS7, width = 18, height = 15)
ggsave("./results/FigureS7.png", plot = plotS7, width = 18, height = 15, dpi = 300)

####后面用Bcell批次校正后的数据进行分析
####FigureS8------
DefaultAssay(B_cell_subset_flu_mouse.integrated) <- "RNA"
B_cell_subset_flu_mouse.integrated <- NormalizeData(B_cell_subset_flu_mouse.integrated)
B_cell_subset_flu_mouse.integrated <- FindVariableFeatures(B_cell_subset_flu_mouse.integrated, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(B_cell_subset_flu_mouse.integrated)
B_cell_subset_flu_mouse.integrated <- ScaleData(B_cell_subset_flu_mouse.integrated, features = scale.genes)
B_cell_subset_flu_mouse.integrated <- RunPCA(B_cell_subset_flu_mouse.integrated, features = VariableFeatures(B_cell_subset_flu_mouse.integrated))
plotS8a <- ElbowPlot(B_cell_subset_flu_mouse.integrated, ndims=30, reduction="pca") 
plotS8b <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "pca")
plotS8 <- plotS8a+plotS8b
print(plotS8)
ggsave("./results/FigureS8.pdf", plot = plotS8, width = 12, height = 6)
ggsave("./results/FigureS8.png", plot = plotS8, width = 12, height = 6, dpi = 300)

####FigureS9------
B_cell_subset_flu_mouse.integrated <- FindNeighbors(B_cell_subset_flu_mouse.integrated, dims = 1:20)
B_cell_subset_flu_mouse.integrated <- FindClusters(B_cell_subset_flu_mouse.integrated, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
plotS9 <- clustree(B_cell_subset_flu_mouse.integrated)
print(plotS9)
ggsave("./results/FigureS9.pdf", plot = plotS9, width = 12, height = 9)
ggsave("./results/FigureS9.png", plot = plotS9, width = 12, height = 9, dpi = 300)

####FigureS10------
umap_data_bcell_subset <- data.frame(
  UMAP_1 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 2],
  bcell_cluster = B_cell_subset_flu_mouse.integrated$seurat_clusters, 
  orig.ident = B_cell_subset_flu_mouse.integrated$orig.ident
)
num_bcell_clusters <- length(unique(umap_data_bcell_subset$bcell_cluster))
bcell_cluster_colors <- viridis(num_bcell_clusters) 
sample_list <- unique(umap_data_bcell_subset$orig.ident)

plot_list_bcell <- list()
for (sample_id in sample_list) {
  print(paste("Generating B cell UMAP for sample:", sample_id))
  current_sample_data <- filter(umap_data_bcell_subset, orig.ident == sample_id)
  background_data <- filter(umap_data_bcell_subset, orig.ident != sample_id)
  
  p_sample <- ggplot() +
    geom_point(data = background_data, aes(x = UMAP_1, y = UMAP_2), color = "grey85", size = 0.5) +
    geom_point(data = current_sample_data, aes(x = UMAP_1, y = UMAP_2, color = bcell_cluster), size = 1) +
    
    # 在图例中添加标题，并确保点的尺寸更大更清晰
    scale_color_manual(
      values = bcell_cluster_colors, 
      name = "B Cell\nSubcluster", # 图例标题
      guide = guide_legend(override.aes = list(size = 4)) # 增大图例中点的大小
    ) +
    
    labs(title = sample_id) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right", # <-- 在这里显示图例
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  plot_list_bcell[[sample_id]] <- p_sample
}

plotS10 <- wrap_plots(plot_list_bcell, ncol = 3)

plotS10 <- plotS10 +
  plot_annotation(
    title = 'B Cell Subset UMAP Projections per Sample',
    theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
  )
print(plotS10)
ggsave("./results/FigureS10.pdf", plot = plotS10, width = 18, height = 12)
ggsave("./results/FigureS10.png", plot = plotS10, width = 18, height = 12, dpi = 300)

####FiguerS11------
B_cell_subset_flu_mouse.integrated <- FindNeighbors(B_cell_subset_flu_mouse.integrated, dims= 1:20)
B_cell_subset_flu_mouse.integrated <- FindClusters(B_cell_subset_flu_mouse.integrated, resolution= 0.6)
table(B_cell_subset_flu_mouse.integrated@meta.data$seurat_clusters)
metadata <- B_cell_subset_flu_mouse.integrated@meta.data
B_cell_subset_flu_mouse.cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(B_cell_subset_flu_mouse.cell_cluster,'./results/TableS3.csv',row.names = F)

plotS11a <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "umap", label = TRUE)+ 
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15))
plotS11b <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "umap",  group.by = "orig.ident")+ 
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15))
plotS11<- plotS11a + plotS11b
print(plotS11)
ggsave("./results/FigureS11.pdf", plot = plotS11, width = 12, height = 6)
ggsave("./results/FigureS11.png", plot = plotS11, width = 12, height = 6)

B_cell_subset_flu_mouse.integrated <- JoinLayers(B_cell_subset_flu_mouse.integrated, assay = "RNA")
dim(B_cell_subset_flu_mouse.integrated[["RNA"]]$data) 
B_cell_subset_flu_mouse.integrated_markers.genes <- FindAllMarkers(B_cell_subset_flu_mouse.integrated, logfc.threshold = 0.25, test.use = "wilcox", assay = "RNA",
                                      min.pct = 0.25,  only.pos = FALSE)
saveRDS(B_cell_subset_flu_mouse.integrated, "./results/B_cell_subset_flu_mouse.integrated.rds")
saveRDS(B_cell_subset_flu_mouse.integrated_markers.genes, "./results/B_cell_subset_flu_mouse.integrated_markers.genes.rds")
B_cell_subset_flu_mouse.integrated <- readRDS(
  "./results/B_cell_subset_flu_mouse.integrated.rds"
)
#通过查文献marker进行注释
#B_MZ = c("Cr2","Cd1d1","S1pr3")#
#B_RC = c("Fcrl5","Zbtb20","Ptpn22")#
#B_GC_Pre = c("Eif4a1","Mif","Ran","Eif5a","Npm1")#
#B_GC_LZ = c("Cd83","Cd69")#
#B_GC_DZ = c("Rps15a","Rps26","Rps20")#
#B_GC = c("Hmgb2","Tuba1b","Top2a","Tubb5")#
#B_PC = c("Sdc1","Slpi","Jchain","Xbp1","Slc3a2","Ly6a","Mzb1")#
#B_Native = c("Ccr7","Cd79a","Fcmr")#
#B_MC_Pre = c("Bach2","Ebf1")#
#B_MC = c("Spib","Vpreb3","Fam129c")#
# cluster14_markers <- B_cell_subset_flu_mouse.integrated_markers.genes %>%
#   filter(cluster == 14)
# head(cluster14_markers)
# cluster16_markers <- B_cell_subset_flu_mouse.integrated_markers.genes %>%
#   filter(cluster == 16)
# head(cluster16_markers)
Bcell_genes_all1 = c( "Cr2","Cd1d1","Cd9",#MZ
                      #"Ccr7","Ebf1","Btg1","Cd79a",#Native
                      #"Cd69","Cd83","Nr4a1",#Native/Activated
                "Eif4a1","Mif","Ran","Eif5a","Npm1",#PreGC
                "Mki67","Hmgb2", "Tuba1b", "Top2a","Tubb5",#GC
                "Sdc1","Slpi","Jchain","Prdm1","Xbp1","Slc3a2","Ly6a","Mzb1",#PB
                "Bach2", "Spib","Cd19","Vpreb3","Fam129c")#Bmem

#mz:0,7
#GC:1,4,11,13,14,19,20,21,23,24
#Bmem:2,22
#PB:3,5,6,8,9,10,12,15,16,17,18

plotS12 <- DotPlot(B_cell_subset_flu_mouse.integrated, features = Bcell_genes_all1)+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
print(plotS12)
ggsave("./results/FigureS12.pdf", plot = plotS12, width = 18, height = 9)
ggsave("./results/FigureS12.png", plot = plotS12, width = 18, height = 9)

##### Figure1d: B cell heatmap ------
table(Idents(B_cell_subset_flu_mouse.integrated))
Idents(B_cell_subset_flu_mouse.integrated) <- "seurat_clusters"
B_cell_subset_flu_mouse.integrated <- RenameIdents(object = B_cell_subset_flu_mouse.integrated,
                                   "0" = "MZ",
                                   "7" = "MZ",
                                   "1" = "GC",
                                   "4" = "GC",
                                   "11" = "GC",
                                   "13" = "GC",
                                   "14" = "GC",
                                   "19" = "GC",
                                   "20" = "GC",
                                   "21" = "GC",
                                   "23" = "GC",
                                   "24" = "GC",
                                   "3" = "PB",
                                   "5" = "PB",
                                   "6" = "PB",
                                   "8" = "PB",
                                   "9" = "PB",
                                   "10" = "PB",
                                   "12" = "PB",
                                   "15" = "PB",
                                   "16" = "PB",
                                   "17" = "PB",
                                   "18" = "PB",
                                   "2" = "Bmem",
                                   "22" = "Bmem")
table(Idents(B_cell_subset_flu_mouse.integrated))
#DefaultAssay(flu_mouse.combined) <- "RNA"
B_cell_subset_flu_mouse.integrated <- JoinLayers(B_cell_subset_flu_mouse.integrated, assay = "RNA")
dim(B_cell_subset_flu_mouse.integrated[["RNA"]]$data) 
B_cell_subset_flu_mouse.integrated_markers.gene <- FindAllMarkers(B_cell_subset_flu_mouse.integrated, logfc.threshold = 0.25, test.use = "wilcox", assay = "RNA",
                                                  min.pct = 0.1, only.pos =  TRUE)
top15_c <- B_cell_subset_flu_mouse.integrated_markers.gene %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

DoHeatmap(B_cell_subset_flu_mouse.integrated,top15_c$gene,size=8)

Idents(B_cell_subset_flu_mouse.integrated) <- factor(
  Idents(B_cell_subset_flu_mouse.integrated),
  levels = c("MZ","GC","PB","Bmem")
)
id.cols <- c(MZ="#e57373", GC="#8bc34a", PB="#26c6da", Bmem="#d59cf5")

p <- DoHeatmap(
  object      = B_cell_subset_flu_mouse.integrated,
  features    = unique(top15_c$gene),
  size        = 8,
  draw.lines  = FALSE,
  group.colors= id.cols
) +
  scale_fill_gradientn(
    colors = c("#6fb2e4","#eee461","#c66526"),
    name   = "Expression"
  ) +
  theme(
    axis.text.x  = element_text(family="Arial", size=20),
    axis.title.x = element_text(family="Arial", size=20),
    axis.text.y  = element_text(family="Arial", size=20),
    legend.title = element_text(family="Arial", size=20),
    legend.text  = element_text(family="Arial", size=20),
    legend.key.height = unit(0.8,"cm")
  )

legend_df <- data.frame(
  x = 1:length(id.cols),
  y = 1,
  Identity = factor(names(id.cols), levels = names(id.cols))
)

plot1c <- p +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = Identity),
    inherit.aes = FALSE, size = 0  # 大小设为0，不会在图上看到
  ) +
  scale_color_manual(
    values = id.cols,
    breaks = names(id.cols),
    drop   = FALSE,
    name   = "Identity"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),  # 让图例点看得见
    fill  = guide_colorbar(barwidth = 1, barheight = 10, title.position = "top")
  )

print(plot1c)

ggsave(filename="./results/Figure1c.pdf", plot = plot1c, width = 9, height = 16)
ggsave(filename="./results/Figure1c.png", plot = plot1c, width = 9, height = 16)
#mz:0,7
#GC:1,4,11,13,14,19,20,21,23,24
#Bmem:2,22
#PB:3,5,6,8,9,10,12,15,16,17,18

##### Figure1a: B cell umap ------
#B_cell_subset_flu_mouse.integrated@meta.data$B_cell_subpopulations <- NULL
B_cell_subset_flu_mouse.integrated@meta.data$B_cell_subpopulations <- plyr::mapvalues(
                                                                           from = c (0,7,1,4,11,13,14,19,20,21,23,24,2,22,3,5,6,8,9,10,12,15,16,17,18),
                                                                           to = c(rep("MZ",2), 
                                                                                  rep("GC",10),
                                                                                  rep("Bmem",2),
                                                                                  rep("PB",11)
                                                                                 ),
                                                                           x = B_cell_subset_flu_mouse.integrated@meta.data$seurat_clusters)

head(B_cell_subset_flu_mouse.integrated@meta.data)
# (可选) 验证结果
table(B_cell_subset_flu_mouse.integrated@meta.data$B_cell_subpopulations)


umap_data <- data.frame(
  UMAP_1 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 2],
  cell_type = B_cell_subset_flu_mouse.integrated$B_cell_subpopulations # 假设这是你的注释列
)

cell_type_medians <- umap_data %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

allcolour <- c("#3A8EBA","#E69F00","#009E73","#D55E00", "#999999", "#CC79A7","#CC8897")

axis_origin_x <- min(umap_data$UMAP_1) - 1
axis_origin_y <- min(umap_data$UMAP_2) - 1
arrow_length <- 4 

plot1a <- ggplot() +
  geom_point(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type), size = 1) +
  geom_label_repel(
    data = cell_type_medians, 
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    fontface = "bold",
    size = 7,
    box.padding = 0.5,
    point.padding = 0.5,
    label.size = 0.5,
    segment.color = 'grey50',
    family = "Arial" 
  ) +
  
  geom_segment(
    aes(x = axis_origin_x, y = axis_origin_y, xend = axis_origin_x + arrow_length, yend = axis_origin_y),
    arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1
  ) +
  geom_segment(
    aes(x = axis_origin_x, y = axis_origin_y, xend = axis_origin_x, yend = axis_origin_y + arrow_length),
    arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1
  ) +
  
  annotate("text", x = axis_origin_x + arrow_length / 2, y = axis_origin_y - 0.5, label = "UMAP 1", size = 6, fontface = "bold") +
  annotate("text", x = axis_origin_x - 0.5, y = axis_origin_y + arrow_length / 2, label = "UMAP 2", size = 6, fontface = "bold", angle = 90) +
  scale_color_manual(values = allcolour) +
  theme_void() +
  theme(
    legend.position = "none", # 因为有标签，所以不需要图例
    plot.background = element_rect(fill = "white", color = NA) # 确保背景是白色无边框
  )
print(plot1a)
ggsave("./results/Figure1a.pdf", plot = plot1a, width = 9, height = 8)
ggsave("./results/Figure1a.png", plot = plot1a, width = 9, height = 8, dpi = 300)

##### Figure1b: Cell counts and percentage, Antigenic stimulation drives robust expansion of antigen-experienced B cell populations ------------------------------------------------

meta_df <- B_cell_subset_flu_mouse.integrated@meta.data
write.csv(meta_df,'./results/TableS4.csv',row.names = F)

plot_data_grouped <- meta_df %>%
  # 在管道的最开始就进行 unlist，确保后续操作的输入是干净的
  mutate(
    orig.ident = as.character(unlist(orig.ident)),
    B_cell_subpopulations = as.character(unlist(B_cell_subpopulations))
  ) %>%
  mutate(
    experiment_group = case_when(
      grepl("naive", orig.ident, ignore.case = TRUE) ~ "Naive",
      grepl("H5", orig.ident, ignore.case = TRUE)    ~ "H5N1",
      grepl("H1", orig.ident, ignore.case = TRUE)    ~ "H1N1",
      TRUE                                           ~ "Other"
    )
  ) %>%
  # 明确使用 dplyr::count
  dplyr::count(experiment_group, B_cell_subpopulations, name = "n")

total_counts_grouped <- plot_data_grouped %>%
  group_by(experiment_group) %>%
  summarise(total_n = sum(n))

cell_type_colors <- c(
  "MZ" = "#3A8EBA", 
  "GC" = "#E69F00", 
  "PB" = "#D55E00", 
  "Bmem" = "#009E73")

plot1b <- ggplot(plot_data_grouped, aes(x = experiment_group, y = n, fill = B_cell_subpopulations)) +
  geom_col(position = "fill") +

  geom_text(
    data = total_counts_grouped,
    aes(x = experiment_group, y = 1, label = paste("n =", total_n)), 
    hjust = -0.1, 
    size = 6, # 可以适当增大字体
    fontface = "bold",
    inherit.aes = FALSE 
  ) +
  
  coord_flip() +
  scale_fill_manual(values = cell_type_colors) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, .25))) + # 增大了扩展空间

  labs(
    title = "B Cell Proportions by Group",
    x = "Experimental Group",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 18, face = "bold"), # Y轴标签（组名）可以加粗
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )
print(plot1b)
ggsave("./results/Figure1b.pdf", plot = plot1b, width = 9, height = 4.5)
ggsave("./results/Figure1b.png", plot = plot1b, width = 9, height = 4.5, dpi = 300)


##### Figure1d: Dotplot, identification of distinct B cell functional subsets -----------------------------------------------
Bcell_genes_all1 = c( "Cr2","Cd1d1","Cd9",#MZ
                      #"Ccr7","Ebf1","Btg1","Cd79a",#Native
                      #"Cd69","Cd83","Nr4a1",#Native/Activated
                      "Eif4a1","Mif","Ran","Eif5a","Npm1",#PreGC
                      "Mki67","Hmgb2", "Tuba1b", "Top2a","Tubb5",#GC
                      "Bach2", "Spib","Cd19","Vpreb3","Fam129c",#Bmem
                      "Sdc1","Slpi","Jchain","Prdm1","Xbp1","Slc3a2","Ly6a","Mzb1")#PB

plot1d <- DotPlot(B_cell_subset_flu_mouse.integrated, features = Bcell_genes_all1,group.by = "B_cell_subpopulations") +
   scale_color_gradientn(
    colors = c("#313695", "#FFFFBF", "#A50026"), # 低值(蓝) - 中值(白) - 高值(红)
    name = "Average\nExpression" # 自定义图例标题
  ) +
  theme(
    axis.text.x = element_text(family = "Arial", angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(family = "Arial", size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(family = "Arial"),
    legend.text = element_text(family = "Arial", size = 20),
    legend.title = element_text(family = "Arial", size = 20, face = "bold")
  ) +
  guides(
    color = guide_colorbar(
      title = "Expression", 
      title.theme = element_text(family = "Arial", size = 20),
      label.theme = element_text(family = "Arial", size = 20)
    ),
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(family = "Arial", size = 20),
      label.theme = element_text(family = "Arial", size = 20)
    )
  )
print(plot1d)
ggsave("./results/Figure1d.pdf", plot = plot1d, width = 18, height = 6)
ggsave("./results/Figure1d.png", plot = plot1d, width = 18, height = 6)

##### Figure1e: Split UMAPs,Naive、H1N1、H5N1 -----------------------------------------------

obj <- B_cell_subset_flu_mouse.integrated

obj$group_name <- dplyr::case_when(
  grepl("^naive[_-]?mouse", obj$orig.ident, ignore.case = TRUE) ~ "naive",
  grepl("^flu[_-]?h1",      obj$orig.ident, ignore.case = TRUE) ~ "flu_H1",
  grepl("^flu[_-]?h5",      obj$orig.ident, ignore.case = TRUE) ~ "flu_H5",
  TRUE ~ "other"   # 兜底，避免 NA；也方便你检查是否还有未覆盖的样本名
)

# 设定想要的展示顺序（如果不需要“other”，后面可以再丢弃）
obj$group_name <- factor(obj$group_name, levels = c("naive","flu_H1","flu_H5","other"))

# 复查
table(obj$group_name, useNA = "ifany")

umap_data <- data.frame(
  UMAP_1 = obj@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = obj@reductions$umap@cell.embeddings[, 2],
  cell_type = obj$B_cell_subpopulations,
  group = obj$group_name # <-- 使用我们刚刚创建的新列
)

# 分组计算标签位置
cell_type_medians <- umap_data %>%
  group_by(group, cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = 'drop'
  )

# 颜色设置 (使用命名向量是一个好习惯，可以确保颜色和细胞类型正确对应)
allcolour <- c("MZ" = "#3A8EBA", "GC" = "#E69F00", "Bmem" = "#009E73", "PB" = "#D55E00")


# --- 步骤 3: 绘制最终的分面 UMAP 图 ---
plot1e <- ggplot() +
  # 使用背景点来显示所有细胞的轮廓，让每个分面图的形状保持一致
  geom_point(data = umap_data[, c("UMAP_1", "UMAP_2")], aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  # 在背景之上，绘制每个组的有色细胞点
  geom_point(data = umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type), size = 0.5) +
  geom_label_repel(
    data = cell_type_medians, 
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    fontface = "bold",
    size = 8,
    box.padding = 0.3,
    point.padding = 0.3,
    label.size = 0.3,
    segment.color = 'grey50',
    family = "Arial" 
  ) +
  scale_color_manual(values = allcolour, name = "Cell Type") + # 添加图例标题
  facet_wrap(~ group) + # 核心分面代码，现在会按 "naive", "flu_H1", "flu_H5" 分面
  theme_void() + # 使用一个简洁的主题
  theme(
    legend.position = "right", # 将图例放在右侧
    plot.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 20, face = "bold", family = "Arial"), # 分面标题样式
    panel.border = element_rect(colour = "black", fill=NA, size=1), # 为每个分面添加边框
    legend.title = element_text(size = 18, face = "bold", family = "Arial"),  # 增大图例标题字体
    legend.text = element_text(size = 18, family = "Arial"),  # 增大图例文本的字体
    legend.key.size = unit(1.5, "cm"),  # 增大图例符号的大小
    legend.key.height = unit(1.2, "cm"),  # 调整图例项的高度
    legend.key.width = unit(1.2, "cm")  # 调整图例项的宽度
     ) +
  guides(color = guide_legend(override.aes = list(size = 8))) # 让图例中的点变大，更易读

# 打印最终的图像
print(plot1e)
ggsave("./results/Figure1e.pdf", plot = plot1e, width = 18, height = 6)
ggsave("./results/Figure1e.png", plot = plot1e, width = 18, height = 6)

##### Figure1f: Pseudotime -----------------------------------------------
# install.packages("sf", dependencies = TRUE)
# install.packages("spdep")
# library(sf)
# library(spdep)
# Sys.setenv(PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/lib/pkgconfig:/usr/share/pkgconfig")
# Sys.setenv(H5CC="/usr/bin/h5cc")
# remotes::install_github("bnprks/BPCells/r")
# 
# remotes::install_github("cole-trapnell-lab/monocle3")
# remotes::install_github("satijalab/seurat-wrappers")
library(BPCells)
packageVersion("BPCells")
library(monocle3)
packageVersion("monocle3")
library(SeuratWrappers)
packageVersion("SeuratWrappers")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                      'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                      'SummarizedExperiment', 'batchelor', 'HDF5Array',
                      'ggrastr'))



#构建 Monocle3 CellDataSet
mat <- GetAssayData(B_cell_subset_flu_mouse.integrated, slot = "data")  
cellInfo <- B_cell_subset_flu_mouse.integrated@meta.data  
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))  

cds <- new_cell_data_set(expression_data = mat,  
                         cell_metadata = cellInfo,  
                         gene_metadata = geneInfo)  

cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(B_cell_subset_flu_mouse.integrated, "umap")

# 轨迹学习
cds = cluster_cells(cds, cluster_method = 'louvain')
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len=10
))

#设定轨迹起始点
start = c("Cluster_A")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name

flag = closest_vertex[as.character(colData(cds)$T_cell_annotation) %in% start,]
flag = as.numeric(names(which.max(table( flag ))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds)

#绘制轨迹图
plot1f = plot_cells(cds,
                    color_cells_by = "pseudotime",
                    label_cell_groups=F,
                    label_groups_by_cluster=F,
                    label_roots=F,
                    label_leaves=F,
                    label_branch_points=F,
                    cell_size=0.8,
                    group_label_size=3,
                    rasterize=F)+ 
  # 专业主题 - 增加坐标轴字体大小
  theme_classic(base_family = "Arial", base_size = 20) +
  theme(
    # 坐标轴线和刻度
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    
    # ★★★ 坐标轴标题 ★★★
    axis.title.x = element_text(  # 横坐标标题
      size = 20, 
      face = "bold", 
      color = "black",
      margin = margin(t = 15)  # 增加标题上方间距
    ),
    axis.title.y = element_text(  # 纵坐标标题
      size = 20,
      face = "bold", 
      color = "black",
      margin = margin(r = 15)   # 增加标题右侧间距
    ),
    
    # ★★★ 坐标轴刻度文字 ★★★
    axis.text.x = element_text(  # 横坐标刻度
      size = 20, 
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(  # 纵坐标刻度
      size = 20, 
      color = "black",
      face = "bold"
    )
  )

#数据导出与保存
time = pseudotime(cds)
time[time==Inf] = 0
time = data.frame(cell=names(time), pseudotime=time)
write.table(time, file="./results/pseudotime_data.txt", sep="\t", quote=F, col.names=T, row.names=F)
saveRDS(cds, file="./results/monocle3_result.rds")
ggsave(filename="./results/Figure1f.pdf", plot = plot1f, width = 10, height = 8)
ggsave(filename="./results/Figure1f.png", plot = plot1f, width = 10, height = 8)

##### Figure1g: DEG for PB-----------------------------------------------

###1
PB_cells <- subset(B_cell_subset_flu_mouse.integrated, 
                   subset = B_cell_subpopulations == "PB")
PB_group_names <- case_when(
  PB_cells$orig.ident %in% c("naive_mouse1", "naive_mouse2") ~ "naive", # 理论上PB_cells里不会有naive组，但写上更稳健
  PB_cells$orig.ident %in% c("flu_H1_mouse1", "flu_H1_mouse2") ~ "flu_H1",
  PB_cells$orig.ident == "flu_H5_mouse" ~ "flu_H5"
)
PB_cells$group_name <- PB_group_names
head(PB_cells@meta.data)
table(PB_cells$group_name)
#
# 确保您的分组信息存储在 'group_name' 列中
# logfc.threshold = 0.25 和 min.pct = 0.1 是常用的过滤参数
deg_PB_standard <- FindMarkers(PB_cells, 
                               ident.1 = "flu_H5", 
                               ident.2 = "flu_H1", 
                               group.by = "group_name",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)

# 查看结果
head(deg_PB_standard)

###2
# 1. 获取两组的细胞ID和数量
h1_PB_cells <- WhichCells(PB_cells, expression = group_name == "flu_H1")
h5_PB_cells <- WhichCells(PB_cells, expression = group_name == "flu_H5")

n_h1 <- length(h1_PB_cells)
n_h5 <- length(h5_PB_cells)

target_n <- min(n_h1, n_h5) # 在您的例子中，这里会是约600

# 2. 进行多次降采样和差异分析
n_iterations <- 50 # 建议至少50次，100次更佳
PB_deg_results_list <- list() # 这次我们保存完整的数据框

for (i in 1:n_iterations) {
  set.seed(i)
  h5_sampled_cells <- sample(h5_PB_cells, size = target_n)
  combined_cells <- c(h1_PB_cells, h5_sampled_cells)
  temp_subset <- subset(PB_cells, cells = combined_cells)
  
  # 我们运行FindMarkers，但这次获取所有基因的结果，而不仅仅是显著的
  # 将 logfc.threshold 和 min.pct 设为0，以返回所有测试过的基因
  PB_deg_temp <- FindMarkers(temp_subset,
                          ident.1 = "flu_H5",
                          ident.2 = "flu_H1",
                          group.by = "group_name",
                          logfc.threshold = 0, # 获取所有基因的FC
                          min.pct = 0,         # 获取所有基因
                          verbose = FALSE)     # 关闭过程打印，让循环更快
  
  # 将基因名作为一列，方便后续处理
  PB_deg_temp$gene <- rownames(PB_deg_temp)
  
  # 存储每次的完整结果数据框
  PB_deg_results_list[[i]] <- PB_deg_temp
  
  print(paste("Completed iteration:", i))
}

### **第二步：聚合所有迭代的结果**

library(dplyr)
library(tibble)

# 1. 将列表中的所有数据框合并成一个大的数据框
PB_all_iterations_df <- bind_rows(PB_deg_results_list)

# 2. 按基因名分组，并计算聚合统计量
# 我们计算 avg_log2FC 的均值，以及 p_val_adj 的均值
# 我们还计算一个基因在多少次迭代中是显著的，作为稳健性的衡量
PB_aggregated_results <- PB_all_iterations_df %>%
  group_by(gene) %>%
  summarise(
    mean_avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    mean_p_val_adj = mean(p_val_adj, na.rm = TRUE),
    # 计算一个基因在多少次迭代中是显著的 (例如 p.adj < 0.05)
    n_significant = sum(p_val_adj < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# 3. 准备用于火山图的数据
# 计算 -log10 of the mean adjusted p-value
PB_volcano_data <- PB_aggregated_results %>%
  mutate(
    neg_log10_padj = -log10(mean_p_val_adj)
  )

# 查看最终的聚合结果
head(PB_volcano_data)

library(ggplot2)
library(ggrepel)

# 1. 设置筛选阈值
log2FC_threshold <- 1  # log2FC > 0.5 或 < -0.5
padj_threshold <- 0.05   # 对应 -log10(0.05) 约为 1.3

# 2. 添加一列来标记基因是上调、下调还是不显著
PB_volcano_data <- PB_volcano_data %>%
  mutate(
    significance = case_when(
      mean_avg_log2FC > log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Upregulated in H5",
      mean_avg_log2FC < -log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Downregulated in H5",
      TRUE ~ "Not Significant"
    )
  )

# (可选) 查看各类基因的数量
table(PB_volcano_data$significance)

# 3. 筛选出我们想要在图上标记的基因
genes_to_label <- PB_volcano_data %>%
  filter(abs(mean_avg_log2FC) > log2FC_threshold & mean_p_val_adj < padj_threshold) %>%
  arrange(mean_p_val_adj) %>%
  head(30)

# 4. 绘制火山图
plot1g <- ggplot(PB_volcano_data, aes(x = mean_avg_log2FC, y = neg_log10_padj)) +
  # 绘制所有基因的点
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  
  # 设置颜色
  scale_color_manual(values = c("Upregulated in H5" = "#D55E00", 
                                "Downregulated in H5" = "#3A8EBA", 
                                "Not Significant" = "grey80")) +
  
  # 添加阈值线
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey50") +
  
  # 使用 ggrepel 添加基因标签，避免重叠
  geom_text_repel(data = genes_to_label, 
                  aes(label = gene),
                  size = 4,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  max.overlaps = Inf) + # 尽可能多地显示标签
  labs(
    title = "Volcano Plot of PB cells: H5 vs H1",
    x = "Mean Log2 Fold Change",
    y = "-log10(Mean Adjusted P-value)",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial"),
    axis.title.x = element_text(size = 18, family = "Arial"),  # x 轴标题字体变大
    axis.title.y = element_text(size = 18, family = "Arial"),  # y 轴标题字体变大
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # 图例标题字体大且为Arial
    legend.text = element_text(size = 16, family = "Arial"),  # 图例文本字体为Arial且较大
    legend.margin = margin(t = 8),  
    text = element_text(size = 18, family = "Arial")  # 默认字体为Arial
  )

print(plot1g)
ggsave(filename="./results/Figure1g.pdf", plot = plot1g, width = 8.5, height = 7)
ggsave(filename="./results/Figure1g.png", plot = plot1g, width = 8.5, height = 7)

##### Figure1h: DEG for GC-----------------------------------------------
# MZ_cells <- subset(B_cell_subset_flu_mouse.integrated, 
#                    subset = B_cell_subpopulations == "MZ")
# MZ_group_names <- case_when(
#   MZ_cells$orig.ident %in% c("naive_mouse1", "naive_mouse2") ~ "naive", # 理论上PB_cells里不会有naive组，但写上更稳健
#   MZ_cells$orig.ident %in% c("flu_H1_mouse1", "flu_H1_mouse2") ~ "flu_H1",
#   MZ_cells$orig.ident == "flu_H5_mouse" ~ "flu_H5"
# )
# MZ_cells$group_name <- MZ_group_names
# head(MZ_cells@meta.data)
# table(MZ_cells$group_name)
# #
# # 确保您的分组信息存储在 'group_name' 列中
# # logfc.threshold = 0.25 和 min.pct = 0.1 是常用的过滤参数
# deg_MZ_standard <- FindMarkers(MZ_cells, 
#                                ident.1 = "flu_H5", 
#                                ident.2 = "flu_H1", 
#                                group.by = "group_name",
#                                logfc.threshold = 0.25,
#                                min.pct = 0.1)
# 
# # 查看结果
# head(deg_MZ_standard)
# 
# 
# ###2
# # 1. 获取两组的细胞ID和数量
# h1_MZ_cells <- WhichCells(MZ_cells, expression = group_name == "flu_H1")
# h5_MZ_cells <- WhichCells(MZ_cells, expression = group_name == "flu_H5")
# 
# n_h1 <- length(h1_MZ_cells)
# n_h5 <- length(h5_MZ_cells)
# 
# target_n <- min(n_h1, n_h5) # 在您的例子中，这里会是约600
# 
# # 2. 进行多次降采样和差异分析
# n_iterations <- 50 # 建议至少50次，100次更佳
# MZ_deg_results_list <- list() # 这次我们保存完整的数据框
# 
# for (i in 1:n_iterations) {
#   set.seed(i)
#   h5_sampled_cells <- sample(h5_MZ_cells, size = target_n)
#   combined_cells <- c(h1_MZ_cells, h5_sampled_cells)
#   temp_subset <- subset(MZ_cells, cells = combined_cells)
#   
#   # 我们运行FindMarkers，但这次获取所有基因的结果，而不仅仅是显著的
#   # 将 logfc.threshold 和 min.pct 设为0，以返回所有测试过的基因
#   MZ_deg_temp <- FindMarkers(temp_subset,
#                           ident.1 = "flu_H5",
#                           ident.2 = "flu_H1",
#                           group.by = "group_name",
#                           logfc.threshold = 0, # 获取所有基因的FC
#                           min.pct = 0,         # 获取所有基因
#                           verbose = FALSE)     # 关闭过程打印，让循环更快
#   
#   # 将基因名作为一列，方便后续处理
#   MZ_deg_temp$gene <- rownames(MZ_deg_temp)
#   
#   # 存储每次的完整结果数据框
#   MZ_deg_results_list[[i]] <- MZ_deg_temp
#   
#   print(paste("Completed iteration:", i))
# }
# ### **第二步：聚合所有迭代的结果**
# 
# library(dplyr)
# library(tibble)
# 
# # 1. 将列表中的所有数据框合并成一个大的数据框
# MZ_all_iterations_df <- bind_rows(MZ_deg_results_list)
# 
# # 2. 按基因名分组，并计算聚合统计量
# # 我们计算 avg_log2FC 的均值，以及 p_val_adj 的均值
# # 我们还计算一个基因在多少次迭代中是显著的，作为稳健性的衡量
# MZ_aggregated_results <- MZ_all_iterations_df %>%
#   group_by(gene) %>%
#   summarise(
#     mean_avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
#     mean_p_val_adj = mean(p_val_adj, na.rm = TRUE),
#     # 计算一个基因在多少次迭代中是显著的 (例如 p.adj < 0.05)
#     n_significant = sum(p_val_adj < 0.05, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# # 3. 准备用于火山图的数据
# # 计算 -log10 of the mean adjusted p-value
# MZ_volcano_data <- MZ_aggregated_results %>%
#   mutate(
#     neg_log10_padj = -log10(mean_p_val_adj)
#   )
# 
# # 查看最终的聚合结果
# head(MZ_volcano_data)
# 
# library(ggplot2)
# library(ggrepel)
# 
# # 1. 设置筛选阈值
# log2FC_threshold <- 1  # log2FC > 0.5 或 < -0.5
# padj_threshold <- 0.05   # 对应 -log10(0.05) 约为 1.3
# 
# # 2. 添加一列来标记基因是上调、下调还是不显著
# MZ_volcano_data <- MZ_volcano_data %>%
#   mutate(
#     significance = case_when(
#       mean_avg_log2FC > log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Upregulated in H5",
#       mean_avg_log2FC < -log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Downregulated in H5",
#       TRUE ~ "Not Significant"
#     )
#   )
# 
# # (可选) 查看各类基因的数量
# table(MZ_volcano_data$significance)
# 
# # 3. 筛选出我们想要在图上标记的基因
# MZ_genes_to_label <- MZ_volcano_data %>%
#   filter(abs(mean_avg_log2FC) > log2FC_threshold & mean_p_val_adj < padj_threshold) %>%
#   arrange(mean_p_val_adj) %>%
#   head(30)
# 
# # 4. 绘制火山图
# plot1h1 <- ggplot(MZ_volcano_data, aes(x = mean_avg_log2FC, y = neg_log10_padj)) +
#   # 绘制所有基因的点
#   geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
#   
#   # 设置颜色
#   scale_color_manual(values = c(
#     "Upregulated in H5" = "#D55E00", 
#     "Downregulated in H5" = "#3A8EBA", 
#     "Not Significant" = "grey80"
#   )) +
#   
#   # 添加阈值线
#   geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
#              linetype = "dashed", color = "grey50") +
#   geom_hline(yintercept = -log10(padj_threshold),
#              linetype = "dashed", color = "grey50") +
#   
#   # 使用 ggrepel 添加基因标签（✅ 修复了 x / y 缺失的问题）
#   geom_text_repel(
#     data = genes_to_label, 
#     aes(
#       x = mean_avg_log2FC, 
#       y = neg_log10_padj, 
#       label = gene
#     ),
#     size = 4,
#     box.padding = 0.5,
#     point.padding = 0.5,
#     max.overlaps = Inf
#   ) +
#   
#   labs(
#     title = "Volcano Plot of GC cells: H5 vs H1",
#     x = "Mean Log2 Fold Change",
#     y = "-log10(Mean Adjusted P-value)",
#     color = "Significance"
#   ) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial"),
#     axis.title.x = element_text(size = 18, family = "Arial"),  # x 轴标题字体变大
#     axis.title.y = element_text(size = 18, family = "Arial"),  # y 轴标题字体变大
#     legend.position = "top",
#     legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # 图例标题字体大且为Arial
#     legend.text = element_text(size = 16, family = "Arial"),  # 图例文本字体为Arial且较大
#     legend.margin = margin(t = 8),  
#     text = element_text(size = 18, family = "Arial")  # 默认字体为Arial
#   )
# 
# print(plot1h1)
# ggsave(filename="./results/Figure1h1.pdf", plot = plot1h, width = 8.5, height = 7)
# ggsave(filename="./results/Figure1h1.png", plot = plot1h, width = 8.5, height = 7)
GC_cells <- subset(B_cell_subset_flu_mouse.integrated, 
                   subset = B_cell_subpopulations == "GC")
GC_group_names <- case_when(
  GC_cells$orig.ident %in% c("naive_mouse1", "naive_mouse2") ~ "naive", # 理论上PB_cells里不会有naive组，但写上更稳健
  GC_cells$orig.ident %in% c("flu_H1_mouse1", "flu_H1_mouse2") ~ "flu_H1",
  GC_cells$orig.ident == "flu_H5_mouse" ~ "flu_H5"
)
GC_cells$group_name <- GC_group_names
head(GC_cells@meta.data)
table(GC_cells$group_name)
#
# 确保您的分组信息存储在 'group_name' 列中
# logfc.threshold = 0.25 和 min.pct = 0.1 是常用的过滤参数
deg_GC_standard <- FindMarkers(GC_cells, 
                               ident.1 = "flu_H5", 
                               ident.2 = "flu_H1", 
                               group.by = "group_name",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)

# 查看结果
head(deg_GC_standard)


###2
# 1. 获取两组的细胞ID和数量
h1_GC_cells <- WhichCells(GC_cells, expression = group_name == "flu_H1")
h5_GC_cells <- WhichCells(GC_cells, expression = group_name == "flu_H5")

n_h1 <- length(h1_GC_cells)
n_h5 <- length(h5_GC_cells)

target_n <- min(n_h1, n_h5) # 在您的例子中，这里会是约600

# 2. 进行多次降采样和差异分析
n_iterations <- 50 # 建议至少50次，100次更佳
GC_deg_results_list <- list() # 这次我们保存完整的数据框

for (i in 1:n_iterations) {
  set.seed(i)
  h5_sampled_cells <- sample(h5_GC_cells, size = target_n)
  combined_cells <- c(h1_GC_cells, h5_sampled_cells)
  temp_subset <- subset(GC_cells, cells = combined_cells)
  
  # 我们运行FindMarkers，但这次获取所有基因的结果，而不仅仅是显著的
  # 将 logfc.threshold 和 min.pct 设为0，以返回所有测试过的基因
  GC_deg_temp <- FindMarkers(temp_subset,
                             ident.1 = "flu_H5",
                             ident.2 = "flu_H1",
                             group.by = "group_name",
                             logfc.threshold = 0, # 获取所有基因的FC
                             min.pct = 0,         # 获取所有基因
                             verbose = FALSE)     # 关闭过程打印，让循环更快
  
  # 将基因名作为一列，方便后续处理
  GC_deg_temp$gene <- rownames(GC_deg_temp)
  
  # 存储每次的完整结果数据框
  GC_deg_results_list[[i]] <- GC_deg_temp
  
  print(paste("Completed iteration:", i))
}
### **第二步：聚合所有迭代的结果**

library(dplyr)
library(tibble)

# 1. 将列表中的所有数据框合并成一个大的数据框
GC_all_iterations_df <- bind_rows(GC_deg_results_list)

# 2. 按基因名分组，并计算聚合统计量
# 我们计算 avg_log2FC 的均值，以及 p_val_adj 的均值
# 我们还计算一个基因在多少次迭代中是显著的，作为稳健性的衡量
GC_aggregated_results <- GC_all_iterations_df %>%
  group_by(gene) %>%
  summarise(
    mean_avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    mean_p_val_adj = mean(p_val_adj, na.rm = TRUE),
    # 计算一个基因在多少次迭代中是显著的 (例如 p.adj < 0.05)
    n_significant = sum(p_val_adj < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# 3. 准备用于火山图的数据
# 计算 -log10 of the mean adjusted p-value
GC_volcano_data <- GC_aggregated_results %>%
  mutate(
    neg_log10_padj = -log10(mean_p_val_adj)
  )

# 查看最终的聚合结果
head(GC_volcano_data)

library(ggplot2)
library(ggrepel)

# 1. 设置筛选阈值
log2FC_threshold <- 1  # log2FC > 0.5 或 < -0.5
padj_threshold <- 0.05   # 对应 -log10(0.05) 约为 1.3

# 2. 添加一列来标记基因是上调、下调还是不显著
GC_volcano_data <- GC_volcano_data %>%
  mutate(
    significance = case_when(
      mean_avg_log2FC > log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Upregulated in H5",
      mean_avg_log2FC < -log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Downregulated in H5",
      TRUE ~ "Not Significant"
    )
  )

# (可选) 查看各类基因的数量
table(GC_volcano_data$significance)

# 3. 筛选出我们想要在图上标记的基因
GC_genes_to_label <- GC_volcano_data %>%
  filter(abs(mean_avg_log2FC) > log2FC_threshold & mean_p_val_adj < padj_threshold) %>%
  arrange(mean_p_val_adj) %>%
  head(30)

# 4. 绘制火山图
plot1h <- ggplot(GC_volcano_data, aes(x = mean_avg_log2FC, y = neg_log10_padj)) +
  # 绘制所有基因的点
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  
  # 设置颜色
  scale_color_manual(values = c(
    "Upregulated in H5" = "#D55E00", 
    "Downregulated in H5" = "#3A8EBA", 
    "Not Significant" = "grey80"
  )) +
  
  # 添加阈值线
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_threshold),
             linetype = "dashed", color = "grey50") +
  
  # 使用 ggrepel 添加基因标签（✅ 修复了 x / y 缺失的问题）
  geom_text_repel(
    data = GC_genes_to_label,
    aes(x = mean_avg_log2FC, y = neg_log10_padj, label = gene),
    size = 4, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf
  )+
  labs(
    title = "Volcano Plot of GC cells: H5 vs H1",
    x = "Mean Log2 Fold Change",
    y = "-log10(Mean Adjusted P-value)",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial"),
    axis.title.x = element_text(size = 18, family = "Arial"),  # x 轴标题字体变大
    axis.title.y = element_text(size = 18, family = "Arial"),  # y 轴标题字体变大
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # 图例标题字体大且为Arial
    legend.text = element_text(size = 16, family = "Arial"),  # 图例文本字体为Arial且较大
    legend.margin = margin(t = 8),  
    text = element_text(size = 18, family = "Arial")  # 默认字体为Arial
  )

print(plot1h)
ggsave(filename="./results/Figure1h.pdf", plot = plot1h, width = 10, height = 8)
ggsave(filename="./results/Figure1h.png", plot = plot1h, width = 10, height = 8.5,dpi = 300)

### Figure 3: ------------------------------------------------

##### Figure3a------
obj <- B_cell_subset_flu_mouse.integrated  

# 1) 统一一个分组列 condition
if (!"condition" %in% colnames(obj@meta.data)) {
  obj$condition <- NA_character_
  obj$condition[grepl("^naive_", obj$orig.ident)]   <- "naive"
  obj$condition[grepl("^flu_H1_", obj$orig.ident)]  <- "flu_H1"
  obj$condition[grepl("^flu_H5",  obj$orig.ident)]  <- "flu_H5"
  obj$condition <- factor(obj$condition, levels = c("naive","flu_H1","flu_H5"))
}

# 2) 标记是否含有重/轻链（尽量兼容不同列名）
md <- obj@meta.data
getcol <- function(nm) if (nm %in% colnames(md)) md[[nm]] else rep(NA, nrow(md))

vh <- getcol("v_call_10x");      vl <- getcol("light_v_call_10x")
vh2<- getcol("junction_aa");  vl2<- getcol("light_junction_10x_aa")      # 备选列
prod <- getcol("productive");     inframe <- getcol("vj_in_frame")
clon <- getcol("clone_id")

has_heavy <- (!is.na(vh)  & grepl("^IGH", vh))  | (!is.na(vh2)  & grepl("^IGH", vh2))
has_light <- (!is.na(vl)  & grepl("^IG[KL]", vl))| (!is.na(vl2) & grepl("^IG[KL]", vl2))

obj$BCR_paired <- has_heavy & has_light
# 可选的更严格过滤：
obj$BCR_paired_strict <- obj$BCR_paired &
  (is.na(prod)   | prod) &
  (is.na(inframe)| inframe) &
  (is.na(clon)   | !is.na(clon))

table(obj$BCR_paired, useNA = "ifany")
table(obj$BCR_paired_strict, useNA = "ifany")

#(二)
# 只保留 RNA–BCR 共捕获细胞（任选其一）
obj_bcr <- subset(obj, subset = BCR_paired)           # 宽松
# obj_bcr <- subset(obj, subset = BCR_paired_strict)  # 严格


# 如果没有 UMAP 嵌入，补跑（整合后一般都有）
if (!"umap" %in% SeuratObject::Reductions(obj_bcr)) {
  obj_bcr <- RunUMAP(obj_bcr, dims = 1:30)  # 用你的整合维度
}

# ① 单图上按分组上色
DimPlot(obj_bcr, reduction = "umap", group.by = "condition", pt.size = 0.25)

# ② 三面板对比（同一坐标系，分面展示）
DimPlot(obj_bcr, reduction = "umap", split.by = "condition", group.by = "B_cell_subpopulations", ncol = 3, pt.size = 0.25)


# 加载库
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 定义您的细胞亚群颜色 (来自您的问题)
cell_type_colors <- c(
  "MZ" = "#3A8EBA", 
  "GC" = "#E69F00", 
  "PB" = "#D55E00", 
  "Bmem" = "#009E73"
)
# 假设您的 Seurat 对象名为 obj_bcr
# 从元数据中提取信息并进行分组统计
conditions_ordered <- c("naive", "flu_H1", "flu_H5")
obj_bcr$condition <- factor(obj_bcr$condition, levels = conditions_ordered)
summary_df <- obj_bcr@meta.data %>%
  # 按分组 (condition) 和细胞亚群 (B_cell_subpopulations) 进行分组
  group_by(condition, B_cell_subpopulations) %>%
  # 计算每个小组的细胞数
  summarise(count = n(), .groups = 'drop_last') %>%
  # 在每个 condition 内部计算比例
  mutate(proportion = count / sum(count)) %>%
  # ungroup() 是一个好习惯，防止后续意外的分组操作
  ungroup()

# 查看一下我们生成的统计数据框
print(summary_df)

# 创建一个空列表，用于存储我们生成的组合图块
plot_list <- list()

# 使用我们定义好的顺序进行循环
for (cond in conditions_ordered) {
  
  # --- a. 为当前 condition 生成 UMAP 图 ---
  obj_subset <- subset(obj_bcr, subset = condition == cond)
  # 计算每个群体的细胞数量
  cell_counts <- table(obj_subset$B_cell_subpopulations)
  cell_counts_df <- data.frame(
    B_cell_subpopulations = names(cell_counts),
    count = as.numeric(cell_counts)
  )
  umap_plot <- DimPlot(obj_subset, 
                       reduction = "umap", 
                       group.by = "B_cell_subpopulations", 
                       pt.size = 1,
                       cols = cell_type_colors) +
    ggtitle(cond) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial")) + # 统一字体
    NoLegend()
  
  # --- b. 为当前 condition 生成竖向条形图 ---
  summary_subset <- summary_df %>% filter(condition == cond)
  
  barplot <- ggplot(summary_subset, aes(x = " ", y = proportion, fill = B_cell_subpopulations)) +
    # <-- 1. 移除 coord_flip() 使柱形图竖向展示 -->
    geom_col(width = 0.3) + # 绘制堆叠的竖向条形图
    
    # <-- 2. 调整比例数字的字体和大小 -->
    geom_text(aes(label = scales::percent(proportion, accuracy = 1)),
              position = position_stack(vjust = 0.2), 
              color = "black", 
              fontface = "bold",
              size = 4,         # 增大字体尺寸
              family = "Arial"  # 设置字体为 Arial
    ) +
    
    scale_fill_manual(values = cell_type_colors, name = "B Cell") + # 为图例添加标题
    theme_void() +
    theme(legend.position = "right", # 将图例放在右侧
          legend.text = element_text(family = "Arial", size = 16),
          legend.title = element_text(family = "Arial", size = 16, face = "bold"))
  
  # --- c. 组合 UMAP 和条形图 ---
  # 调整宽度比例，让 UMAP 占 3 份，条形图占 1 份
  combined_panel <- umap_plot + barplot + plot_layout(widths = c(3, 1), heights = c(1, 0.3))+
    theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
  
  plot_list[[cond]] <- combined_panel
}
# 将列表中的所有图块组合成一个大图，按 3 列排列
plot3a <- wrap_plots(plot_list, ncol = 3)+plot_layout(guides = "collect")

# 打印最终的图形
print(plot3a)
ggsave("./results/Figure3a.pdf", plot = plot3a, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure3a.png", plot = plot3a, width = 15, height = 5, dpi = 300, device = "png")

##### Figure3b------
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  # 步骤 1: 创建分组列 (和之前一样)
  mutate(treatment_group = case_when(
    grepl("naive", orig.ident, ignore.case = TRUE) ~ "Naive",
    grepl("H1", orig.ident, ignore.case = TRUE)    ~ "Flu_H1",
    grepl("H5", orig.ident, ignore.case = TRUE)    ~ "Flu_H5",
    TRUE                                          ~ "Other"
  )) %>%
  # 步骤 2: 将其转换为指定顺序的因子 (这是关键优化！)
  mutate(treatment_group = factor(treatment_group, levels = c("Naive", "Flu_H1", "Flu_H5")))

# (强烈推荐) 验证顺序是否正确
cat("实验分组的因子水平顺序 (应为 Naive, Flu_H1, Flu_H5):\n")
print(levels(obj_bcr$treatment_group))

head(obj_bcr@meta.data$c_call)
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(isotype = case_when(
    c_call == "IGHM" ~ "IgM",
    c_call == "IGHD" ~ "IgD",
    c_call == "IGHG1" ~ "IgG1",
    c_call == "IGHG2B" ~ "IgG2B", # 将小鼠的亚型合并
    c_call == "IGHG2C" ~ "IgG2C",
    c_call == "IGHG3" ~ "IgG3",
    c_call == "IGHA" ~ "IgA",   # 如果只有 IGHA
    c_call == "IGHE" ~ "IgE",
    TRUE ~ "Unknown" # 将所有其他情况 (包括 NA) 都归为 "Unknown"
  ))
print(table(obj_bcr$isotype))
isotype_colors <- c(
  "IgM" = "#5A5A5A",    # 深灰色/紫色
  "IgD" = "#0072B2",    # 蓝色
  "IgG1" = "#F0E442",   # 黄色
  "IgG2B" = "#E69F00",   # 橙黄色
  "IgG2C" = "#D55E00",   # 橙色
  "IgG3" = "#CC79A7",   # 粉紫色
  "IgA" = "#d62728",   # 红色
  "IgG3" = "#009E73",    # 绿色
  "Unknown" = "#D3D3D3" # 浅灰色
)
# 使用 DimPlot 并应用我们的自定义设置
plot3b <- DimPlot(
  obj_bcr,
  reduction = "umap",
  group.by = "isotype",           # 告诉 Seurat 按 Isotype 为点着色
  split.by = "treatment_group",   # 关键！告诉 Seurat 按实验分组拆分图形
  cols = isotype_colors,          # 应用我们的自定义颜色
  pt.size = 1,
  ncol = 3                        # 将三个组的图排成一行 (3列)
) +
  # --- 您的 ggplot2 美化代码可以继续使用 ---
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold") # 美化拆分视图的标题 (e.g., "Naive")
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 5),
    ncol = 1,
    title = NULL
  ))+
  labs(title = "") 

# 打印最终的图形
print(plot3b)
ggsave("./results/Figure3b.pdf", plot = plot3b, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure3b.png", plot = plot3b, width = 15, height = 5, dpi = 300, device = "png")

#####figure3c------
naive_mouse.bcr$barcode <- sub("_contig_\\d+", "", naive_mouse.bcr$sequence_id)

library(dplyr)
library(stringr)

## 0) 先给每个细胞打上 treatment_group（可按需调整）
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(treatment_group = case_when(
    str_detect(orig.ident, regex("^naive",  ignore_case = TRUE)) ~ "Naive",
    str_detect(orig.ident, regex("^flu_H1", ignore_case = TRUE)) ~ "Flu_H1",
    str_detect(orig.ident, regex("^flu_H5", ignore_case = TRUE)) ~ "Flu_H5",
    TRUE ~ "Other"
  ))

## 1) meta 侧 join key：直接用 barcode_seurat（与你截图一致，保留 -1）
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(barcode_key = as.character(barcode_seurat))

## 2) BCR（naive）侧 join key：去掉 _contig_数字，保留 -1
##    若一个 barcode 有多条 contig，聚合成一个细胞一个值（这里用 max，可换 mean/median）
bcr_naive_keyed <- naive_mouse.bcr %>%
  mutate(barcode_key = barcode) %>%
  group_by(barcode_key) %>%
  summarise(MU_FREQ_naive = max(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE), .groups = "drop")

## 3) 非破坏式合并到临时列，再只在 naive 行写回；其他组完全不动
tmp <- obj_bcr@meta.data %>%
  left_join(bcr_naive_keyed, by = "barcode_key")

idx_naive <- tmp$treatment_group == "Naive"

tmp <- tmp %>%
  mutate(
    MU_FREQ_HEAVY_TOTAL =
      if_else(idx_naive, coalesce(MU_FREQ_naive, MU_FREQ_HEAVY_TOTAL), MU_FREQ_HEAVY_TOTAL)
  ) %>%
  select(-MU_FREQ_naive)

obj_bcr@meta.data <- tmp

## 4) 合并质量检查（应当只在 naive 里出现更多非 NA）
cat("Non-NA SHM counts by group:\n")
print(obj_bcr@meta.data %>%
        group_by(treatment_group) %>%
        summarise(non_NA = sum(!is.na(MU_FREQ_HEAVY_TOTAL)), .groups = "drop"))

## 5) 看看 naive 组前几行是否正确写入
obj_bcr@meta.data %>%
  filter(treatment_group == "Naive") %>%
  select(orig.ident, barcode_seurat, MU_FREQ_HEAVY_TOTAL) %>%
  slice_head(n = 15)

## 6) 如需汇总各组 SHM（跳过 NA）
shm_summary <- obj_bcr@meta.data %>%
  group_by(treatment_group) %>%
  summarise(
    cell_count = n(),
    min_shm    = min(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE),
    max_shm    = max(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE),
    mean_shm   = mean(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE),
    median_shm = median(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE),
    .groups = "drop"
  )
print(shm_summary)
library(tidyr)
# --- 1. 完整的数据准备流程 (这是最终修正版) ---
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  # 步骤 1a: 创建并排序实验分组列
  mutate(treatment_group = case_when(
    grepl("^naive", orig.ident, ignore.case = TRUE) ~ "Naive",
    grepl("^flu_H1", orig.ident, ignore.case = TRUE) ~ "Flu_H1",
    grepl("^flu_H5", orig.ident, ignore.case = TRUE) ~ "Flu_H5",
    TRUE ~ "Other"
  )) %>%
  mutate(treatment_group = factor(treatment_group, levels = c("Naive", "Flu_H1", "Flu_H5"))) %>%
  
  # 步骤 1b: 清洗 SHM 数据：将 NA 替换为 0
  mutate(MU_FREQ_HEAVY_TOTAL = replace_na(MU_FREQ_HEAVY_TOTAL, 0))

# --- 2. 创建一个绝对唯一的元数据列名 (关键步骤！) ---
# 我们使用 "shm_freq_meta" 这个名字，它几乎不可能与任何基因名冲突
obj_bcr$shm_freq_meta <- obj_bcr$MU_FREQ_HEAVY_TOTAL

# (可选) 验证新列的数据是否正确
summary(obj_bcr$shm_freq_meta)
# 1. 提取 UMAP 坐标
umap_coords <- as.data.frame(Embeddings(obj_bcr, reduction = "umap"))

# 2. 提取我们需要的元数据
metadata_to_plot <- obj_bcr@meta.data %>%
  select(treatment_group, shm_freq_meta) # 只选择我们需要的列

# 3. 将坐标和元数据合并成一个用于绘图的 data.frame
plotting_df <- cbind(umap_coords, metadata_to_plot)

head(plotting_df)
plotting_df <- plotting_df %>%
  arrange(shm_freq_meta)

# 开始绘图
plot3c <- ggplot(
  plotting_df, 
  aes(x = umap_1, y = umap_2, color = shm_freq_meta) # 设置 x, y 轴和颜色映射
) +
  # 绘制散点图
  geom_point(size = 1.0) +
  
  # 关键！使用 facet_wrap 来实现按分组拆分视图
  facet_wrap(~ treatment_group, ncol = 3) +
  
  # 应用与之前完全相同的颜色梯度
  scale_color_gradientn(
    colors = c("grey85", "#FEE0D2", "#DE2D26", "#A50F15"),
    name = "SHM Freq."
  ) +
  
  # 应用主题和美化
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"), # 美化分组标题
    strip.background = element_blank(), # 移除分组标题的背景框
    legend.position = "right", # 将图例放在右侧
    
    # 增大图例字体并将字体设为 Arial
    legend.title = element_text(size = 16, family = "Arial", face = "bold"),  # 设置图例标题字体为 Arial，并增大字体
    legend.text = element_text(size = 16, family = "Arial")  # 设置图例文本字体为 Arial，并增大字体
  )

# 打印最终的图形
print(plot3c)
ggsave("./results/Figure3c.pdf", plot = plot3c, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure3c.png", plot = plot3c, width = 15, height = 5, dpi = 300, device = "png")

#####figure3d------
library(dplyr)
library(tibble)
library(ggplot2)


CLONE_ID_COLUMN <- "clone_id.y"
GROUP_COLUMN    <- "condition"
SUBPOP_COLUMN   <- "B_cell_subpopulations" 
TOP_N           <- 10
MIN_CLONE_SIZE  <- 2    


clean_meta_df <- obj@meta.data %>%
  rownames_to_column("barcode") %>%
  transmute(
    barcode,
    clone_id.y  = as.character(.data[[CLONE_ID_COLUMN]]),
    condition = as.character(.data[[GROUP_COLUMN]]),
    subpop    = as.character(.data[[SUBPOP_COLUMN]])
  )

clone_counts <- clean_meta_df %>%
  filter(!is.na(clone_id.y), !is.na(condition)) %>%
  dplyr::count(condition, clone_id.y, name = "clone_size") %>%
  arrange(condition, desc(clone_size))

top_by_group <- clone_counts %>%
  group_by(condition) %>%
  arrange(desc(clone_size), .by_group = TRUE) %>%
  filter(clone_size >= MIN_CLONE_SIZE) %>%
  slice_head(n = TOP_N) %>%
  ungroup()


missing_groups <- setdiff(unique(clone_counts$condition), unique(top_by_group$condition))
if (length(missing_groups) > 0) {
  top_fallback <- clone_counts %>%
    filter(condition %in% missing_groups) %>%
    group_by(condition) %>%
    arrange(desc(clone_size), .by_group = TRUE) %>%
    slice_head(n = TOP_N) %>%
    ungroup()
  top_by_group <- bind_rows(top_by_group, top_fallback)
}

highlight_df <- clean_meta_df %>%
  left_join(top_by_group %>% select(condition, clone_id.y, clone_size), 
            by = c("condition","clone_id.y")) %>%
  mutate(
    top_clone_highlight = case_when(
      condition == "flu_H1" & !is.na(clone_size) ~ "H1",
      condition == "flu_H5" & !is.na(clone_size) ~ "H5",
      condition == "naive"  & !is.na(clone_size) ~ "Naive",
      TRUE ~ "Other"
    )
  ) %>%
  select(barcode, clone_id.y, subpop, top_clone_highlight)

# 合并回 meta.data，注意避免 .x/.y
obj@meta.data <- obj@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(highlight_df, by = "barcode") %>%
  column_to_rownames("barcode")

emb <- obj@reductions$umap@cell.embeddings
plot_df <- tibble(
  UMAP_1   = emb[,1],
  UMAP_2   = emb[,2],
  group    = obj$condition,
  tag      = obj$top_clone_highlight,
  clone_id = obj$clone_id.y.x,
  subpop   = obj$subpop.x
) %>%
  filter(group %in% c("naive","flu_H1","flu_H5")) %>%
  mutate(
    group = factor(group, levels = c("naive","flu_H1","flu_H5")),
    tag   = factor(tag,   levels = c("Other","Naive","H1","H5"))
  )

cols <- c(
  "Other"     = "grey80",
  "Naive" = "#E69F00",
  "H1"    = "#56B4E9",
  "H5"    = "#009E73"
)

plot3d <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = filter(plot_df, tag == "Other"),
             color = cols["Other"], size = 1, alpha = 0.3) +
  geom_point(data = filter(plot_df, tag != "Other"),
             aes(color = tag, shape = subpop),
             size = 1.8, alpha = 0.9) +
  scale_color_manual(values = cols,
                     breaks = c("Naive","H1","H5"),
                     name = "Top Clones") +
  labs(title = "",
       x = "umap_1", y = "umap_2") +
  facet_wrap(~ group, ncol = 3) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial"),
    strip.text = element_text(face = "bold", size = 18, family = "Arial"),
    legend.text = element_text(size = 16, family = "Arial"),
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),
    strip.background = element_blank(),   # 去掉 facet 边框
    legend.position = "right",
    axis.title  = element_text(size = 16, family = "Arial", face = "bold"), # 坐标轴标题
    axis.text   = element_text(size = 16, family = "Arial")                 # 坐标轴刻度
  )
print(plot3d)

table_summary <- plot_df %>%
  filter(tag != "Other") %>%
  dplyr::count(group, tag, clone_id, subpop, name = "n_cells") %>%
  arrange(group, tag, desc(n_cells))

print(table_summary)
ggsave("./results/Figure3d.pdf", plot = plot3d, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure3d.png", plot = plot3d, width = 15, height = 5, dpi = 300, device = "png")


### Figure 2: ------------------------------------------------

flu_H1_mouse1.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH1_mouse1.merge.final.bcr.shm.tsv")
flu_H1_mouse2.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH1_mouse2.merge.final.bcr.shm.tsv")
flu_H5_mouse.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH5_mouse.merge.final.bcr.shm.tsv")
naive_mouse1.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/naive_mouse1.merge.final.bcr.shm.tsv")
naive_mouse2.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/naive_mouse2.merge.final.bcr.shm.tsv")

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)


####Figure 2a VH_VL usage------

library(readr)
library(dplyr)
library(stringr)
library(tidyr) 
library(ggplot2)
library(treemapify)
library(RColorBrewer)

# --- 步骤 1: 设置输出目录 ---
# !!! 重要: 请将下面的路径修改为您希望保存图片的文件夹 !!!
results_directory <- "./results/Chord_Diagrams_Final_Corrected"

# 检查输出文件夹是否存在，如果不存在则创建它
if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE)
  message("已创建结果文件夹: ", results_directory)
}

# --- 步骤 2: 将所有待处理的数据框放入一个命名的列表中 ---
# 确保这些数据框对象在您的R环境中已经存在
list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

# --- 步骤 3: 定义最终修正版的绘图函数 ---
generate_chord_diagram_final_corrected <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- 开始处理样本: ", sample_name, " ---")
  message("----------------------------------------------------")
  
  # 3.1 数据处理
  pairing_freq <- sample_data %>%
    select(heavy_v = v_call_10x, light_v = light_v_call_10x) %>%
    filter(!is.na(heavy_v) & heavy_v != "" & !is.na(light_v) & light_v != "") %>%
    mutate(
      heavy_v = str_remove(heavy_v, "\\*.*"),
      light_v = str_remove(light_v, "\\*.*")
    ) %>%
    group_by(heavy_v, light_v) %>%
    summarise(count = n(), .groups = 'drop') %>%
    dplyr::rename(from = heavy_v, to = light_v, value = count)
  
  if (nrow(pairing_freq) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的配对数据，跳过绘图。")
    return(NULL)
  }
  
  # 3.2 准备绘图数据
  data_matrix <- xtabs(value ~ from + to, data = pairing_freq)
  heavy_genes <- sort(unique(pairing_freq$from))
  light_genes <- sort(unique(pairing_freq$to))
  gene_order <- c(heavy_genes, light_genes)
  all_genes <- unique(gene_order)
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  grid_colors <- setNames(colors[1:length(all_genes)], all_genes)
  
  # 3.3 定义绘图命令函数
  create_plot_commands <- function() {
    circos.clear()
    circos.par(gap.after = c(rep(0.5, length(heavy_genes) - 1), 10, rep(0.5, length(light_genes) - 1), 10))
    
    chordDiagram(data_matrix, order = gene_order, grid.col = grid_colors,
                 annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(get.cell.meta.data("xcenter"), get.cell.meta.data("ylim")[1] + mm_y(5),
                  get.cell.meta.data("sector.index"), facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA)
    
    # --- !!! 核心修正点 !!! ---
    # 创建一个空轨道时，必须明确指定ylim
    circos.track(track.height = 0.05, bg.border = NA, ylim = c(0, 1))

    highlight.sector(sector.index = heavy_genes, track.index = 2,
                     col = "black", border = NA,
                     text = "Heavy Chain V Gene Locus", cex = 1.3,
                     text.vjust = -0.5, facing = "bending.inside")
    
    highlight.sector(sector.index = light_genes, track.index = 2,
                     col = "darkgreen", border = NA,
                     text = "Light Chain V Gene Locus", cex = 1.3,
                     text.vjust = -0.5, facing = "bending.inside")
  }
  
  # 3.4 保存文件
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Chord_Diagram.pdf"))
  message("正在生成 PDF 文件: ", pdf_file_name)
  pdf(pdf_file_name, width = 12, height = 12)
  create_plot_commands()
  dev.off()
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Chord_Diagram.png"))
  message("正在生成 PNG 文件: ", png_file_name)
  png(png_file_name, width = 1800, height = 1800, res = 150)
  create_plot_commands()
  dev.off()
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}

# --- 步骤 4: 循环处理列表中的所有样本 ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_chord_diagram_final_corrected(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}

message("\n\n所有样本处理完毕！请检查 '", results_directory, "' 文件夹。")


###top 30

# --- 步骤 0: 加载必要的R包 ---
library(readr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(circlize) # 明确加载circlize包

# --- 步骤 1: 设置输出目录 ---
results_directory <- "./results/Chord_Diagrams_Top30_Heavy"

if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE)
  message("已创建结果文件夹: ", results_directory)
}

# --- 步骤 2: 将所有待处理的数据框放入一个命名的列表中 ---

list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

# --- 步骤 3: 定义优化后的绘图函数 (筛选Top 30重链) ---
generate_chord_diagram_top30_heavy <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- 开始处理样本: ", sample_name, " ---")
  message("----------------------------------------------------")
  
  # 数据预处理：计算所有基因对的频率
  pairing_freq <- sample_data %>%
    select(heavy_v = v_call_10x, light_v = light_v_call_10x) %>%
    filter(!is.na(heavy_v) & heavy_v != "" & !is.na(light_v) & light_v != "") %>%
    mutate(
      heavy_v = str_remove(heavy_v, "\\*.*"),
      light_v = str_remove(light_v, "\\*.*")
    ) %>%
    group_by(heavy_v, light_v) %>%
    summarise(count = n(), .groups = 'drop') %>%
    dplyr::rename(from = heavy_v, to = light_v, value = count)
  
  if (nrow(pairing_freq) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的配对数据，跳过绘图。")
    return(NULL)
  }
  
  # 筛选Top 30重链基因
  heavy_gene_usage <- pairing_freq %>%
    group_by(from) %>%
    summarise(total_value = sum(value), .groups = 'drop') %>%
    arrange(desc(total_value))
  
  top_heavy_genes <- head(heavy_gene_usage$from, 20)
  
  message(sprintf("在 '%s' 中发现 %d 个独特的重链基因。将筛选出 Top %d 进行可视化。",
                  sample_name, nrow(heavy_gene_usage), length(top_heavy_genes)))
  
  pairing_freq_filtered <- pairing_freq %>%
    filter(from %in% top_heavy_genes)
  
  if (nrow(pairing_freq_filtered) == 0) {
    message("警告: 样本 '", sample_name, "' 筛选Top 30重链后无数据，跳过绘图。")
    return(NULL)
  }
  
  # 准备绘图数据
  data_matrix <- xtabs(value ~ from + to, data = pairing_freq_filtered)
  heavy_genes <- sort(unique(pairing_freq_filtered$from))
  light_genes <- sort(unique(pairing_freq_filtered$to))
  gene_order <- c(heavy_genes, light_genes)
  all_genes <- unique(gene_order)
  
  # 设置颜色
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  grid_colors <- setNames(rep_len(colors, length.out = length(all_genes)), all_genes)
  
  # 定义绘图命令函数
  create_plot_commands <- function() {
    circos.clear()
    circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2), start.degree = 90)
    
    # 绘制和弦图 (彩色弧)
    chordDiagram(data_matrix, order = gene_order, grid.col = grid_colors,
                 annotationTrack = "grid", 
                 preAllocateTracks = list(
                   list(track.height = 0.05), # 轨道1 (彩色弧)
                   list(track.height = 0.15)  # 轨道3 (基因标签)
                 ))
    
    # 步骤 A: 在最内圈 (轨道1) 绘制彩色背景弧
    circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {})
    #highlight.sector(sector.index = heavy_genes, track.index = 1, col = "#2171b5")
    #highlight.sector(sector.index = light_genes, track.index = 1, col = "#238b45")
    
    # 步骤 C: 在最外圈 (轨道3) 绘制基因标签
    circos.track(track.index = 2, bg.border = NA, panel.fun = function(x, y) {
      circos.text(get.cell.meta.data("xcenter"), 
                  get.cell.meta.data("ylim")[1] + mm_y(4),
                  get.cell.meta.data("sector.index"), 
                  facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), 
                  cex = 0.6) # 调整标签字体大小
    })
  }
  
  # 保存文件
  save_plot_with_titles <- function(file_path, device, width, height, res = NA) {
    message("正在生成文件: ", file_path)
    
    if (device == "pdf") {
      pdf(file_path, width = width, height = height)
    } else {
      png(file_path, width = width, height = height, res = res)
    }
    
    par(mar = c(1, 1, 4, 1)) # 预留顶部空间
    create_plot_commands()
    title(paste("V-Gene Pairing Frequencies for", sample_name), cex.main = 2.0)
    dev.off()
  }
  
  save_plot_with_titles(
    file_path = file.path(output_dir, paste0(sample_name, "_Top30_Heavy_Chord.pdf")),
    device = "pdf", width = 9, height = 9
  )
  save_plot_with_titles(
    file_path = file.path(output_dir, paste0(sample_name, "_Top30_Heavy_Chord.png")),
    device = "png", width = 1500, height = 1500, res = 150
  )
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}

# 循环处理列表中的所有样本
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_chord_diagram_top30_heavy(
      sample_data = current_data, 
      sample_name = name, 
      output_dir = results_directory
    )
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}


####Figure 2b VH and VL density------
# --- 步骤 0: 安装和加载所有必要的R包 ---
# 这段代码会自动检查并安装所有需要的包

packages_to_install <- c("dplyr", "readr", "tidyr", "stringr", "ggplot2", 
                         "treemapify", "RColorBrewer", "showtext")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(treemapify)
library(RColorBrewer)


# --- 步骤 2: 加载、合并并准备所有数据 ---

# 2.1 加载5个原始文件
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"))

# 2.2 创建2个合并文件
flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 2.3 将所有7个数据框放入一个命名的列表中
list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

message("所有7个数据框已成功加载并准备就绪。")


# --- 步骤 3: 设置输出目录 ---
results_directory <- "./results/Treemap"

tryCatch({
  if (!dir.exists(results_directory)) {
    dir.create(results_directory, recursive = TRUE)
    message("成功创建文件夹: ", results_directory)
  } else {
    message("文件夹已存在: ", results_directory)
  }
}, error = function(e) {
  # 如果即使有 recursive=TRUE 仍然失败，那几乎肯定是权限问题
  message("创建文件夹失败！这很可能是由于网络共享的写入权限不足。")
  message("原始错误信息: ", e$message)
})

# --- 步骤 4: 定义绘图函数 ---
# (与之前相同，但现在它能找到Arial字体了)
generate_publication_treemap <- function(sample_data, sample_name, output_dir) {
  message("\n--- 开始处理样本: ", sample_name, " ---")
  
  processed_H_fre <- sample_data %>%
    tidyr::unnest(v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = v_call_10x) %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    count(heavy_v, name = "total_value") %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH")) %>%
    arrange(desc(total_value))
  
  if (nrow(processed_H_fre) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的重链数据，跳过绘图。")
    return(NULL)
  }
  
  unique_labels <- processed_H_fre$label
  num_colors <- length(unique_labels)
  high_sat_pals <- c("Dark2", "Set1", "Accent", "Set2")
  all_colors <- unlist(sapply(high_sat_pals, function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)))
  unique_colors <- unique(all_colors)
  if (length(unique_colors) < num_colors) {
    color_vector <- rep(unique_colors, length.out = num_colors)
  } else {
    color_vector <- unique_colors[1:num_colors]
  }
  color_palette <- setNames(sample(color_vector, num_colors), unique_labels)
  
  treemap_plot <- ggplot(processed_H_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6, family = "Arial"
    ) +
    scale_fill_manual(values = color_palette) +
    labs(title = paste("Heavy Chain V Gene Usage:", sample_name)) +
    theme_minimal(base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.pdf"))
  message("正在生成 PDF 文件: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 8, device = cairo_pdf) # 使用cairo_pdf确保字体嵌入
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.png"))
  message("正在生成 PNG 文件: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 8, dpi = 300, device = "png")
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}

# --- 步骤 5: 循环处理 ---
# (与之前相同)
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_publication_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}

####VL density
# --- 步骤 3: 设置输出目录 ---
results_directory <- "./results/Treemap_Light"

tryCatch({
  if (!dir.exists(results_directory)) {
    dir.create(results_directory, recursive = TRUE)
    message("成功创建文件夹: ", results_directory)
  } else {
    message("文件夹已存在: ", results_directory)
  }
}, error = function(e) {
  # 如果即使有 recursive=TRUE 仍然失败，那几乎肯定是权限问题
  message("创建文件夹失败！这很可能是由于网络共享的写入权限不足。")
  message("原始错误信息: ", e$message)
})

# --- 步骤 4: 定义绘图函数 ---
# (与之前相同，但现在它能找到Arial字体了)
generate_publication_treemap <- function(sample_data, sample_name, output_dir) {
  message("\n--- 开始处理样本: ", sample_name, " ---")
  
  processed_H_fre <- sample_data %>%
    tidyr::unnest(light_v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = light_v_call_10x) %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    count(heavy_v, name = "total_value") %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH")) %>%
    arrange(desc(total_value))
  
  if (nrow(processed_H_fre) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的重链数据，跳过绘图。")
    return(NULL)
  }
  
  unique_labels <- processed_H_fre$label
  num_colors <- length(unique_labels)
  high_sat_pals <- c("Dark2", "Set1", "Accent", "Set2")
  all_colors <- unlist(sapply(high_sat_pals, function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)))
  unique_colors <- unique(all_colors)
  if (length(unique_colors) < num_colors) {
    color_vector <- rep(unique_colors, length.out = num_colors)
  } else {
    color_vector <- unique_colors[1:num_colors]
  }
  color_palette <- setNames(sample(color_vector, num_colors), unique_labels)
  
  treemap_plot <- ggplot(processed_H_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6, family = "Arial"
    ) +
    scale_fill_manual(values = color_palette) +
    labs(title = paste("Heavy Chain V Gene Usage:", sample_name)) +
    theme_minimal(base_family = "Arial") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.pdf"))
  message("正在生成 PDF 文件: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 8, device = cairo_pdf) # 使用cairo_pdf确保字体嵌入
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.png"))
  message("正在生成 PNG 文件: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 8, dpi = 300, device = "png")
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}

# --- 步骤 5: 循环处理 ---
# (与之前相同)
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_publication_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}


####top30
# --- Heavy top 30 ---
# --- 步骤 0: 检查并加载必要的R包 ---
packages_to_install <- c("dplyr", "readr", "stringr", "tidyr", 
                         "ggplot2", "treemapify", "RColorBrewer")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(treemapify)
library(RColorBrewer)


# --- 步骤 1: 加载、合并并准备所有数据 ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- readr::read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"))
flu_H1_mouse2.bcr <- readr::read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"))
flu_H5_mouse.bcr  <- readr::read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"))
naive_mouse1.bcr  <- readr::read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"))
naive_mouse2.bcr  <- readr::read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"))

flu_H1_mouse.bcr <- dplyr::bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr  <- dplyr::bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

list_of_samples <- list(
  flu_H1_mouse1         = flu_H1_mouse1.bcr,
  flu_H1_mouse2         = flu_H1_mouse2.bcr,
  flu_H5_mouse          = flu_H5_mouse.bcr,
  naive_mouse1          = naive_mouse1.bcr,
  naive_mouse2          = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined  = naive_mouse.bcr
)
message("所有7个数据框已成功加载并准备就绪。")


# --- 步骤 2: 设置输出目录 ---
# 建议使用一个更具体的文件夹名，以反映分析内容
results_directory <- "./results/Treemap_Heavy_Top30"
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)
message("输出目录已准备就绪: ", results_directory)


# --- 步骤 3: 定义优化后的绘图函数 (仅绘制Top 30重链) ---
generate_heavy_chain_top30_treemap <- function(sample_data, sample_name, output_dir) {
  
  message("\n--- 开始处理样本: ", sample_name, " ---")
  
  # 3.1 数据处理：计算频率并筛选Top 30
  processed_H_fre <- sample_data %>%
    tidyr::unnest(v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = v_call_10x) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    count(heavy_v, name = "total_value") %>%
    arrange(desc(total_value)) %>%
    # --- !!! 核心修改点: 仅选取频率最高的30个基因 !!! ---
    slice_head(n = 20) %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH"))
  
  if (nrow(processed_H_fre) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的重链数据，跳过绘图。")
    return(NULL)
  }
  
  # 3.2 颜色管理
  num_colors <- nrow(processed_H_fre)
  color_vector <- RColorBrewer::brewer.pal.info %>%
    filter(category == "qual") %>%
    row.names() %>%
    lapply(function(pal) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal)) %>%
    unlist() %>%
    unique()
  color_palette <- setNames(rep_len(color_vector, length.out = num_colors), processed_H_fre$label)
  
  # 3.3 创建绘图对象
  treemap_plot <- ggplot(processed_H_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6
    ) +
    scale_fill_manual(values = color_palette) +
    # 更新标题以反映Top 30的筛选
    labs(title = paste("Top 30 Heavy Chain V Gene Usage:", sample_name)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 3.4 保存文件
  # 使用更具描述性的文件名
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Heavy_Top30_Treemap.pdf"))
  message("正在生成 PDF 文件: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 4, device = "pdf")
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Heavy_Top30_Treemap.png"))
  message("正在生成 PNG 文件: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 4, dpi = 300, device = "png")
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}


# --- 步骤 4: 循环处理 ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    # 调用优化后的新函数
    generate_heavy_chain_top30_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}

# --- Light top 30 ---
# --- 步骤 2: 设置轻链分析的输出目录 ---
results_directory <- "./results/Treemap_Light_Top30"
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)
message("轻链分析的输出目录已准备就绪: ", results_directory)


# --- 步骤 3: 定义优化后的绘图函数 (仅绘制Top 30轻链) ---
generate_light_chain_top30_treemap <- function(sample_data, sample_name, output_dir) {
  
  message("\n--- 开始处理样本: ", sample_name, " (Light Chain) ---")
  
  # 3.1 数据处理：计算轻链频率并筛选Top 30
  processed_L_fre <- sample_data %>%
    # --- !!! 核心修改点 1: 从 light_v_call_10x 列提取数据 !!! ---
    tidyr::unnest(light_v_call_10x, keep_empty = TRUE) %>%
    transmute(light_v = light_v_call_10x) %>%
    filter(!is.na(light_v) & light_v != "") %>%
    mutate(light_v = str_remove(light_v, "\\*.*")) %>%
    count(light_v, name = "total_value") %>%
    arrange(desc(total_value)) %>%
    slice_head(n = 20) %>%
    # --- !!! 核心修改点 2: 美化轻链基因标签 (IGKV -> VK, IGLV -> VL) !!! ---
    mutate(label = case_when(
      grepl("IGKV", light_v) ~ str_replace(light_v, "IGKV", "VK"),
      grepl("IGLV", light_v) ~ str_replace(light_v, "IGLV", "VL"),
      TRUE ~ light_v # 保留原始名称以防意外情况
    ))
  
  if (nrow(processed_L_fre) == 0) {
    message("警告: 样本 '", sample_name, "' 没有有效的轻链数据，跳过绘图。")
    return(NULL)
  }
  
  # 3.2 颜色管理
  num_colors <- nrow(processed_L_fre)
  color_vector <- RColorBrewer::brewer.pal.info %>%
    filter(category == "qual") %>%
    row.names() %>%
    lapply(function(pal) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal)) %>%
    unlist() %>%
    unique()
  color_palette <- setNames(rep_len(color_vector, length.out = num_colors), processed_L_fre$label)
  
  # 3.3 创建绘图对象
  treemap_plot <- ggplot(processed_L_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6
    ) +
    scale_fill_manual(values = color_palette) +
    # --- !!! 核心修改点 3: 更新图表标题 !!! ---
    labs(title = paste("Top 30 Light Chain V Gene Usage:", sample_name)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 3.4 保存文件
  # 使用更具描述性的文件名
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Light_Top30_Treemap.pdf"))
  message("正在生成 PDF 文件: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 5, device = "pdf")
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Light_Top30_Treemap.png"))
  message("正在生成 PNG 文件: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 5, dpi = 300, device = "png")
  
  message("--- 样本 '", sample_name, "' 处理完成 ---")
}


# --- 步骤 4: 循环处理 ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    # 调用新的、专门用于轻链分析的函数
    generate_light_chain_top30_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("处理样本 '", name, "' 时发生错误: ", e$message)
  })
}

####Figure 2c V-D-J------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(readr)

# 设置工作目录（请修改为您的数据所在目录）
base_path <- "./cell_reports/VDJ/annotation_results/"
setwd(base_path)

# 创建输出目录
output_dir <- "./results/VDJ_Sankey_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 读取所有样本数据
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# 创建合并样本
flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 创建样本列表
sample_list <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

# 数据预处理函数
prepare_sankey_data <- function(data) {
  # 选择需要的列
  sankey_data <- data %>%
    select(v_call_10x, j_call_10x, light_v_call_10x, light_j_call_10x) %>%
    # 清理基因名称，去除等位基因信息
    mutate(
      VH = gsub("\\*.*", "", v_call_10x),
      JH = gsub("\\*.*", "", j_call_10x),
      VL = gsub("\\*.*", "", light_v_call_10x),
      JL = gsub("\\*.*", "", light_j_call_10x)
    ) %>%
    # 过滤掉缺失值
    filter(!is.na(VH) & !is.na(JH) & !is.na(VL) & !is.na(JL)) %>%
    # 创建唯一克隆ID
    mutate(clone_id = paste0("clone_", row_number()))
  
  return(sankey_data)
}

# 桑基图绘制函数
generate_sankey_plot <- function(sample_data, sample_name, output_dir) {
  # 准备桑基图数据
  sankey_data <- prepare_sankey_data(sample_data)
  
  # 如果数据为空，跳过
  if (nrow(sankey_data) == 0) {
    message("样本 '", sample_name, "' 没有有效的V-J配对数据，跳过绘图。")
    return(NULL)
  }
  
  # 创建长格式数据
  sankey_long <- sankey_data %>%
    select(clone_id, VH, JH, VL, JL) %>%
    pivot_longer(
      cols = c(VH, JH, VL, JL),
      names_to = "axis",
      values_to = "gene"
    ) %>%
    mutate(
      gene_type = ifelse(axis %in% c("VH", "VL"), "V", "J"),
      chain = ifelse(axis %in% c("VH", "JH"), "Heavy", "Light")
    )
  
  # 创建V基因和J基因的调色板
  v_genes <- unique(sankey_long$gene[sankey_long$gene_type == "V"])
  j_genes <- unique(sankey_long$gene[sankey_long$gene_type == "J"])
  
  v_palette <- colorRampPalette(brewer.pal(8, "Set1"))(length(v_genes))
  j_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(j_genes))
  
  # 合并调色板
  gene_colors <- setNames(c(v_palette, j_palette), c(v_genes, j_genes))
  
  # 绘制桑基图
  sankey_plot <- ggplot(sankey_long,
                        aes(x = axis, stratum = gene, alluvium = clone_id,
                            y = 1, label = gene, fill = gene)) +
    # 绘制流和节点
    geom_flow(alpha = 0.7, color = "white", width = 1/3) +
    geom_stratum(width = 1/3, color = "black", linewidth = 0.2) +
    
    # 设置颜色
    scale_fill_manual(values = gene_colors) +
    
    # 添加文字标签
    geom_text(stat = "stratum", size = 2.5) +
    
    # 设置坐标轴
    scale_x_discrete(
      limits = c("VH", "JH", "VL", "JL"),
      labels = c("VH" = "Heavy V", "JH" = "Heavy J", 
                 "VL" = "Light V", "JL" = "Light J"),
      expand = c(0.15, 0.05)
    ) +
    
    # 设置标题和主题
    labs(
      title = paste("BCR V-J Gene Pairing:", sample_name),
      x = "Gene Segment",
      y = "Frequency",
      fill = "Gene"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",  # 删除图例
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
      axis.text.x = element_text(hjust = 1, size = 15),
      axis.text.y = element_blank(),  # 隐藏y轴文本
      axis.title.y = element_blank(),  # 隐藏y轴标题
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)  # 调整边距
    )
  
  # 保存图形
  pdf_file <- file.path(output_dir, paste0(sample_name, "_sankey_plot.pdf"))
  png_file <- file.path(output_dir, paste0(sample_name, "_sankey_plot.png"))
  
  ggsave(pdf_file, plot = sankey_plot, width = 7, height = 14, device = "pdf")
  ggsave(png_file, plot = sankey_plot, width = 7, height = 14, dpi = 300, device = "png", bg = "white")
  
  message("已为样本 '", sample_name, "' 生成桑基图")
  
  return(sankey_plot)
}

# 批量处理所有样本
for (sample_name in names(sample_list)) {
  message("正在处理样本: ", sample_name)
  tryCatch({
    generate_sankey_plot(sample_list[[sample_name]], sample_name, output_dir)
  }, error = function(e) {
    message("处理样本 '", sample_name, "' 时出错: ", e$message)
  })
}


# --- 步骤 0: 设置环境并加载必要的R包 ---
# --- 步骤 0: 设置环境并加载必要的R包 ---
# 确保已安装这些包: install.packages(c("dplyr", "tidyr", "ggplot2", "ggalluvial", "RColorBrewer", "readr"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(readr)


# --- 步骤 1: 加载、合并并准备所有数据 ---

# 1.1 定义路径和创建输出目录
base_path <- "./cell_reports/VDJ/annotation_results/"
output_dir <- "./results/VDJ_Sankey_Top30_VH_Plots_Large" # 使用新目录以防覆盖
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1.2 高效读取所有样本数据
message("正在加载数据...")
flu_H1_mouse1.bcr <- readr::read_tsv(file.path(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- readr::read_tsv(file.path(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr  <- readr::read_tsv(file.path(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr  <- readr::read_tsv(file.path(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr  <- readr::read_tsv(file.path(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# 1.3 创建合并样本
flu_H1_mouse.bcr <- dplyr::bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr  <- dplyr::bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 1.4 将所有数据框放入一个命名的列表中
sample_list <- list(
  flu_H1_mouse1         = flu_H1_mouse1.bcr,
  flu_H1_mouse2         = flu_H1_mouse2.bcr,
  flu_H5_mouse          = flu_H5_mouse.bcr,
  naive_mouse1          = naive_mouse1.bcr,
  naive_mouse2          = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined  = naive_mouse.bcr
)
message("所有数据已加载并准备就绪。")


prepare_sankey_data <- function(data, top_n = 20) {
  processed_data <- data %>%
    select(v_call_10x, j_call_10x, light_v_call_10x, light_j_call_10x) %>%
    mutate(
      VH = gsub("\\*.*", "", v_call_10x),
      JH = gsub("\\*.*", "", j_call_10x),
      VL = gsub("\\*.*", "", light_v_call_10x),
      JL = gsub("\\*.*", "", light_j_call_10x)
    ) %>%
    filter(!is.na(VH) & VH != "" & !is.na(JH) & JH != "" & 
             !is.na(VL) & VL != "" & !is.na(JL) & JL != "")
  
  if (nrow(processed_data) == 0) return(processed_data)
  
  vh_usage <- processed_data %>% count(VH, sort = TRUE)
  top_vh_genes <- head(vh_usage$VH, top_n)
  
  message(sprintf("共发现 %d 个独特的VH基因，将筛选出 Top %d 用于可视化。", 
                  nrow(vh_usage), length(top_vh_genes)))
  
  sankey_data <- processed_data %>%
    filter(VH %in% top_vh_genes) %>%
    mutate(clone_id = paste0("clone_", row_number()))
  
  return(sankey_data)
}



generate_sankey_plot <- function(
    sample_data, sample_name, output_dir,
    # --- !!! 修改点 1: 在这里调整图片尺寸 !!! ---
    output_width = 10,        # 输出图片的宽度 (单位: 英寸)
    base_height = 1.5,          # 图片的基础高度 (单位: 英寸)
    height_per_gene = 0.15,    # 每增加一个基因，图片高度额外增加的量
    # --- !!! 修改点 2: 在这里调整字体大小 !!! ---
    title_size = 22,          # 图表标题的字体大小
    axis_text_size = 18,      # 坐标轴标签 (Heavy V, etc.) 的字体大小
    gene_label_size = 4     # 基因名称标签的字体大小
) {
  
  if (nrow(sample_data) == 0) {
    message("样本 '", sample_name, "' 在筛选后没有有效的V-J配对数据，跳过绘图。")
    return(NULL)
  }
  
  # 数据转换为长格式
  sankey_long <- sample_data %>%
    select(clone_id, VH, JH, VL, JL) %>%
    pivot_longer(cols = -clone_id, names_to = "axis", values_to = "gene") %>%
    mutate(axis = factor(axis, levels = c("VH", "JH", "VL", "JL")))
  
  # 调色板管理
  all_genes <- unique(sankey_long$gene)
  num_colors <- length(all_genes)
  color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
  gene_colors <- setNames(color_palette, all_genes)
  
  # 绘制桑基图
  sankey_plot <- ggplot(data = sankey_long,
                        aes(x = axis, stratum = gene, alluvium = clone_id,
                            y = 1, fill = gene, label = gene)) +
    geom_flow(alpha = 0.6, color = "white", width = 0.4) +
    geom_stratum(width = 0.4, color = "grey30") +
    # 使用参数化的字体大小
    geom_text(stat = "stratum", size = gene_label_size, color = "black") + # <-- 使用 gene_label_size
    scale_fill_manual(values = gene_colors) +
    scale_x_discrete(
      limits = c("VH", "JH", "VL", "JL"),
      labels = c("Heavy V", "Heavy J", "Light V", "Light J")
    ) +
    labs(
      title = paste("BCR V-J Pairing:", sample_name),
      x = NULL, y = ""
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      legend.position = "none",
      # 使用参数化的字体大小
      plot.title = element_text(hjust = 0.2, size = title_size, face = "bold"), 
      axis.text.x = element_text(
        size = axis_text_size, 
        face = "bold",
        # 使用 margin() 设置边距 (上, 右, 下, 左)
        # 将上边距 (t) 设置为负值，即可向上移动文本，减小间距
        margin = margin(t = -30, unit = "pt") # -10是一个不错的起始值，您可以根据需要调整
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  
  # 动态计算高度并保存图形
  num_strata <- length(unique(sankey_long$gene))
  # 使用参数化的尺寸
  plot_height <- max(base_height, num_strata * height_per_gene) # <-- 使用 base_height 和 height_per_gene
  
  pdf_file <- file.path(output_dir, paste0(sample_name, "_Top30VH_Sankey.pdf"))
  png_file <- file.path(output_dir, paste0(sample_name, "_Top30VH_Sankey.png"))
  
  # 使用参数化的宽度
  ggsave(pdf_file, plot = sankey_plot, width = output_width, height = plot_height, device = "pdf", limitsize = FALSE) # <-- 使用 output_width
  ggsave(png_file, plot = sankey_plot, width = output_width, height = plot_height, dpi = 300, device = "png", bg = "white", limitsize = FALSE) # <-- 使用 output_width
  
  message("已为样本 '", sample_name, "' 生成桑基图。")
}


for (name in names(sample_list)) {
  message("\n====================================================")
  message("--- 开始处理样本: ", name, " ---")
  message("====================================================")
  
  tryCatch({
    data_for_plotting <- prepare_sankey_data(sample_list[[name]], top_n = 20)
    
    # 调用绘图函数 (所有尺寸和字体调整都在函数内部完成)
    generate_sankey_plot(data_for_plotting, name, output_dir)
    
  }, error = function(e) {
    message("!!!!!! 处理样本 '", name, "' 时发生严重错误: ", e$message, " !!!!!!")
  })
}


### Figure 4: ------------------------------------------------
#####Figure4a: Isotype ------------------------------------------------

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "RColorBrewer", "scales")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(dplyr); library(readr); library(stringr); library(ggplot2); library(RColorBrewer); library(scales)

base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

list_of_samples <- list(
  "Flu H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

isotype_proportions <- all_data %>%
  dplyr::filter(locus == "IGH", !is.na(.data$c_call), .data$c_call != "") %>%
  dplyr::rename(isotype = c_call) %>%
  # 按组别和同种型计数
  dplyr::count(group, isotype, name = "count") %>%
  # 按组别计算比例
  dplyr::group_by(group) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(proportion))

# 定义同种型的顺序（在图例中显示的顺序）
all_isotypes_found <- unique(isotype_proportions$isotype)
isotype_levels <- sort(all_isotypes_found) 
isotype_proportions$isotype <- factor(isotype_proportions$isotype, levels = rev(isotype_levels))

# 定义组别的顺序
group_levels <- c("Naive", "Flu H1", "Flu H5")
isotype_proportions$group <- factor(isotype_proportions$group, levels = group_levels)

print(isotype_proportions) 

unique_isotypes <- levels(isotype_proportions$isotype)
color_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique_isotypes)),
  unique_isotypes
)

plot4a <- ggplot(isotype_proportions, aes(x = group, y = proportion, fill = isotype)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = color_palette, name = "Isotype") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell Fraction (%)") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(hjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(plot4a)
ggsave("./results/Figure4a.pdf", plot = plot4a, width = 6, height = 5, device = "pdf")
ggsave("./results/Figure4a.png", plot = plot4a, width = 6, height = 5, dpi = 300, device = "png")

#####FigureS13----

# 1.1  创建一个包含5个独立样本数据框的列表
list_of_samples <- list(
  "Flu_H1_Mouse1" = flu_H1_mouse1.bcr,
  "Flu_H1_Mouse2" = flu_H1_mouse2.bcr,
  "Flu_H5_Mouse" = flu_H5_mouse.bcr,
  "Naive_Mouse1" = naive_mouse1.bcr,
  "Naive_Mouse2" = naive_mouse2.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

isotype_proportions <- all_data %>%
  dplyr::filter(locus == "IGH", !is.na(.data$c_call), .data$c_call != "") %>%
  dplyr::rename(isotype = c_call) %>%
  # 按组别和同种型计数
  dplyr::count(group, isotype, name = "count") %>%
  # 按组别计算比例
  dplyr::group_by(group) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(proportion))

# 定义同种型的顺序（在图例中显示的顺序）
all_isotypes_found <- unique(isotype_proportions$isotype)
isotype_levels <- sort(all_isotypes_found) 
isotype_proportions$isotype <- factor(isotype_proportions$isotype, levels = rev(isotype_levels))

group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")
isotype_proportions$group <- factor(isotype_proportions$group, levels = group_levels)
print(isotype_proportions) 

# --- 步骤 3: 绘制出版级100%堆叠柱状图 ---
unique_isotypes <- levels(isotype_proportions$isotype)
color_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique_isotypes)),
  unique_isotypes
)

plotS13 <- ggplot(isotype_proportions, aes(x = group, y = proportion, fill = isotype)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = color_palette, name = "Isotype") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell Fraction (%)") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(plotS13)
ggsave("./results/FigureS13.pdf", plot = plotS13, width = 7, height = 7, device = "pdf")
ggsave("./results/FigureS13.png", plot = plotS13, width = 7, height = 7, dpi = 300, device = "png")

#####Figure4b Shannon's Entropy -------------------

library(vegan)

# 读入5个样本
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr  <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"),  col_types = cols(.default = "c"))
naive_mouse1.bcr  <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr  <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# 合并得到2个新样本
naive_merge <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)
flu_merge   <- bind_rows(flu_H1_mouse2.bcr, flu_H5_mouse.bcr)


sample_list <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse  = flu_H5_mouse.bcr,
  naive_mouse1  = naive_mouse1.bcr,
  naive_mouse2  = naive_mouse2.bcr
)

calc_diversity <- function(df) {
  df %>%
    mutate(umi_count = as.numeric(umi_count)) %>%
    group_by(clone_id) %>%
    summarise(count = sum(umi_count, na.rm = TRUE), .groups = "drop") %>%
    { 
      x <- .$count
      tibble(
        shannon = diversity(x, index = "shannon"),
        simpson = diversity(x, index = "simpson")
      )
    }
}

diversity_results <- map_df(sample_list, calc_diversity, .id = "sample") %>%
  mutate(group = case_when(
    grepl("naive", sample) ~ "Naive",
    grepl("flu_H1", sample) ~ "Flu_H1",
    grepl("flu_H5", sample) ~ "Flu_H5"
  ))


cell_counts <- map_df(sample_list, ~tibble(n_cells = nrow(.x)), .id = "sample") %>%
  mutate(group = case_when(
    grepl("naive", sample) ~ "Naive",
    grepl("flu_H1", sample) ~ "Flu_H1",
    grepl("flu_H5", sample) ~ "Flu_H5"
  )) %>%
  group_by(group) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")


color_palette <- c(
  "Naive" = "#0072B2", 
  "Flu_H1" = "#E69F00", 
  "Flu_H5" = "#D55E00"
)
diversity_results <- diversity_results %>%
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))
cell_counts <- cell_counts %>%
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))


plot4b <- ggplot(diversity_results, aes(x = group, y = shannon, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, shape = 21, color = "black", fill = "white") +
  scale_fill_manual(values = color_palette) +
  labs(title = "BCR Diversity",
       x = NULL, y = "Shannon entropy")  +
  geom_text(data = cell_counts,
            aes(x = group, y = min(diversity_results$shannon) - 0.5, 
                label = paste0("n=", n_cells)),
            inherit.aes = FALSE,
            size = 4, family = "Arial", fontface = "bold")


plot4b<- plot4b + theme_minimal(base_size = 16, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # 背景透明
    plot.background = element_rect(fill = "transparent", colour = NA),  # 图整体透明
    panel.grid = element_blank(),       # 去掉网格线
    axis.line = element_line(color = "black"), # 只保留横纵轴
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    axis.text = element_text(size = 16, face = "bold",family = "Arial"),
    axis.title = element_text(size = 18, face = "bold", family = "Arial"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial")
  )
print(plot4b)
ggsave("./results/Figure4b.pdf", plot = plot4b, width = 6, height = 5)
ggsave("./results/Figure4b.png", plot = plot4b, width = 6, height = 5)

#####Figure4c CDR3 clones share in H1 and H5 -------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggalluvial)   # 用于飘带连接


CDR3_AA_COLUMN <- "junction_10x_aa"          # CDR3 AA 列
SUBSET_COL     <- "B_cell_subpopulations"    # 亚群列（这里不再用于汇总，只做过滤用）
COND_LEVELS    <- c("flu_H1","flu_H5")       # 只保留 H1 与 H5（不含 naive）
TOP_N          <- 10                          # 画前多少条共享 CDR3；不设 Other

meta_df1 <- obj_bcr@meta.data %>%
  tibble::rownames_to_column(var = "barcode") %>%
  tibble::as_tibble()

dat <- meta_df1 %>%
  dplyr::filter(!is.na(.data[[CDR3_AA_COLUMN]]), .data[[CDR3_AA_COLUMN]] != "") %>%
  dplyr::transmute(
    condition = factor(.data$condition),                 # 原始组别
    subset    = .data[[SUBSET_COL]],                     # 这里先保留（不参与汇总）
    cdr3      = toupper(.data[[CDR3_AA_COLUMN]])         # 统一成大写
  ) %>%
  dplyr::filter(condition %in% COND_LEVELS) %>%          # 只保留 H1 与 H5
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))


shared_h1h5 <- dat %>%
  dplyr::distinct(condition, cdr3) %>%
  dplyr::add_count(cdr3, name = "n_groups") %>%
  dplyr::filter(n_groups == length(COND_LEVELS)) %>%     # 两组都出现
  dplyr::pull(cdr3) %>%
  unique()


df_flow <- dat %>%
  dplyr::filter(cdr3 %in% shared_h1h5) %>%
  dplyr::count(condition, cdr3, name = "count") %>%      # 每组该 CDR3 的细胞数
  dplyr::group_by(condition) %>%
  dplyr::mutate(percentage = 100 * count / sum(count)) %>%
  dplyr::ungroup()

#
top_shared <- df_flow %>%
  dplyr::group_by(cdr3) %>%
  dplyr::summarise(total = sum(count), .groups = "drop") %>%
  dplyr::slice_max(total, n = TOP_N, with_ties = FALSE) %>%
  dplyr::pull(cdr3)

df_flow_top <- df_flow %>%
  dplyr::filter(cdr3 %in% top_shared) %>%
  # 保证每条 CDR3 在两个组都有一行（若极少数缺失，补 0，通常不会发生，因为已筛共享）
  tidyr::complete(condition, cdr3, fill = list(count = 0, percentage = 0)) %>%
  dplyr::mutate(
    condition = factor(condition, levels = COND_LEVELS),
    cdr3 = forcats::fct_reorder(cdr3, percentage, .fun = max, .desc = TRUE) # 图例按占比排序
  )
pal_example <- c(
  "#9ACD66", "#8FCAD0", "#F39A89", "#A7B7D9", "#CFE6A0",
  "#F4D57E", "#D9B07E", "#AFBE77", "#D6C7E8", "#B390C6"
)

# 按当前因子水平数量截取（不依赖 CDR3 文本；TOP_N≦10时正好用前 n 色）
n_lv <- nlevels(df_flow_top$cdr3)
pal_use <- if (n_lv <= length(pal_example)) {
  pal_example[seq_len(n_lv)]
} else {
  # 极少数需要 >10 色时，用原风格插值扩展
  grDevices::colorRampPalette(pal_example)(n_lv)
}

plot4c <- ggplot(
  df_flow_top,
  aes(x = condition, y = percentage, alluvium = cdr3, stratum = cdr3, fill = cdr3)
) +
  geom_alluvium(alpha = 0.65, color = NA, width = 0.45) +
  geom_stratum(width = 0.45, color = "grey30", size = 0.2) +
  scale_x_discrete(labels = c("flu_H1" = "H1", "flu_H5" = "H5")) +
  scale_fill_manual(values = pal_use) +
  labs(
    x = NULL,
    y = "Percentage of Total CDR3 (%)",
    fill = "CDR3 Sequence",
    title = "Shared CDR3 between H1 and H5"
  ) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    # 横纵坐标文字加粗
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    # 坐标刻度文字也加粗（可选）
    axis.text = element_text(face = "bold"),
    # 图标题加粗并居中
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(plot4c)

ggsave("./results/Figure4c.pdf", plot = plot4c, width = 8, height = 5)
ggsave("./results/Figure4c.png", plot = plot4c, width = 8, height = 5)

#####Figure4d HCDR3 Length ------------------------------------------------
# --- 步骤 0: 安装和加载所有必要的R包 ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(purrr)

base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)
# 创建一个包含文件名和对应组别名的列表
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

# 提取HCDR3长度，并进行必要的清理
hcdr3_lengths <- all_data %>%
  # 只保留重链数据 (locus == "IGH")
  filter(locus == "IGH") %>%
  # 选择我们需要的列
  select(group, cdr3_aa) %>%
  # 移除任何无效的行
  filter(!is.na(cdr3_aa) & cdr3_aa != "") %>%
  # 计算HCDR3的氨基酸长度
  mutate(hcdr3_length = str_length(cdr3_aa)) %>%
  # 将组别转换为因子，并按指定顺序排列
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

# 计算每个组的样本量，用于在图上标注
sample_sizes <- hcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3长度计算完成。")
print(head(hcdr3_lengths))
print(sample_sizes)

stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = hcdr3_lengths,
  method = "wilcox.test"
)
message("统计检验结果:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\n显著的组间差异 (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("在 p.adj < 0.05 的水平上，未发现显著的组间差异。")
}

# 定义组间的比较，用于显著性检验
my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

# 定义颜色方案
color_palette <- c("Naive" = "#0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") # 蓝、橙、红

plot4d <- ggplot(hcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  # 绘制小提琴图层
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  # 在小提琴内部添加箱线图
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  # 设置颜色方案
  scale_fill_manual(values = color_palette) +
  
  # 添加显著性检验和P值
  # method = "wilcox.test" (非参数检验) 或 "t.test" (参数检验)
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", # 显示 *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  # 添加样本量标签
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 5, fontface = "bold", color = "black") +
  
  # 设置坐标轴标签和范围
  labs(
    title = "HCDR3 Length Distribution",
    y = "HCDR3 Length (aa)",
    x = "" # 隐藏X轴标题
  ) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 30, 10)) +
  
  # 设置一个简洁、专业的主题
  theme_classic() +
  theme(
    legend.position = "none", # 隐藏图例
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot4d)
ggsave("./results/Figure4d.pdf", plot = plot4d, width = 6, height = 5)
ggsave("./results/Figure4d.png", plot = plot4d, width = 6, height = 5,dpi = 300)

#####FigureS14: 5 samples for HCDR3 Length-------
# --- 步骤 1: 定义文件和组别，并读取所有5个样本的数据 ---

# 创建一个包含文件名和对应组别名的列表
files_to_process <- list(
  "Flu_H1_Mouse1" = "fluH1_mouse1.merge.final.bcr.shm.tsv",
  "Flu_H1_Mouse2" = "fluH1_mouse2.merge.final.bcr.shm.tsv",
  "Flu_H5_Mouse" = "fluH5_mouse.merge.final.bcr.shm.tsv",
  "Naive_Mouse1" = "naive_mouse1.merge.final.bcr.shm.tsv",
  "Naive_Mouse2" = "naive_mouse2.merge.final.bcr.shm.tsv"
)

# 使用 purrr::map_dfr 自动读取所有文件并合并成一个数据框
all_data <- map_dfr(files_to_process, ~ read_tsv(.x, col_types = cols(.default = "c")), .id = "group")


# --- 步骤 2: 数据处理 - 计算HCDR3长度 ---

# 定义组别的顺序，以便在X轴上按逻辑排列
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

hcdr3_lengths <- all_data %>%
  filter(locus == "IGH") %>%
  select(group, cdr3_aa) %>%
  filter(!is.na(cdr3_aa) & cdr3_aa != "") %>%
  mutate(hcdr3_length = str_length(cdr3_aa)) %>%
  mutate(group = factor(group, levels = group_levels))

# 计算每个组的样本量
sample_sizes <- hcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3长度计算完成。")

# --- 步骤 2.5: 执行并解读显著性检验 ---
message("\n----------------------------------------------------")
message("--- 执行5个样本间的HCDR3长度显著性检验 ---")
message("----------------------------------------------------")
stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = hcdr3_lengths,
  method = "wilcox.test"
)
message("统计检验结果:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\n显著的组间差异 (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("在 p.adj < 0.05 的水平上，未发现显著的组间差异。")
}
message("----------------------------------------------------\n")


# --- 步骤 3: 使用 ggplot2 和 ggpubr 绘制小提琴图 ---

message("开始绘制小提琴图...")

# 定义组间的比较，用于在图上显示

my_comparisons <- list( 
  c("Naive_Mouse1", "Naive_Mouse2"),
  c("Naive_Mouse2", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Flu_H5_Mouse")
)

# 定义颜色方案
color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plotS14 <- ggplot(hcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  
  # 添加显著性检验和P值
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, size = 5) +
  
  # 添加样本量标签
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 6, fontface = "bold", color = "black") +
  
  # 设置坐标轴标签和范围
  labs(
    title = "HCDR3 Length Distribution (5 Samples)",
    y = "HCDR3 Length (aa)",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 30, 10)) +
  
  # 设置一个简洁、专业的主题
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    # 调整X轴标签，使其不重叠
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plotS14)
ggsave("./results/FigureS14.pdf", plot = plotS14, width = 8, height = 7)
ggsave("./results/FigureS14.png", plot = plotS14, width = 8, height = 7,dpi = 300)


#####Figure4e LCDR3 Length ------------------------------------------------

# --- 步骤 0: 安装和加载所有必要的R包 ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(purrr)

# --- 步骤 1: 定义文件和组别，并读取所有数据 ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)
# 创建一个包含文件名和对应组别名的列表
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")


# --- 步骤 2: 数据处理 - 计算HCDR3长度 ---

# 提取LCDR3长度，并进行必要的清理
lcdr3_lengths <- all_data %>%
  # 只保留重链数据 (locus == "IGH")
  #filter(locus == "IGH") %>%
  # 选择我们需要的列
  select(group, light_cdr3_aa) %>%
  # 移除任何无效的行
  filter(!is.na(light_cdr3_aa) & light_cdr3_aa != "") %>%
  # 计算HCDR3的氨基酸长度
  mutate(hcdr3_length = str_length(light_cdr3_aa)) %>%
  # 将组别转换为因子，并按指定顺序排列
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

# 计算每个组的样本量，用于在图上标注
sample_sizes <- lcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3长度计算完成。")
print(head(lcdr3_lengths))
print(sample_sizes)

# --- 步骤 2.5: (新增) 执行并解读显著性检验 ---
message("\n----------------------------------------------------")
message("--- 执行组间HCDR3长度的显著性检验 ---")
message("----------------------------------------------------")
stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = lcdr3_lengths,
  method = "wilcox.test"
)
message("统计检验结果:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\n显著的组间差异 (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("在 p.adj < 0.05 的水平上，未发现显著的组间差异。")
}
message("----------------------------------------------------\n")

# --- 步骤 3: 使用 ggplot2 和 ggpubr 绘制出版级小提琴图 ---

# 定义组间的比较，用于显著性检验
my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

# 定义颜色方案
color_palette <- c("Naive" = "#0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") # 蓝、橙、红

plot4e <- ggplot(lcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  # 绘制小提琴图层
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  # 在小提琴内部添加箱线图
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  # 设置颜色方案
  scale_fill_manual(values = color_palette) +
  
  # 添加显著性检验和P值
  # method = "wilcox.test" (非参数检验) 或 "t.test" (参数检验)
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", # 显示 *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  # 添加样本量标签
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 5, fontface = "bold", color = "black") +
  
  # 设置坐标轴标签和范围
  labs(
    title = "LCDR3 Length Distribution",
    y = "LCDR3 Length (aa)",
    x = "" # 隐藏X轴标题
  ) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 30, 10)) +
  
  # 设置一个简洁、专业的主题
  theme_classic() +
  theme(
    legend.position = "none", # 隐藏图例
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot4e)
ggsave("./results/Figure4e.pdf", plot = plot4e, width = 6, height = 5)
ggsave("./results/Figure4e.png", plot = plot4e, width = 6, height = 5,dpi = 300)

#####Figure4f H_SHM ----------------------
# --- 步骤 0: 安装和加载所有必要的R包 ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr", "scales")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr) 
library(purrr)  
library(scales) 

# --- 步骤 1: 定义文件和组别，并读取所有数据 ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 创建一个包含文件名和对应组别名的列表
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

# --- 步骤 2: 数据处理 - 提取和转换SHM频率 ---
# 检查必需的列是否存在
if (!"MU_FREQ_HEAVY_TOTAL" %in% names(all_data)) {
  stop("错误：合并后的数据中找不到 'MU_FREQ_HEAVY_TOTAL' 列。请确保您的输入文件是正确的、已计算过SHM的文件。")
}

shm_frequencies <- all_data %>%
  # 只保留重链数据 (locus == "IGH")
  #filter(locus == "IGH") %>%
  # 选择我们需要的列
  select(group, MU_FREQ_HEAVY_TOTAL) %>%
  # 移除任何无效的行
  filter(!is.na(MU_FREQ_HEAVY_TOTAL)) %>%
  # 将频率列从字符型转换为数值型
  mutate(shm_freq = as.numeric(MU_FREQ_HEAVY_TOTAL)) %>%
  # 将小数频率转换为百分比
  mutate(shm_freq_percent = shm_freq * 100) %>%
  # 将组别转换为因子，并按指定顺序排列
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

print(head(shm_frequencies))

# --- 步骤 3: 使用 ggplot2 和 ggpubr 绘制出版级小提琴图 ---

# 定义组间的比较，用于显著性检验
my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

# 定义颜色方案
color_palette <- c("Naive" = "#0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") # 蓝、橙、红

plot4f <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
  # 绘制小提琴图层
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  # 在小提琴内部添加箱线图
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  # 设置颜色方案
  scale_fill_manual(values = color_palette) +
  
  # 添加显著性检验和P值
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", # 显示 *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  # 设置坐标轴标签和范围
  labs(
    title = "Heavy Chain Mutation Frequency",
    y = "Heavy chain mutation frequency (%)",
    x = "" # 隐藏X轴标题
  ) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 30, 10)) +
  
  # 设置一个简洁、专业的主题
  theme_classic() +
  theme(
    legend.position = "none", # 隐藏图例
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

# 打印图表到R的绘图窗口
print(plot4f)
ggsave("./results/Figure4f.pdf", plot = plot4f, width = 6, height = 5)
ggsave("./results/Figure4f.png", plot = plot4f, width = 6, height = 5,dpi = 300)


##### FigureS15 ----------

# --- 步骤 0: 安装和加载所有必要的R包 ---
# (与之前相同)
packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr", "scales")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(dplyr); library(readr); library(stringr); library(ggplot2); library(ggpubr); library(purrr); library(scales)

# --- 步骤 1: 加载原始文件 ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# --- 步骤 2: (已修正) 将5个独立样本的数据框放入一个列表中 ---
# 我们为每个样本赋予一个清晰的、用于绘图的名称
list_of_samples <- list(
  "Flu_H1_Mouse1" = flu_H1_mouse1.bcr,
  "Flu_H1_Mouse2" = flu_H1_mouse2.bcr,
  "Flu_H5_Mouse" = flu_H5_mouse.bcr,
  "Naive_Mouse1" = naive_mouse1.bcr,
  "Naive_Mouse2" = naive_mouse2.bcr
)

# 使用 bind_rows 将列表中的数据框合并，并自动创建 'group' 列
all_data <- bind_rows(list_of_samples, .id = "group")

message("所有5个样本的数据已成功合并。")

# --- 步骤 3: 数据处理 - 提取和转换SHM频率 ---
if (!"MU_FREQ_HEAVY_TOTAL" %in% names(all_data)) {
  stop("错误：数据中找不到 'MU_FREQ_HEAVY_TOTAL' 列。")
}

# 定义组别的顺序，以便在X轴上按逻辑排列
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

shm_frequencies <- all_data %>%
  # 恢复locus筛选，确保只分析重链
  filter(locus == "IGH") %>%
  select(group, MU_FREQ_HEAVY_TOTAL) %>%
  filter(!is.na(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM频率数据处理完成。")

# --- 步骤 4: (已修正) 使用 ggplot2 和 ggpubr 绘制5个样本的小提琴图 ---
message("开始绘制小提琴图...")

# 定义组间的比较
# 为了避免图形过于混乱，我们只显示部分关键比较
my_comparisons <- list( 
  c("Naive_Mouse1", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Naive_Mouse1", "Flu_H5_Mouse"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Naive_Mouse2")
)

# 定义新的颜色方案
color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plotS15 <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, size = 5) +
  labs(
    title = "Heavy Chain Mutation Frequency (5 Samples)",
    y = "Heavy chain mutation frequency (%)",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 30, 10)) + # 您可以根据需要调整Y轴范围
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plotS15)
ggsave("./results/FigureS15.pdf", plot = plotS15, width = 8, height = 7)
ggsave("./results/FigureS15.png", plot = plotS15, width = 8, height = 7,dpi = 300)

#####Figure4g L_SHM -----------------------

base_path <- "./cell_reports/VDJ/annotation_results/light_shm_results/"
flu_H1_mouse1.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr.light <- read_tsv(paste0(base_path, "fluH5_mouse_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr.light <- read_tsv(paste0(base_path, "naive_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr.light <- read_tsv(paste0(base_path, "naive_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))

# 1.1 (新增) 在R内存中创建合并后的数据框
flu_H1_mouse.bcr.light <- bind_rows(flu_H1_mouse1.bcr.light, flu_H1_mouse2.bcr.light)
naive_mouse.bcr.light <- bind_rows(naive_mouse1.bcr.light, naive_mouse2.bcr.light)

# 1.2 (已修改) 将3个合并后的数据框放入一个列表中
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr.light,
  "Flu_H5" = flu_H5_mouse.bcr.light,
  "Naive" = naive_mouse.bcr.light
)

# 1.3 使用 bind_rows 将列表中的数据框合并，并自动创建 'group' 列
all_data <- bind_rows(list_of_samples, .id = "group")

# --- 步骤 2: 数据处理 - 提取和转换SHM频率 ---
if (!"MU_FREQ_LIGHT_TOTAL" %in% names(all_data)) {
  stop("错误：数据中找不到 'MU_FREQ_LIGHT_TOTAL' 列。")
}

# 2.1 (已修改) 定义组别的顺序
group_levels <- c("Naive", "Flu_H1", "Flu_H5")

shm_frequencies <- all_data %>%
  select(group, MU_FREQ_LIGHT_TOTAL) %>%
  filter(!is.na(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM频率数据处理完成。")

# --- 步骤 3: (已修正) 使用 ggplot2 和 ggpubr 绘制3个组别的小提琴图 ---
message("开始绘制小提琴图...")

# 3.1 (已修改) 定义组间的比较
my_comparisons <- list( 
  c("Naive", "Flu_H1"), 
  c("Naive", "Flu_H5"),
  c("Flu_H1", "Flu_H5")
)

# 3.2 (已修改) 定义新的颜色方案
color_palette <- c(
  "Naive" = "#0072B2", 
  "Flu_H1" = "#E69F00", 
  "Flu_H5" = "#D55E00"
)

plot4g <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, size = 6) +
  labs(
    title = "Light Chain Mutation Frequency",
    y = "Light chain mutation frequency (%)",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) + # 调整Y轴范围和刻度
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text( hjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot4g)
ggsave("./results/Figure4g.pdf", plot = plot4g, width = 6, height = 5)
ggsave("./results/Figure4g.png", plot = plot4g, width = 6, height = 5,dpi = 300)


##### FigureS16 ------------------

base_path <- "./cell_reports/VDJ/annotation_results/light_shm_results/"
flu_H1_mouse1.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr.light <- read_tsv(paste0(base_path, "fluH5_mouse_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr.light <- read_tsv(paste0(base_path, "naive_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr.light <- read_tsv(paste0(base_path, "naive_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))

# 我们为每个样本赋予一个清晰的、用于绘图的名称
list_of_samples <- list(
  "Flu_H1_Mouse1" = flu_H1_mouse1.bcr.light,
  "Flu_H1_Mouse2" = flu_H1_mouse2.bcr.light,
  "Flu_H5_Mouse" = flu_H5_mouse.bcr.light,
  "Naive_Mouse1" = naive_mouse1.bcr.light,
  "Naive_Mouse2" = naive_mouse2.bcr.light
)

# 使用 bind_rows 将列表中的数据框合并，并自动创建 'group' 列
all_data <- bind_rows(list_of_samples, .id = "group")

message("所有5个样本的数据已成功合并。")

# --- 步骤 3: 数据处理 - 提取和转换SHM频率 ---
if (!"MU_FREQ_LIGHT_TOTAL" %in% names(all_data)) {
  stop("错误：数据中找不到 'MU_FREQ_LIGHT_TOTAL' 列。")
}

# 定义组别的顺序，以便在X轴上按逻辑排列
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

shm_frequencies <- all_data %>%
  # 恢复locus筛选，确保只分析重链
  #filter(locus == "IGH") %>%
  select(group, MU_FREQ_LIGHT_TOTAL) %>%
  filter(!is.na(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM频率数据处理完成。")

# --- 步骤 4: (已修正) 使用 ggplot2 和 ggpubr 绘制5个样本的小提琴图 ---
message("开始绘制小提琴图...")

# 定义组间的比较
# 为了避免图形过于混乱，我们只显示部分关键比较
my_comparisons <- list( 
  c("Naive_Mouse1", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Naive_Mouse1", "Flu_H5_Mouse"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Naive_Mouse2")
)

# 定义新的颜色方案
color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plotS16 <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, size = 5) +
  labs(
    title = "Light Chain Mutation Frequency (5 Samples)",
    y = "Light chain mutation frequency (%)",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 30, 10)) + # 您可以根据需要调整Y轴范围
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text( hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plotS16)
ggsave("./results/FigureS16.pdf", plot = plotS16, width = 8, height = 7)
ggsave("./results/FigureS16.png", plot = plotS16, width = 8, height = 7,dpi = 300)

#####Figure4h------
## 0) 取对象与 meta.data
obj <- B_cell_subset_flu_mouse.integrated

## 1) 在 meta.data 中新增 SHM 分组（直接写回对象）
#    阈值：>0.02 = high；(0,0.02] = low；==0 = none
obj$MU_FREQ_HEAVY_TOTAL <- suppressWarnings(as.numeric(obj$MU_FREQ_HEAVY_TOTAL))

obj$SHM_group <- dplyr::case_when(
  obj$MU_FREQ_HEAVY_TOTAL > 0.02 ~ "high",
  obj$MU_FREQ_HEAVY_TOTAL > 0    & obj$MU_FREQ_HEAVY_TOTAL <= 0.02 ~ "low",
  obj$MU_FREQ_HEAVY_TOTAL == 0   ~ "none",
  TRUE ~ NA_character_
) %>% factor(levels = c("high","low","none"))

table(obj$SHM_group)  # 总览一下数量

h1_GC <- subset(
  obj,
  subset = grepl("^flu_H1", orig.ident)&
    B_cell_subpopulations %in% c("GC")
)

## 3) 仅取 high/low 两组做 DEG
h1_GC <- subset(h1_GC, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h1_GC$SHM_group)

Idents(h1_GC) <- "SHM_group"

#（可选）确保使用 RNA assay 的归一化数据
# DefaultAssay(h1) <- "RNA"   # 若你是在 integrated assay 上做，也可保持默认

## 4) 差异表达：high vs low
deg_high_vs_low <- FindMarkers(
  h1_GC,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  # 有多个 flu_H1_* 样本时建议保留
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 阈值（可改）
LFC_cut  <- 0.25    # |log2FC| 阈值（可改）
LABEL_N  <- 25       # Up/Down 各标注多少个基因（可改）

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 的列名
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), # 避免 -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

# 要标注的代表基因（各取最显著的 LABEL_N 个）
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)  # 让标签位置可复现
# 挑选要标注的基因（各取 FDR 最小的和 |log2FC| 最大的）
LABEL_N <- 25
labels_df <- dplyr::bind_rows(
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2))
) %>% dplyr::distinct(gene, .keep_all = TRUE)

plot4h <- ggplot(df, aes(x = log2FC, y = negLog10)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.9) +
  geom_hline(yintercept = -log10(FDR_cut), linetype = "dotted") +
  geom_vline(xintercept = c(-LFC_cut, LFC_cut), linetype = "dotted") +
  geom_label_repel(
    data = labels_df, aes(label = gene),
    size = 4, label.r = unit(0.15, "lines"), label.size = 0.2,
    box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, max.overlaps = Inf   # 强制不丢标签
  ) +
  scale_color_manual(values = c("Up"="#E68A00","Down"="#3C6EB4","NS"="grey80")) +
  labs(title = "High vs Low SHM (flu_H1 GC)",
       x = "Log2 Foldchange", y = expression(-log[10]~FDR), color = NULL) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(plot.title = element_text(face="bold", hjust=0.5),
        axis.title = element_text(face="bold"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 28, 10, 10)) +
  coord_cartesian(clip = "off") +
  annotate("text",
           x = min(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 0, vjust = -0.5, label = paste0(n_down, " DOWN"), fontface = 2) +
  annotate("text",
           x = max(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 1, vjust = -0.5, label = paste0(n_up, " UP"), fontface = 2)

print(plot4h)
ggsave("./results/Figure4h.pdf", plot = plot4h, width = 8, height = 7)
ggsave("./results/Figure4h.png", plot = plot4h, width = 8, height = 7,dpi = 300)

##### Figure4i------
h1_PB <- subset(
  obj,
  subset = grepl("^flu_H1", orig.ident)&
    B_cell_subpopulations %in% c("PB")
)

## 3) 仅取 high/low 两组做 DEG
h1_PB <- subset(h1_PB, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h1_PB$SHM_group)

Idents(h1_PB) <- "SHM_group"

#（可选）确保使用 RNA assay 的归一化数据
# DefaultAssay(h1) <- "RNA"   # 若你是在 integrated assay 上做，也可保持默认

## 4) 差异表达：high vs low
deg_high_vs_low <- FindMarkers(
  h1_PB,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  # 有多个 flu_H1_* 样本时建议保留
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 阈值（可改）
LFC_cut  <- 0.25    # |log2FC| 阈值（可改）
LABEL_N  <- 25       # Up/Down 各标注多少个基因（可改）

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 的列名
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), # 避免 -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

# 要标注的代表基因（各取最显著的 LABEL_N 个）
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)  # 让标签位置可复现
# 挑选要标注的基因（各取 FDR 最小的和 |log2FC| 最大的）
LABEL_N <- 25
labels_df <- dplyr::bind_rows(
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2))
) %>% dplyr::distinct(gene, .keep_all = TRUE)

plot4i <- ggplot(df, aes(x = log2FC, y = negLog10)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.9) +
  geom_hline(yintercept = -log10(FDR_cut), linetype = "dotted") +
  geom_vline(xintercept = c(-LFC_cut, LFC_cut), linetype = "dotted") +
  geom_label_repel(
    data = labels_df,
    aes(label = gene),
    size = 4, label.r = unit(0.15, "lines"), label.size = 0.2,
    box.padding = 0.35, point.padding = 0.2, min.segment.length = 0
  ) +
  scale_color_manual(values = c("Up" = "#E68A00", "Down" = "#3C6EB4", "NS" = "grey80")) +
  labs(
    title = "High vs Low SHM (flu_H1 PB)",
    x = "Log2 Foldchange",
    y = expression(-log[10]~FDR),
    color = NULL
  ) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  annotate("text",
           x = min(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 0, vjust = -0.5, label = paste0(n_down, " DOWN"), fontface = 2) +
  annotate("text",
           x = max(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 1, vjust = -0.5, label = paste0(n_up, " UP"), fontface = 2)

print(plot4i)
ggsave("./results/Figure4i.pdf", plot = plot4i, width = 8, height = 7)
ggsave("./results/Figure4i.png", plot = plot4i, width = 8, height = 7,dpi = 300)

##### Figure4j------
h5_GC <- subset(
  obj,
  subset = grepl("^flu_H5", orig.ident)&
    B_cell_subpopulations %in% c("GC")
)

## 3) 仅取 high/low 两组做 DEG
h5_GC <- subset(h5_GC, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h5_GC$SHM_group)

Idents(h5_GC) <- "SHM_group"

## 4) 差异表达：high vs low
deg_high_vs_low <- FindMarkers(
  h5_GC,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 阈值（可改）
LFC_cut  <- 0.25    # |log2FC| 阈值（可改）
LABEL_N  <- 25       # Up/Down 各标注多少个基因（可改）

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 的列名
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), # 避免 -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

# 要标注的代表基因（各取最显著的 LABEL_N 个）
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)  # 让标签位置可复现
# 挑选要标注的基因（各取 FDR 最小的和 |log2FC| 最大的）
LABEL_N <- 25
labels_df <- dplyr::bind_rows(
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2))
) %>% dplyr::distinct(gene, .keep_all = TRUE)

plot4j <- ggplot(df, aes(x = log2FC, y = negLog10)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.9) +
  geom_hline(yintercept = -log10(FDR_cut), linetype = "dotted") +
  geom_vline(xintercept = c(-LFC_cut, LFC_cut), linetype = "dotted") +
  geom_label_repel(
    data = labels_df,
    aes(label = gene),
    size = 4, label.r = unit(0.15, "lines"), label.size = 0.2,
    box.padding = 0.35, point.padding = 0.2, min.segment.length = 0
  ) +
  scale_color_manual(values = c("Up" = "#E68A00", "Down" = "#3C6EB4", "NS" = "grey80")) +
  labs(
    title = "High vs Low SHM (flu_h5 GC)",
    x = "Log2 Foldchange",
    y = expression(-log[10]~FDR),
    color = NULL
  ) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  annotate("text",
           x = min(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 0, vjust = -0.5, label = paste0(n_down, " DOWN"), fontface = 2) +
  annotate("text",
           x = max(df$log2FC, na.rm = TRUE), y = max(df$negLog10, na.rm = TRUE),
           hjust = 1, vjust = -0.5, label = paste0(n_up, " UP"), fontface = 2)

print(plot4j)
ggsave("./results/Figure4j.pdf", plot = plot4j, width = 8, height = 7)
ggsave("./results/Figure4j.png", plot = plot4j, width = 8, height = 7,dpi = 300)

##### Figure4k------
h5_PB <- subset(
  obj,
  subset = grepl("^flu_H5", orig.ident)&
    B_cell_subpopulations %in% c("PB")
)

## 3) 仅取 high/low 两组做 DEG
h5_PB <- subset(h5_PB, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h5_PB$SHM_group)

Idents(h5_PB) <- "SHM_group"

## 4) 差异表达：high vs low
deg_high_vs_low <- FindMarkers(
  h5_PB,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  # 有多个 flu_h5_* 样本时建议保留
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 阈值（可改）
LFC_cut  <- 0.25    # |log2FC| 阈值（可改）
LABEL_N  <- 25       # Up/Down 各标注多少个基因（可改）

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 的列名
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), # 避免 -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

# 要标注的代表基因（各取最显著的 LABEL_N 个）
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N))

set.seed(1)  # 让标签位置可复现

# 挑选要标注的基因（各取 FDR 最小的和 |log2FC| 最大的）
LABEL_N <- 25
labels_df <- dplyr::bind_rows(
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Up")   %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(FDR)        %>% head(ceiling(LABEL_N/2)),
  df %>% dplyr::filter(sig == "Down") %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>% head(floor(LABEL_N/2))
) %>% dplyr::distinct(gene, .keep_all = TRUE)

plot4k <- ggplot(df, aes(x = log2FC, y = negLog10)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 1.9) +
  geom_hline(yintercept = -log10(FDR_cut), linetype = "dotted") +
  geom_vline(xintercept = c(-LFC_cut, LFC_cut), linetype = "dotted") +
  geom_label_repel(
    data = labels_df, aes(label = gene),
    size = 4, label.r = unit(0.15, "lines"), label.size = 0.2,
    box.padding = 0.4, point.padding = 0.25,
    min.segment.length = 0, segment.size = 0.3,
    max.overlaps = Inf                     # ← 关键：不因重叠而丢标签
  ) +
  scale_color_manual(values = c("Up"="#E68A00","Down"="#3C6EB4","NS"="grey80")) +
  labs(title = "High vs Low SHM (flu_h5 PB)",
       x = "Log2 Foldchange", y = expression(-log[10]~FDR), color = NULL) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 24, 10, 10)  # 右侧留白给标签
  ) +
  coord_cartesian(clip = "off") +         # ← 防止标签被裁切
  annotate("text",
           x = min(df$log2FC, na.rm=TRUE), y = max(df$negLog10, na.rm=TRUE),
           hjust = 0, vjust = -0.5, label = paste0(n_down, " DOWN"), fontface = 2) +
  annotate("text",
           x = max(df$log2FC, na.rm=TRUE), y = max(df$negLog10, na.rm=TRUE),
           hjust = 1, vjust = -0.5, label = paste0(n_up, " UP"), fontface = 2)

print(plot4k)
ggsave("./results/Figure4k.pdf", plot = plot4k, width = 8, height = 7)
ggsave("./results/Figure4k.png", plot = plot4k, width = 8, height = 7,dpi = 300)

###Figure5 ------

write.csv(md,'./results/md.csv',row.names = T)
md1 <- md %>%
  mutate(sequence_id = as.character(sequence_id)) %>%     # 保险：转成字符
  filter(!is.na(sequence_id) & trimws(sequence_id) != "" & toupper(sequence_id) != "NA")
md1 <- md1 %>%
  mutate(v_call_10x = as.character(v_call_10x)) %>%     # 保险：转成字符
  filter(!is.na(v_call_10x) & trimws(v_call_10x) != "" & toupper(v_call_10x) != "NA")
write.csv(md1,'./results/md1.csv',row.names = T)
write_tsv(md1,'./results/md1.tsv')

# 建议在 conda/mamba 的 immcantation 环境或 R>=4.2 下
install.packages(c("tidyverse","ggplot2","ggtree","cowplot"))
install.packages("alakazam")   # CRAN 有包
install.packages("shazam")
install.packages("dowser")     # CRAN 有包
# 任选：多色标叠加（如果想把SHM分组画成外环颜色）
# install.packages("ggnewscale")

library(alakazam)
library(shazam)
library(dowser)
library(ggtree)
library(cowplot)
#R 里检查可执行文件路径
igphyml_exec <- Sys.which("igphyml")
if (igphyml_exec == "" || !file.exists(igphyml_exec)) {
  # 如果没在 PATH 里，请手动填绝对路径，例如：
  igphyml_exec <- "/usr/local/bin/igphyml"   # <- 改成你的实际路径
}
system2(igphyml_exec, args = "-h")   # 看看是否能正常打印帮助

#(一)
# 1.1 读取
db <- readr::read_tsv("./results/md1.tsv", progress = TRUE, show_col_types = FALSE)

# 1.2 样本特异性分组（从 orig.ident 解析）
library(dplyr)
library(stringr)

db <- db %>%
  mutate(
    # 从 orig.ident 解析样本分组
    specificity = case_when(
      str_detect(orig.ident, "^flu[_-]?H1") ~ "flu_H1",
      str_detect(orig.ident, "^flu[_-]?H5") ~ "flu_H5",
      str_detect(orig.ident, "^naive")      ~ "naive",
      TRUE ~ "other"
    ),
    
    # 计算突变比例并分箱
    MU_pct = MU_FREQ_HEAVY_TOTAL * 100,
    MU_bin = case_when(
      is.na(MU_pct) ~ NA_character_,
      MU_pct > 5    ~ ">5%",
      MU_pct >= 2   ~ "2%-5%",
      TRUE          ~ "<2%"
    )
  ) %>%
  
  # 常规质量过滤（可根据你表格的列名进一步调整）
  filter(
    !is.na(sequence_alignment),
    !is.na(germline_alignment),
    productive %in% c(TRUE, 1)
  )
write.csv(db,'./results/db.csv',row.names = T)

#(二)
# 1) 生成 presence 矩阵（每个 clone_id 是否出现在各组）
pres <- db %>%
  filter(specificity %in% c("flu_H1","flu_H5","naive")) %>%
  distinct(clone_id, specificity) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = specificity, values_from = present, values_fill = 0L)

# 2) 组内“特异”克隆
clones_h1_only    <- pres %>% filter(flu_H1==1, flu_H5==0, naive==0) %>% pull(clone_id)
clones_h5_only    <- pres %>% filter(flu_H1==0, flu_H5==1, naive==0) %>% pull(clone_id)
clones_naive_only <- pres %>% filter(flu_H1==0, flu_H5==0, naive==1) %>% pull(clone_id)

# 3) H1–H5 “共享”克隆
#    3.1 严格：H1==1 & H5==1 且 naive==0
clones_h1h5_shared_strict   <- pres %>% filter(flu_H1==1, flu_H5==1, naive==0) %>% pull(clone_id)
#    3.2 宽松：只要 H1==1 & H5==1（不管 naive）
clones_h1h5_shared_inclusive <- pres %>% filter(flu_H1==1, flu_H5==1) %>% pull(clone_id)

lengths(list(
  H1_only    = clones_h1_only,
  H5_only    = clones_h5_only,
  naive_only = clones_naive_only,
  H1H5_shared_strict = clones_h1h5_shared_strict,
  H1H5_shared_inclusive = clones_h1h5_shared_inclusive
))

# 4) filter: H1–H5 “共享”克隆
min_seqs_per_clone <- 5  # 需要更严格就调大
clone_sizes <- db %>% dplyr::count(clone_id, name="n")

keep_by_size <- function(ids) {
  tibble(clone_id = ids) %>% inner_join(clone_sizes, by="clone_id") %>%
    filter(n >= min_seqs_per_clone) %>% pull(clone_id)
}

clones_h1_only    <- keep_by_size(clones_h1_only)
clones_h5_only    <- keep_by_size(clones_h5_only)
clones_naive_only <- keep_by_size(clones_naive_only)
clones_h1h5_shared_strict    <- keep_by_size(clones_h1h5_shared_strict)
clones_h1h5_shared_inclusive <- keep_by_size(clones_h1h5_shared_inclusive)

#(三) formatClones()

library(dplyr)
library(stringr)

keep_traits <- c("c_call","B_cell_subpopulations","specificity","MU_bin","MU_pct","orig.ident")
keep_text   <- c("c_call","B_cell_subpopulations","specificity","MU_bin","orig.ident")

format_for <- function(db, ids) {
  db %>%
    filter(clone_id %in% ids) %>%
    # 1) 选用 D 区屏蔽后的胚系；若列不存在则回退到 germline_alignment
    mutate(germ_use = dplyr::coalesce(germline_alignment_d_mask, germline_alignment)) %>%
    # 2) 规范化：大小写、去空白、统一缺口符号
    mutate(germ_use = stringr::str_replace_all(stringr::str_to_upper(trimws(germ_use)), "-", ".")) %>%
    # text_fields 必须为 character
    mutate(across(all_of(keep_text), as.character)) %>%
    dowser::formatClones(
      seq        = "sequence_alignment",
      germ       = "germ_use",        # <— 用规范化后的列
      clone      = "clone_id",
      id         = "sequence_id",
      v_call     = "v_call",
      j_call     = "j_call",
      junc_len   = "junction_length",
      traits     = keep_traits,
      text_fields= keep_text,
      minseq     = 2,
      collapse   = TRUE
    )
}

clones_H1_only    <- format_for(db, clones_h1_only)
clones_H5_only    <- format_for(db, clones_h5_only)
clones_naive_only <- format_for(db, clones_naive_only)
clones_H1H5_shared_strict    <- format_for(db, clones_h1h5_shared_strict)
clones_H1H5_shared_inclusive <- format_for(db, clones_h1h5_shared_inclusive)

#(四)  IgPhyML
build_with_igphyml <- function(clones_tbl, title="", exec=igphyml_exec,
                               locus="IGH", region="V", model="HLP19",
                               nproc=parallel::detectCores()) {
  if (nrow(clones_tbl) == 0) return(NULL)
  dowser::getTrees(
    clones_tbl,
    build     = "igphyml",
    exec      = exec,
    locus     = locus,
    region    = region,
    model     = model,
    reorient  = TRUE,   # 以胚系为根方向
    ref       = "Germline",
    nproc     = nproc,
    rm_temp   = TRUE,   # 清理临时文件
    collapse  = TRUE
  )
}

trees_H1_only    <- build_with_igphyml(clones_H1_only,    "H1-specific")
trees_H5_only    <- build_with_igphyml(clones_H5_only,    "H5-specific")
#trees_naive_only <- build_with_igphyml(clones_naive_only, "naive-specific")
#trees_H1H5_shared_strict    <- build_with_igphyml(clones_H1H5_shared_strict,    "H1–H5 shared (strict)")
#trees_H1H5_shared_inclusive <- build_with_igphyml(clones_H1H5_shared_inclusive, "H1–H5 shared (inclusive)")

#(五) plot

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(purrr)
  library(tibble)
  library(fs)
})

## 兼容旧版：某些依赖会用到 is.waive()
if (!exists("is.waive", mode = "function")) {
  is.waive <- function(x) inherits(x, "waiver")
}

isotype_pal <- c(
  "IGHM"   = "#3333aa",
  "IGHD"   = "#0072B2",
  "IGHG1"  = "#E41A1C",
  "IGHG2B" = "#8E24AA",
  "IGHG2C" = "#FF7F0E",
  "IGHG3"  = "#2ECC40",
  "IGHA"   = "#FFD700",
  "IGHE"   = "#A0522D",
  "NA"     = "#BDBDBD"
)

## 更高对比度的 SHM 颜色（描边）
mu_pal <- c(
  "<2%"   = "#0055FF",
  "2%-5%" = "#FFB000",
  ">5%"   = "#CC0000"
)

## 21-25 是带描边的形状；GC/Bmem圆、PB三角、MZ菱形
shape_pal <- c(GC=21, PB=24, Bmem=21, MZ=22)

alignment_len_for_tree <- function(tr, db) {
  seqs <- db %>% filter(sequence_id %in% tr$tip.label) %>% pull(sequence_alignment)
  s <- seqs[!is.na(seqs)][1]
  if (length(s) == 0 || is.na(s)) return(NA_integer_)
  nchar(gsub("[\\.-]", "", s))
}

## tip 注释（按唯一 cell_id 计数；新增 tip_label）
## tip 注释（按唯一 cell_id 计数；新增 tip_label）
tip_anno_for_tree <- function(tr, db){
  tips <- tr$tip.label
  
  ## 选出该树对应的 clone（用于扩增统计）
  cid <- db %>% 
    dplyr::filter(sequence_id %in% tips) %>%
    dplyr::count(clone_id, sort = TRUE) %>% 
    dplyr::slice_head(n = 1) %>% 
    dplyr::pull(clone_id)
  
  ## 同一序列（alignment）上的 cell_id 去重后计数 -> 扩增数
  exp_tbl <- db %>% 
    dplyr::filter(clone_id == cid) %>%
    dplyr::group_by(sequence_alignment) %>%
    dplyr::summarise(exp_n = dplyr::n_distinct(cell_id), .groups = "drop")
  
  ## —— 关键：确保每个 sequence_id（= label）只保留一行，避免 join 放大
  anno <- db %>% 
    dplyr::filter(sequence_id %in% tips) %>%
    dplyr::group_by(sequence_id) %>%
    dplyr::summarise(
      c_call = dplyr::first(c_call),
      B_cell_subpopulations = dplyr::first(B_cell_subpopulations),
      MU_bin = dplyr::first(MU_bin),
      MU_pct = dplyr::first(MU_pct),
      sequence_alignment = dplyr::first(sequence_alignment),
      cell_id = dplyr::first(cell_id),
      .groups = "drop"
    ) %>%
    dplyr::left_join(exp_tbl, by = "sequence_alignment") %>%
    dplyr::mutate(
      exp_n = ifelse(is.na(exp_n), 1L, exp_n),
      B_cell_subpopulations = forcats::fct_drop(
        factor(B_cell_subpopulations, levels = c("GC","PB","Bmem","MZ"))
      ),
      label = as.character(sequence_id),
      c_call = as.character(c_call),
      MU_bin = factor(as.character(MU_bin),
                      levels = c("<2%","2%-5%",">5%"), ordered = TRUE),
      tip_label = dplyr::coalesce(as.character(cell_id), as.character(sequence_id))
    )
  
  return(anno)
}

## 主绘图函数
lineage_plot <- function(
    tr, title = NULL, db,
    show_tiplab = FALSE,
    show_branch_numbers = TRUE,
    show_tip_numbers = FALSE,
    branch_scale = 1,
    branch_transform = c("none","sqrt","log1p"),
    tree_layout = c("rectangular","slanted")   # 新增：布局可选
){
  branch_transform <- match.arg(branch_transform)
  tree_layout <- match.arg(tree_layout)
  
  stopifnot(inherits(tr, "phylo"))
  has_len <- !is.null(tr$edge.length) && length(tr$edge.length) > 0
  
  ## —— 等比例/非线性缩放 edge.length（若存在）
  if (has_len) {
    el <- tr$edge.length
    if (branch_transform == "sqrt")  el <- sqrt(pmax(el, 0))
    if (branch_transform == "log1p") el <- log1p(pmax(el, 0))
    if (is.numeric(branch_scale) && branch_scale != 1) el <- el * branch_scale
    tr$edge.length <- el
  }
  
  anno <- tip_anno_for_tree(tr, db)
  alnL <- alignment_len_for_tree(tr, db)
  
  ## —— 只画一层枝条；肘线时要明确使用分支长度
  p <- ggtree(
    tr,
    ladderize = TRUE,
    layout = tree_layout,
    branch.length = if (has_len) "branch.length" else "none"
  ) +
    theme_tree2() +
    ggtitle(ifelse(is.null(title), "", title))
  
  ## —— 合并注释；再按 node 去重，防止 join 放大导致“编号重复”
  p$data <- p$data %>%
    dplyr::left_join(anno %>% dplyr::distinct(label, .keep_all = TRUE), by = "label") %>%
    dplyr::group_by(node) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup()
  
  ## tips：填充=isotype；描边=SHM；大小=扩增；形状=亚群
  p <- p + ggtree::geom_point2(
    ggplot2::aes(subset = isTip, size = exp_n, shape = B_cell_subpopulations,
                 fill = c_call, colour = MU_bin),
    stroke = 1.1, na.rm = TRUE
  ) +
    scale_fill_manual(values = isotype_pal, name = "Isotype", na.translate = FALSE) +
    scale_colour_manual(values = mu_pal,    name = "SHM",     na.translate = FALSE) +
    scale_shape_manual(
      values = shape_pal,
      breaks = levels(stats::na.omit(anno$B_cell_subpopulations)),
      drop   = TRUE, na.translate = FALSE,
      name   = "B cell subpop"
    ) +
    scale_size_area(name = "Expansion (cells)", max_size = 7) +
    guides(
      fill   = guide_legend(override.aes = list(shape = 21, size = 4, colour = "grey40")),
      colour = guide_legend(override.aes = list(shape = 21, size = 4, fill = NA)),
      shape  = guide_legend(override.aes = list(size = 4))
    )
  
  # ## —— tip 编号：这里用唯一的 node 编号（不再重复）
  # if (isTRUE(show_tip_numbers)) {
  #   p <- p + ggtree::geom_text2(
  #     ggplot2::aes(subset = isTip, label = .data$node),
  #     size = 5, fontface = "bold", colour = "black",
  #     vjust = 0.5, hjust = 0.5, check_overlap = TRUE
  #   )
  # }
  
  ## 胚系（若存在）
  p <- p + ggtree::geom_point2(
    data = function(df) dplyr::filter(df, label %in% c("Germline","Germline_IGH")),
    ggplot2::aes(x = x, y = y), shape = 21, fill = "black", colour = "black", size = 3
  )
  
  if (isTRUE(show_tiplab)) p <- p + ggtree::geom_tiplab(size = 2.8)
  
  ## 枝上数字（替换位点数）：使用“缩放后的” branch.length × alignment_len
  if (isTRUE(show_branch_numbers) && has_len && is.finite(alnL)) {
    td <- ggtree::fortify(tr) %>%
      dplyr::mutate(
        mut_n = ifelse(is.na(branch.length), NA_real_,
                       round(branch.length * alnL, 0))
      )
    p <- p + ggplot2::geom_text(
      data = td,
      ggplot2::aes(x = x, y = y, label = ifelse(is.na(mut_n), "", mut_n)),
      size = 3.5, vjust = -0.35, hjust = -0.1, colour = "grey30",
      inherit.aes = FALSE
    )
  }
  
  ## 朝向 + 去坐标轴
  p <- p + ggplot2::coord_flip() + ggplot2::scale_x_reverse() +
    ggplot2::theme(axis.line  = element_blank(),
                   axis.text  = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank())
  
  return(p)
}

extract_phylo_col <- function(airr_tbl) {
  if (is.null(airr_tbl) || nrow(airr_tbl) == 0) return(tibble())
  if (!"trees" %in% names(airr_tbl))
    stop("这个 airrTrees 没有 `trees` 列，请用 str(trees_*) 查看真实列名。")
  tibble(
    clone_id = airr_tbl$clone_id,
    tree     = purrr::map(airr_tbl$trees, function(z){
      if (is.list(z) && "tree" %in% names(z) && inherits(z$tree, "phylo")) return(z$tree)
      if (inherits(z, "phylo")) return(z)
      NULL
    })
  ) %>% filter(!purrr::map_lgl(tree, is.null))
}

tH1 <- extract_phylo_col(trees_H1_only) %>% mutate(group = "H1")
tH5 <- extract_phylo_col(trees_H5_only) %>% mutate(group = "H5")
all_trees_df <- bind_rows(tH1, tH5)

clone_meta <- db %>%
  filter(specificity %in% c("flu_H1","flu_H5")) %>%
  count(clone_id, specificity, name = "n") %>%
  mutate(group = recode(specificity, flu_H1 = "H1", flu_H5 = "H5")) %>%
  select(clone_id, group, n)

df <- all_trees_df %>% left_join(clone_meta, by = c("clone_id","group"))

outdir <- "./results/lineage_by_group1"
fs::dir_create(outdir, recurse = TRUE)

make_group_pdf <- function(df_group, grp,
                           show_tiplab = FALSE,
                           show_branch_numbers = TRUE,
                           show_tip_numbers = FALSE,
                           branch_scale = 1,
                           branch_transform = c("none","sqrt","log1p")) {
  
  sub <- df_group %>% filter(group == grp) %>% arrange(desc(n))
  if (nrow(sub) == 0) { message("组 ", grp, " 没有树"); return(invisible(NULL)) }
  
  outfile <- file.path(outdir, sprintf("%s_sorted_by_size.pdf", grp))
  pdf(outfile, width = 7.5, height = 10, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  
  for (i in seq_len(nrow(sub))) {
    cid  <- sub$clone_id[i]
    size <- sub$n[i]
    tr   <- sub$tree[[i]]
    ttl  <- if (!is.na(size)) sprintf("%s | clone %s (n=%d)", grp, cid, size)
    else sprintf("%s | clone %s", grp, cid)
    
    p <- lineage_plot(
      tr, title = ttl, db = db,
      show_tiplab = show_tiplab,
      show_branch_numbers = show_branch_numbers,
      show_tip_numbers = show_tip_numbers,
      branch_scale = branch_scale,
      branch_transform = branch_transform
    )
    print(p)
  }
  message("✅ 写出：", outfile)
}

## 例：默认不缩放；如需“整体压缩 0.7 倍”，把 branch_scale=0.7
make_group_pdf(df, "H1",
               show_tiplab = FALSE, show_branch_numbers = TRUE, show_tip_numbers = FALSE,
               branch_scale = 0.5, branch_transform = "none")
make_group_pdf(df, "H5",
               show_tiplab = FALSE, show_branch_numbers = TRUE, show_tip_numbers = FALSE,
               branch_scale = 0.5, branch_transform = "none")



###Figure6 ------
# data filter
db_igh <- db %>%
  filter(locus == "IGH", productive == TRUE) %>%
  mutate(across(c(clone_id, cell_id, c_call,
                  sequence_alignment, germline_alignment, germline_alignment_d_mask),
                ~ as.character(.)))

# 每个细胞在同一克隆只留一条代表序列（按 umi_count 最大）
db_igh <- db_igh %>%
  group_by(specificity, clone_id, cell_id) %>%
  slice_max(order_by = coalesce(umi_count, 1), n = 1, with_ties = FALSE) %>%
  ungroup()

#（可选）确保同一 clone 的胚系一致；不一致的拆成子克隆
db_igh <- db_igh %>%
  mutate(subclone_id = interaction(clone_id, germline_alignment_d_mask, drop = TRUE))

#去掉没有序列 ID 的行
db_igh <- db_igh %>% filter(!is.na(sequence_id) & sequence_id != "NA")

library(dplyr)
library(tidyr)

db_igh_dedup_bc_clone <- db_igh %>%
  arrange(desc(umi_count), desc(consensus_count)) %>%
  group_by(barcode_seurat, clone_id) %>%
  dplyr::slice(1) %>%         # ← 指定 dplyr::slice
  ungroup()

# 查看结果
head(db_igh_dedup_bc_clone)
library(dplyr)
bad <- db_igh_dedup_bc_clone %>%
  filter(locus == "IGH", productive) %>%
  mutate(across(c(clone_id, sequence_alignment, germline_alignment,
                  germline_alignment_d_mask), as.character)) %>%
  group_by(clone_id) %>%
  summarise(
    n_germ_mask = n_distinct(germline_alignment_d_mask),
    n_germ      = n_distinct(germline_alignment),
    .groups = "drop") %>%
  filter(n_germ_mask > 1 | n_germ > 1)
nrow(bad)      # 有问题的 clone 数
head(bad, 10)  # 看看样子

library(dplyr)

db_igh_fix <- db_igh_dedup_bc_clone %>%
  filter(locus == "IGH", productive) %>%
  mutate(across(c(clone_id, germline_alignment, germline_alignment_d_mask), as.character)) %>%
  # 以“clone_id × D-遮蔽后的胚系”细分 clone
  mutate(clone_id = interaction(clone_id, germline_alignment_d_mask, drop = TRUE)) %>%
  # 用遮蔽 D 的胚系替换原始胚系（BuildTrees 检查这一列的一致性）
  mutate(germline_alignment = germline_alignment_d_mask)

# 自检：每个 clone 的胚系现在是否唯一
db_igh_fix %>%
  group_by(clone_id) %>%
  summarise(n_germ = n_distinct(germline_alignment), .groups = "drop") %>%
  summarise(all(n_germ == 1))
# 应返回 TRUE

db_igh_fix <- db_igh_fix %>%
  mutate(
    # 用细胞条码+克隆号拼一个短ID，并限制长度，避免再超长
    db_short = str_trunc(paste0(barcode_seurat, "_", clone_id), 40, ellipsis = "")
  )

names(db_igh_fix)[names(db_igh_fix)=="sequence_id"]   <- "SEQUENCE_ID"
names(db_igh_fix)[names(db_igh_fix)=="barcode_seurat"]<- "CELL_ID"   # 或你想用的名字
readr::write_tsv(db_igh_fix, "db_igh_fix.tsv")

write_tsv(db_igh_fix,'./results/db_igh_fix.tsv')

##### Figure6a Clonotype counts------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggvenn)   # 如未安装：install.packages("ggvenn")
})

# —— 可选：标准化 CDR3，避免大小写/尾部*影响交集
canon_cdr3 <- function(x){
  x %>% as.character() %>% str_trim() %>% str_replace("\\*+$","") %>% toupper()
}

stopifnot(all(c("specificity","cdr3_aa") %in% names(db_igh_fix)))

db_igh_fix <- db_igh_fix %>%
  mutate(
    specificity   = case_when(
      str_to_lower(specificity) == "naive"  ~ "naive",
      str_to_lower(specificity) == "flu_h1" ~ "flu_H1",
      str_to_lower(specificity) == "flu_h5" ~ "flu_H5",
      TRUE ~ specificity
    ),
    cdr3_aa = canon_cdr3(cdr3_aa)
  )

# 1) 提取三组各自的“唯一 CDR3”
naive_clones  <- db_igh_fix %>% filter(specificity == "naive")  %>% distinct(cdr3_aa) %>% pull(cdr3_aa)
flu_h1_clones <- db_igh_fix %>% filter(specificity == "flu_H1") %>% distinct(cdr3_aa) %>% pull(cdr3_aa)
flu_h5_clones <- db_igh_fix %>% filter(specificity == "flu_H5") %>% distinct(cdr3_aa) %>% pull(cdr3_aa)

# 2) 统计各组唯一 CDR3 数
cat("Naive 组的克隆型数:",  length(naive_clones),  "\n")
cat("Flu_H1 组的克隆型数:", length(flu_h1_clones), "\n")
cat("Flu_H5 组的克隆型数:", length(flu_h5_clones), "\n")

# 3) 画 Venn（集合基数）
venn_list <- list(Naive = naive_clones, Flu_H1 = flu_h1_clones, Flu_H5 = flu_h5_clones)

plot6a <- ggvenn(
  venn_list,
  fill_color    = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size   = 0.6,
  set_name_size = 6,      # 增加集合名称的字体大小
  text_size     = 6,      # 增加韦恩图内部数字的字体大小
  show_percentage = FALSE
) +
  labs(title = "Clonotype counts (Unique CDR3 AA)") +
  theme(
    # 设置标题样式
    plot.title = element_text(hjust = 0.5, size = 18, family = "Arial",face = "bold"),
    # 将图中所有文本元素（包括数字和集合名称）设置为粗体
    text = element_text(family = "Arial", face = "bold")
  )

print(plot6a)
ggsave("./results/Figure6a.pdf",
       plot = plot6a, width = 6, height = 6)
ggsave("./results/Figure6a.png",
       plot = plot6a, width = 6, height = 6, dpi = 300)

##### Figure6b Clone counts------
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(ggplot2)
})

# ---- 1) 计算各区域的“数量”（仍以 CDR3 为键、按行数计数）
A <- "naive"; B <- "flu_H1"; C <- "flu_H5"

canon_cdr3 <- function(x) x |> as.character() |> str_trim() |> str_replace("\\*+$","") |> toupper()

stopifnot(all(c("specificity","cdr3_aa") %in% names(db_igh_fix)))
dat <- db_igh_fix %>%
  mutate(
    specificity   = case_when(
      str_to_lower(specificity) == "naive"  ~ A,
      str_to_lower(specificity) == "flu_h1" ~ B,
      str_to_lower(specificity) == "flu_h5" ~ C,
      TRUE ~ specificity
    ),
    cdr3_aa = canon_cdr3(cdr3_aa)
  ) %>%
  filter(specificity %in% c(A,B,C), !is.na(cdr3_aa), nzchar(cdr3_aa))

present_tbl <- dat %>%
  distinct(specificity, cdr3_aa) %>%
  pivot_wider(names_from = specificity, values_from = specificity, values_fn = length, values_fill = 0) %>%
  mutate(
    !!A := as.integer(.data[[A]] > 0),
    !!B := as.integer(.data[[B]] > 0),
    !!C := as.integer(.data[[C]] > 0),
    pattern = case_when(
      .data[[A]]==1 & .data[[B]]==0 & .data[[C]]==0 ~ "A",
      .data[[A]]==0 & .data[[B]]==1 & .data[[C]]==0 ~ "B",
      .data[[A]]==0 & .data[[B]]==0 & .data[[C]]==1 ~ "C",
      .data[[A]]==1 & .data[[B]]==1 & .data[[C]]==0 ~ "AB",
      .data[[A]]==1 & .data[[B]]==0 & .data[[C]]==1 ~ "AC",
      .data[[A]]==0 & .data[[B]]==1 & .data[[C]]==1 ~ "BC",
      .data[[A]]==1 & .data[[B]]==1 & .data[[C]]==1 ~ "ABC",
      TRUE ~ "OTHER"
    )
  ) %>% select(cdr3_aa, pattern)

dat_with_pattern <- dat %>% left_join(present_tbl, by = "cdr3_aa")

region_counts <- dat_with_pattern %>%
  filter(pattern %in% c("A","B","C","AB","AC","BC","ABC")) %>%
  count(pattern, name = "count") %>%
  complete(pattern = c("A","B","C","AB","AC","BC","ABC"), fill = list(count = 0)) %>%
  arrange(factor(pattern, levels = c("A","B","C","AB","AC","BC","ABC")))

# ---- 2) 画“等半径”三圆韦恩图（纯 ggplot2）
r <- 1.0
# 三个圆的圆心（更居中、留白更足）
pos <- tibble::tibble(
  set = c("A","B","C"),
  x0  = c(-0.55, 1.05, 0.00),
  y0  = c( 0.00, 0.00, -1.15)
)

circle_df <- function(x0, y0, r, n = 361){
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = x0 + r*cos(t), y = y0 + r*sin(t))
}
poly_A <- circle_df(pos$x0[pos$set=="A"], pos$y0[pos$set=="A"], r)
poly_B <- circle_df(pos$x0[pos$set=="B"], pos$y0[pos$set=="B"], r)
poly_C <- circle_df(pos$x0[pos$set=="C"], pos$y0[pos$set=="C"], r)

# --- 各区域数字的位置（更稳健的相对坐标）
cx <- setNames(pos$x0, pos$set)
cy <- setNames(pos$y0, pos$set)

# --- 基准坐标（由三圆圆心推出来）
cx <- setNames(pos$x0, pos$set)
cy <- setNames(pos$y0, pos$set)

base_lab_pos <- tibble::tibble(
  pattern = c("A","B","C","AB","AC","BC","ABC"),
  x = c(
    cx["A"] - 0.35*r,                   # A-only
    cx["B"] + 0.35*r,                   # B-only
    cx["C"],                            # C-only
    (cx["A"]+cx["B"])/2,                # AB
    (cx["A"]+cx["C"])/2 - 0.05,         # AC
    (cx["B"]+cx["C"])/2 + 0.05,         # BC
    (cx["A"]+cx["B"]+cx["C"])/3         # ABC
  ),
  y = c(
    cy["A"],                            # A-only
    cy["B"],                            # B-only
    cy["C"] + 0.55*r,                   # C-only
    (cy["A"]+cy["B"])/2 + 0.12,         # AB
    (cy["A"]+cy["C"])/2 - 0.12,         # AC
    (cy["B"]+cy["C"])/2 - 0.12,         # BC
    (cy["A"]+cy["B"]+cy["C"])/3 - 0.10  # ABC
  )
)

# --- 自定义偏移（想调哪个写哪个；单位与坐标一致）---
# 例子：把 H5 的 6294 稍微上移；把 AB 的 21 往右上挪；把 AC 的 70 往右挪……
offsets <- tibble::tibble(
  pattern = c("C",  "AB",  "AC", "BC", "A", "B", "ABC"),
  dx      = c( 0.00, 0.10, 0.08, 0.00, 0.00, 0.00, 0.10),
  dy      = c( -0.58, 0.04, 0.00, 0.02, 0.00, 0.00, 0.15)
)
# 对没有填写的区域默认偏移 0
offsets <- offsets %>% tidyr::complete(pattern = c("A","B","C","AB","AC","BC","ABC"),
                                       fill = list(dx=0, dy=0))

# --- 合并“基准 + 偏移”并加上计数 ---
lab_pos <- base_lab_pos %>%
  dplyr::left_join(offsets, by = "pattern") %>%
  dplyr::mutate(x = x + dx, y = y + dy) %>%
  dplyr::left_join(region_counts, by = "pattern")

plot6b <- ggplot() +
  geom_polygon(data = poly_A, aes(x, y), fill = "#77AADD", alpha = 0.65, color = "black") +
  geom_polygon(data = poly_B, aes(x, y), fill = "#EFC000", alpha = 0.65, color = "black") +
  geom_polygon(data = poly_C, aes(x, y), fill = "#CD534C", alpha = 0.65, color = "black") +
  # 各区域数字
  geom_text(data = lab_pos, aes(x, y, label = count), size = 6) +
  # 组名（圆外）
  geom_text(data = label_out, aes(x, y, label = lab), fontface = "bold", size = 6) +
  # 画布范围留白，避免出界
  coord_fixed(
    xlim = c(min(pos$x0) - r - 0.2, max(pos$x0) + r + 0.2),
    ylim = c(min(pos$y0) - r - 0.3, max(pos$y0) + r + 0.4),
    clip = "off"
  ) +
  theme_void() +
  ggtitle("Clone counts (Unique CDR3 AA)") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

print(plot6b)

ggsave("./results/Figure6b.pdf", plot6b, width = 6, height = 6)
ggsave("./results/Figure6b.png", plot6b, width = 6, height = 6, dpi = 300)


##### Figure6c clone expansion------
# 必要包
library(dplyr)
library(stringr)
library(forcats)
library(tidyr)
library(ggplot2)
library(packcircles)   # 用于克隆小圆的无重叠打包
library(purrr)
library(ggrepel)
library(scales)
----------
# 交集：Flu_H1 和 Flu_H5 的共有克隆
flu_h1_h5_shared_clones <- intersect(flu_h1_clones, flu_h5_clones)
length(flu_h1_h5_shared_clones)
head(flu_h1_h5_shared_clones)

flu_h1_h5_shared_clones_df <- db_igh_fix %>%
  filter(cdr3_aa %in% flu_h1_h5_shared_clones, specificity %in% c("flu_H1", "flu_H5"))
write.csv(shared_df, "./results/flu_h1_h5_shared_clones.csv")

# flu_h1的特有克隆
# 先算 flu_H1 特有的克隆
unique_flu_h1 <- setdiff(flu_h1_clones, union(flu_h5_clones, naive_clones))
# 查看数量
length(unique_flu_h1)
# 查看前几个
head(unique_flu_h1)
unique_flu_h1_df <- db_igh_fix %>%
  filter(specificity == "flu_H1", cdr3_aa %in% unique_flu_h1)
head(unique_flu_h1_df)
write.csv(unique_flu_h1_df, "./results/unique_flu_h1_df_clones.csv")
#
library(dplyr)
# 1) 先得到 flu_H1 特有的克隆（不在 flu_H5 和 naive 中）
unique_flu_h1 <- setdiff(flu_h1_clones, union(flu_h5_clones, naive_clones))

# 2) 在 flu_H1 里统计每个克隆的扩增大小（细胞数），按降序排，取前10
top10_flu_h1_expansion <- as.data.frame(db_igh_fix) %>%
  filter(specificity == "flu_H1",
         !is.na(cdr3_aa),
         cdr3_aa %in% unique_flu_h1) %>%
  group_by(cdr3_aa) %>%
  summarise(clone_size = n(), .groups = "drop") %>%
  arrange(desc(clone_size)) %>%
  slice_head(n = 10)
print(top10_flu_h1_expansion)


# 3) 提取这些 Top10 克隆在 flu_H1 组内的完整行，形成子集
top10_flu_h1_df <- db_igh_fix %>%
  filter(specificity == "flu_H1",
         cdr3_aa %in% top10_flu_h1_expansion$cdr3_aa)

# 看看子集的前几行
head(top10_flu_h1_df)
unique_flu_h1_df <- db_igh_fix %>%
  filter(cdr3_aa %in% unique_flu_h1)
write.csv(unique_flu_h1_df, "./results/unique_flu_h1.csv")

# flu_h5的特有克隆
unique_flu_h5 <- setdiff(flu_h5_clones, union(flu_h1_clones, naive_clones))
length(unique_flu_h5)
head(unique_flu_h5)
unique_flu_h5 <- db_igh_fix %>%
  filter(specificity == "flu_H5", cdr3_aa %in% unique_flu_h5)
head(unique_flu_h5)
write.csv(unique_flu_h5, "./results/unique_flu_h5_df_clones.csv")
library(dplyr)

# 1) 先得到 flu_H1 特有的克隆（不在 flu_H5 和 naive 中）
unique_flu_h5 <- setdiff(flu_h5_clones, union(flu_h1_clones, naive_clones))

# 2) 在 flu_H1 里统计每个克隆的扩增大小（细胞数），按降序排，取前10
top10_flu_h5_expansion <- as.data.frame(db_igh_fix) %>%
  filter(specificity == "flu_H5",
         !is.na(cdr3_aa),
         cdr3_aa %in% unique_flu_h5) %>%
  group_by(cdr3_aa) %>%
  summarise(clone_size = n(), .groups = "drop") %>%
  arrange(desc(clone_size)) %>%
  slice_head(n = 10)

print(top10_flu_h5_expansion)
unique_flu_h5_df <- db_igh_fix %>%
  filter(cdr3_aa %in% unique_flu_h5)

write.csv(unique_flu_h5_df, "./results/unique_flu_h5.csv")

path_shared <- "./results/flu_h1_h5_shared_clones.csv"
path_h1_unique <- "./results/unique_flu_h1.csv"
path_h5_unique <- "./results/unique_flu_h5.csv"

df_shared <- read.csv(path_shared)    %>% mutate(clone_origin = "Shared")
df_h1     <- read.csv(path_h1_unique) %>% mutate(clone_origin = "Unique Flu_H1")
df_h5     <- read.csv(path_h5_unique) %>% mutate(clone_origin = "Unique Flu_H5")
df_all <- bind_rows(df_shared, df_h1, df_h5) %>% mutate(cell_id_unique = row_number())


df_use <- df_all %>%
  filter(!is.na(clone_id)) %>%
  mutate(
    isotype = str_extract(c_call %||% "", "IGH[AMG]") %>%
      fct_explicit_na(na_level = "Other") %>%
      fct_collapse(IGHA = "IGHA", IGHG = "IGHG", IGHM = "IGHM", Other = "Other"),
    cell_type = as.factor(B_cell_subpopulations),
    shm_bin = case_when(
      is.na(MU_FREQ_HEAVY_TOTAL) ~ "<2%",
      MU_FREQ_HEAVY_TOTAL < 0.02 ~ "<2%",
      MU_FREQ_HEAVY_TOTAL < 0.05 ~ "2–5%",
      TRUE                       ~ ">5%"
    ),
    clone_origin = factor(clone_origin,
                          levels = c("Shared","Unique Flu_H1","Unique Flu_H5")),
    cell_id_unique = if ("cell_id_unique" %in% names(.)) cell_id_unique else row_number()
  ) %>%
  select(clone_id, cell_id_unique, isotype, cell_type, shm_bin, clone_origin, v_call_10x) %>%
  mutate(v_call_10x = as.character(v_call_10x))

stopifnot(nrow(df_use) > 0)


base_r   <- 0.6   # 小圆基准半径（随克隆大小放大）
edge_gap <- 0.92  # 细胞点距离小圆边界的比例

# —— 类型保险（确保字符列）——
df_all <- df_all %>%
  mutate(
    clone_id       = as.character(clone_id),
    cell_id_unique = if ("cell_id_unique" %in% names(.)) as.character(cell_id_unique) else as.character(row_number())
  )

make_layout_one_facet <- function(df_facet, base_r = 0.6, edge_gap = 0.92) {
  df_facet <- as_tibble(df_facet) %>% mutate(clone_id = as.character(clone_id))
  
  sizes <- df_facet %>%
    dplyr::count(clone_id, name = "n") %>%
    mutate(r = base_r * sqrt(pmax(n, 3)))
  
  pack <- packcircles::circleProgressiveLayout(sizes$r)
  centers <- bind_cols(sizes, as_tibble(pack)) %>%
    transmute(clone_id, n, r, cx = x, cy = y)
  
  cells <- df_facet %>%
    left_join(centers, by = "clone_id") %>%
    group_by(clone_id) %>%
    mutate(
      idx   = row_number(),
      theta = 2*pi*(idx-1)/n,
      rad   = r * edge_gap,
      x     = cx + rad * cos(theta),
      y     = cy + rad * sin(theta),
      x0    = cx, y0 = cy
    ) %>% ungroup()
  
  centers_nodes <- centers %>% transmute(clone_id, x0 = cx, y0 = cy, r)
  list(cells = cells, centers = centers_nodes)
}

# 生成布局
split_dfs <- split(df_use, df_use$clone_origin, drop = TRUE)
# 保证三个分面按顺序存在（即使为空）便于后续处理
for (nm in c("Shared","Unique Flu_H1","Unique Flu_H5")) {
  if (!nm %in% names(split_dfs)) split_dfs[[nm]] <- df_use[0,]
}
split_dfs <- split_dfs[c("Shared","Unique Flu_H1","Unique Flu_H5")]

layouts <- imap(split_dfs, ~{
  lay <- make_layout_one_facet(.x, base_r = base_r, edge_gap = edge_gap)
  lay$cells$facet   <- .y
  lay$centers$facet <- .y
  lay
})

cells_xy   <- bind_rows(map(layouts, "cells"))
centers_xy <- bind_rows(map(layouts, "centers"))


if (!all(c("x0","y0") %in% names(cells_xy))) {
  cells_xy <- cells_xy %>%
    left_join(centers_xy %>% select(facet, clone_id, x0, y0), by = c("facet","clone_id"))
}
cells_xy <- cells_xy %>%
  group_by(facet, clone_id) %>%
  mutate(
    x0 = ifelse(is.na(x0), mean(x, na.rm = TRUE), x0),
    y0 = ifelse(is.na(y0), mean(y, na.rm = TRUE), y0)
  ) %>% ungroup()

# 统一 facet 名称，后续所有步骤都用同一套键
cells_xy   <- cells_xy   %>% dplyr::mutate(facet = dplyr::recode(as.character(facet),
                                                                 "Shared" = "Flu_H1_H5_shared"))
centers_xy <- centers_xy %>% dplyr::mutate(facet = dplyr::recode(as.character(facet),
                                                                 "Shared" = "Flu_H1_H5_shared"))


# 3) 二层空间优化 + 组级等比缩放（核心）

# A) 让 H1 / Shared 更松散（团与团之间）
inter_clone_scale <- c(
  "Unique Flu_H5"    = 1.05,
  "Unique Flu_H1"    = 1.65,
  "Flu_H1_H5_shared" = 1.95   # 先 1.80，不要 20。后面可微调
)

# B) 团内松紧（每个克隆内部）
intra_scale <- c(
  "Unique Flu_H5"    = 0.98,
  "Unique Flu_H1"    = 0.92,
  "Flu_H1_H5_shared" = 0.95
)

# —— 计算 r_pack（带“半径下限”），避免 n=1 的 r95=0 不起作用 ——
r_floor_frac <- 0.8
r_min_abs    <- 0.4
pad          <- 1.10

sizes_floor <- centers_xy %>%
  dplyr::transmute(facet = as.character(facet),
                   clone_id, r_floor = pmax(r * r_floor_frac, r_min_abs))

clone_radii <- cells_xy %>%
  dplyr::mutate(facet = as.character(facet)) %>%
  dplyr::group_by(facet, clone_id) %>%
  dplyr::summarise(
    r95 = as.numeric(stats::quantile(sqrt((x - x0)^2 + (y - y0)^2), 0.95, na.rm = TRUE)),
    cx  = mean(x0, na.rm = TRUE),
    cy  = mean(y0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(sizes_floor, by = c("facet","clone_id")) %>%
  dplyr::mutate(
    r_eff  = pmax(coalesce(r95, 0), coalesce(r_floor, r_min_abs)),
    scale_ = coalesce(as.numeric(inter_clone_scale[facet]), 1),   # ★ 防 NA
    r_pack = pmax(r_eff * pad * scale_, 1e-6)
  )

repack_facet <- function(fct_name) {
  rd <- dplyr::filter(clone_radii, facet == fct_name)
  if (!nrow(rd)) return(invisible(NULL))
  pk <- packcircles::circleProgressiveLayout(rd$r_pack)
  offsets <- rd %>%
    dplyr::mutate(dx = pk$x - cx, dy = pk$y - cy) %>%
    dplyr::select(facet, clone_id, dx, dy)
  
  # 平移细胞/中心
  cells_xy  <<- cells_xy %>%
    dplyr::left_join(offsets, by = c("facet","clone_id")) %>%
    dplyr::mutate(
      x  = ifelse(facet == fct_name, x  + coalesce(dx,0), x),
      y  = ifelse(facet == fct_name, y  + coalesce(dy,0), y),
      x0 = ifelse(facet == fct_name, x0 + coalesce(dx,0), x0),
      y0 = ifelse(facet == fct_name, y0 + coalesce(dy,0), y0)
    ) %>% dplyr::select(-dx, -dy)
  
  centers_xy <<- centers_xy %>%
    dplyr::left_join(offsets, by = c("facet","clone_id")) %>%
    dplyr::mutate(
      x0 = ifelse(facet == fct_name, x0 + coalesce(dx,0), x0),
      y0 = ifelse(facet == fct_name, y0 + coalesce(dy,0), y0)
    ) %>% dplyr::select(-dx, -dy)
  
  invisible(NULL)
}
for (fct in names(inter_clone_scale)) repack_facet(fct)

# 团内缩放（一次）
cells_xy <- cells_xy %>%
  dplyr::mutate(facet = as.character(facet)) %>%
  dplyr::group_by(facet, clone_id) %>%
  dplyr::mutate(
    s = coalesce(as.numeric(intra_scale[facet]), 1),
    x = x0 + (x - x0) * s,
    y = y0 + (y - y0) * s
  ) %>% dplyr::ungroup() %>% dplyr::select(-s)

## 3.3 组级直径等比（修正版）：把“想要的松散度”乘进 R_target
# 让 H5 做基准；H1/Shared 额外放大到更松
group_sparsity <- c(
  "Unique Flu_H5"    = 1.00,   # 基准
  "Unique Flu_H1"    = 1.80,   # → 更松
  "Flu_H1_H5_shared" = 10.00    # → 更松
)

centers_stats <- centers_xy %>%
  dplyr::group_by(facet) %>%
  dplyr::summarise(
    n_clones = dplyr::n_distinct(clone_id),
    cx_bar   = mean(x0), cy_bar = mean(y0),
    R_now    = sqrt(max((x0 - cx_bar)^2 + (y0 - cy_bar)^2, na.rm = TRUE)),
    .groups  = "drop"
  ) %>%
  dplyr::filter(is.finite(R_now))

if (nrow(centers_stats)) {
  # 以 H5 当前半径为基准常数（也可以继续用你原来的“最多克隆组”为基准）
  k_const <- (centers_stats %>% dplyr::filter(facet == "Unique Flu_H5"))$R_now /
    sqrt((centers_stats %>% dplyr::filter(facet == "Unique Flu_H5"))$n_clones + 1e-9)
  
  centers_stats <- centers_stats %>%
    dplyr::mutate(
      gscale   = dplyr::coalesce(as.numeric(group_sparsity[facet]), 1),
      R_target = k_const * sqrt(n_clones) * gscale,   # ★ 把“想要更松”乘进来
      s_between = dplyr::if_else(R_now > 0, R_target / R_now, 1)
    )
  
  # 按各自组心各向同性缩放（这一步会保留 3.1 的相对布局，只改变整体“盘子”大小）
  cells_xy <- cells_xy %>%
    dplyr::left_join(centers_stats, by = "facet") %>%
    dplyr::mutate(
      x  = cx_bar + (x  - cx_bar) * s_between,
      y  = cy_bar + (y  - cy_bar) * s_between,
      x0 = cx_bar + (x0 - cx_bar) * s_between,
      y0 = cy_bar + (y0 - cy_bar) * s_between
    ) %>%
    dplyr::select(-cx_bar, -cy_bar, -R_now, -R_target, -s_between, -n_clones, -gscale)
  
  centers_xy <- centers_xy %>%
    dplyr::left_join(centers_stats, by = "facet") %>%
    dplyr::mutate(
      x0 = cx_bar + (x0 - cx_bar) * s_between,
      y0 = cy_bar + (y0 - cy_bar) * s_between
    ) %>%
    dplyr::select(-cx_bar, -cy_bar, -R_now, -R_target, -s_between, -n_clones, -gscale)
}


# 4) 出图：去掉 Origin 图例（描边固定色），保留 Subtype/Isotype/Mutation

pal_origin <- c("Flu_H1_H5_shared" = "#F6C343", "Unique Flu_H1" = "#2C7BE5", "Unique Flu_H5" = "#E63757")
pal_ct     <- c("MZ"="#3A8EBA","GC"="#E69F00","PB"="#D55E00","Bmem"="#009E73")
shape_iso  <- c(IGHA=23, IGHG=22, IGHM=21, Other=25)
size_mut   <- c("<2%"=2.0, "2–5%"=3.0, ">5%"=4.0)

# 统一分面命名
cells_xy   <- cells_xy   %>% dplyr::mutate(facet = dplyr::recode(as.character(facet), "Shared"="Flu_H1_H5_shared"))
centers_xy <- centers_xy %>% dplyr::mutate(facet = dplyr::recode(as.character(facet), "Shared"="Flu_H1_H5_shared"))

# 只保留配色内的亚群，并设定因子水平（保证颜色稳定）
cells_xy <- cells_xy %>%
  dplyr::filter(as.character(cell_type) %in% names(pal_ct)) %>%
  dplyr::mutate(
    facet     = factor(facet, levels = c("Unique Flu_H5","Unique Flu_H1","Flu_H1_H5_shared")),
    cell_type = factor(as.character(cell_type), levels = names(pal_ct)),
    isotype   = factor(as.character(isotype),   levels = names(shape_iso)),
    shm_bin   = factor(as.character(shm_bin),   levels = names(size_mut))
  )

# —— 同一克隆“主亚群”的连接线（中心->细胞），每克隆最多抽样 15 条，避免糊 —— 
set.seed(1)
edges_lines <- cells_xy %>%
  dplyr::group_by(facet, clone_id) %>%
  dplyr::mutate(major_ct = names(sort(table(na.omit(cell_type)), decreasing = TRUE))[1]) %>%
  dplyr::filter(as.character(cell_type) == major_ct) %>%
  dplyr::mutate(.rand = runif(dplyr::n())) %>%
  dplyr::slice_min(order_by = .rand, n = 15, with_ties = FALSE) %>%
  dplyr::ungroup()

# —— 大克隆标签（>20），取每克隆最常见 v_call_10x —— 
# —— 用于打标签的基础聚合：每个克隆的大小 n 与 v_call_10x 众数 —— 
clone_counts <- df_use %>%
  dplyr::mutate(facet = dplyr::recode(as.character(clone_origin),
                                      "Shared" = "Flu_H1_H5_shared")) %>%
  dplyr::group_by(facet, clone_id) %>%
  dplyr::summarise(
    n = dplyr::n(),
    v_top = {
      vt <- stats::na.omit(v_call_10x)
      if (length(vt)) names(sort(table(vt), decreasing = TRUE))[1] else NA_character_
    },
    .groups = "drop"
  )

# 1) H5 组：保留你原逻辑（大克隆 > 20）——标签只显示 V 基因名称
label_min_H5 <- 25
lab_H5 <- clone_counts %>%
  dplyr::filter(facet == "Unique Flu_H5", n > label_min_H5) %>%
  dplyr::left_join(centers_xy %>% dplyr::select(facet, clone_id, x0, y0),
                   by = c("facet","clone_id")) %>%
  dplyr::mutate(label = v_top) %>%
  dplyr::filter(is.finite(x0), is.finite(y0))

## H1：仅显示克隆数 > 20 的克隆
lab_h1 <- clone_counts %>%
  dplyr::filter(facet == "Unique Flu_H1", n > 10) %>%
  dplyr::left_join(centers_xy %>% dplyr::select(facet, clone_id, x0, y0),
                   by = c("facet","clone_id")) %>%
  dplyr::mutate(label = v_top) %>%
  dplyr::filter(is.finite(x0), is.finite(y0))

## Shared：按克隆大小取 Top 4
lab_shared_top4 <- clone_counts %>%
  dplyr::filter(facet == "Flu_H1_H5_shared") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::slice_head(n = 4) %>%      # ← 这里控制 top 数量
  dplyr::left_join(centers_xy %>% dplyr::select(facet, clone_id, x0, y0),
                   by = c("facet","clone_id")) %>%
  dplyr::mutate(label = v_top) %>%
  dplyr::filter(is.finite(x0), is.finite(y0))

## 合并到总标签表
lab_all <- dplyr::bind_rows(lab_H5, lab_h1, lab_shared_top4)


# -------- 单面板绘图函数（无克隆外圈；连主亚群线；右侧统一图例）
draw_panel <- function(facet_name,
                       title = facet_name,
                       title_size = 20,
                       title_bpad = 2,
                       pad_left  = 0.02,   # 左右 padding 比例
                       pad_right = 0.02,
                       pad_top   = 0.005,  # ↑ 上 padding 比例（调小 = 标题更近）
                       pad_bottom= 0.03) { # ↓ 下 padding 比例
  dat_cells <- dplyr::filter(cells_xy, facet == facet_name)
  
  if (nrow(dat_cells) == 0) {
    XLIM <- c(-1, 1); YLIM <- c(-1, 1)
  } else {
    xr <- range(dat_cells$x, na.rm = TRUE)
    yr <- range(dat_cells$y, na.rm = TRUE)
    dx <- diff(xr); dy <- diff(yr)
    if (!is.finite(dx) || dx == 0) dx <- 1
    if (!is.finite(dy) || dy == 0) dy <- 1
    XLIM <- c(xr[1] - dx*pad_left,  xr[2] + dx*pad_right)
    YLIM <- c(yr[1] - dy*pad_bottom, yr[2] + dy*pad_top)  # ← 不对称 padding
  }
  
  ggplot() +
    geom_segment(
      data = dplyr::filter(edges_lines, facet == facet_name),
      aes(x = x0, y = y0, xend = x, yend = y),
      color = scales::alpha("grey55", 0.45), linewidth = 0.22, show.legend = FALSE
    ) +
    geom_point(
      data = dat_cells,
      aes(x = x, y = y, fill = cell_type, color = facet, shape = isotype, size = shm_bin),
      stroke = 0.55, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(lab_all, facet == facet_name),
      aes(x = x0, y = y0, label = label),
      color = "black", fontface = "bold", size = 5,
      min.segment.length = 0, segment.color = "grey60",
      box.padding = 0.25, point.padding = 0.2,
      max.overlaps = Inf, seed = 1, show.legend = FALSE
    ) +
    scale_color_manual(values = pal_origin, name = "Origin (border)", drop = FALSE) +
    scale_fill_manual(values  = pal_ct,     name = "B cell subtype", drop = FALSE) +
    scale_shape_manual(values  = shape_iso, name = "Isotype") +
    scale_size_manual(values   = size_mut,  name = expression(V[H]~"mutation")) +
    coord_fixed(xlim = XLIM, ylim = YLIM, clip = "on") +
    theme_void(base_size = 12) +
    ggtitle(title) +
    theme(
      plot.title = element_text(face = "bold", size = title_size, hjust = 0.5,
                                margin = margin(b = title_bpad)),
      plot.title.position = "panel",
      legend.position = "none",
      plot.margin = margin(t = 2, r = 6, b = 6, l = 6)
    )
}

# 三个面板（左列 H5；右上 H1；右下 Shared）
p_h5 <- draw_panel("Unique Flu_H5",  "Unique Flu_H5",
                   title_size = 20, title_bpad = 1,
                   pad_top = 0.001, pad_bottom = 0.03)

p_h1 <- draw_panel("Unique Flu_H1",  "Unique Flu_H1",
                   title_size = 20, title_bpad = 2,
                   pad_top = 0.012, pad_bottom = 0.03)

p_sh <- draw_panel("Flu_H1_H5_shared","Flu_H1_H5_shared",
                   title_size = 20, title_bpad = 2,
                   pad_top = 0.010, pad_bottom = 0.05)

# -------- 图例（放右侧 & 字体更大）
pal_ct <- c("MZ"="#3A8EBA","GC"="#E69F00","PB"="#D55E00","Bmem"="#009E73")

cells_xy <- cells_xy %>%
  dplyr::mutate(cell_type = as.character(cell_type)) %>%
  dplyr::mutate(
    cell_type = dplyr::case_when(
      cell_type %in% c("MZ","marginal zone","Marginal zone") ~ "MZ",
      cell_type %in% c("GC","germinal center","Germinal center") ~ "GC",
      cell_type %in% c("PB","Plasmablast","plasmablast") ~ "PB",
      cell_type %in% c("Bmem","Memory B cell","B memory","memory B") ~ "Bmem",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = names(pal_ct)))
# --- 构造“独立图例”：每个图例用一层，互不影响 ---
legend_plot <- ggplot() +
  # 1) Origin（描边色）
  geom_point(
    data = dplyr::distinct(cells_xy, facet),
    aes(x = 1, y = 1, color = facet),
    shape = 21, fill = NA, size = 4, stroke = 1.1
  ) +
  # 2) VH mutation（点大小）
  geom_point(
    data = data.frame(shm_bin = factor(names(size_mut), levels = names(size_mut))),
    aes(x = 2, y = 1, size = shm_bin),
    shape = 19, color = "black"
  ) +
  # 3) Isotype（点形状）
  geom_point(
    data = data.frame(isotype = factor(names(shape_iso), levels = names(shape_iso))),
    aes(x = 3, y = 1, shape = isotype),
    size = 4, color = "grey20"
  ) +
  # 4) B cell subtype（填充色，外框固定为灰色，避免受 color 影响）
  geom_point(
    data = data.frame(cell_type = factor(names(pal_ct), levels = names(pal_ct))),
    aes(x = 4, y = 1, fill = cell_type),
    shape = 21, color = "grey30", size = 4, stroke = 0.8
  ) +
  scale_color_manual(values = pal_origin, name = "Origin (border)", drop = FALSE) +
  scale_fill_manual(values  = pal_ct,     name = "B cell subtype",  drop = FALSE) +
  scale_shape_manual(values  = shape_iso,  name = "Isotype") +
  scale_size_manual(values   = size_mut,   name = expression(V[H]~"mutation")) +
  guides(
    color = guide_legend(override.aes = list(fill = NA), order = 1),
    size  = guide_legend(override.aes = list(shape = 19, color = "black"), order = 2),
    shape = guide_legend(order = 3),
    fill  = guide_legend(override.aes = list(shape = 21, colour = "grey30", size = 4), order = 4)
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 13),
    legend.key.size = unit(0.9, "cm"),
    legend.margin = margin(t = 6, r = 8, b = 6, l = 6)
  )

legend_grob <- cowplot::get_legend(legend_plot)

# -------- 三角形排版
right_col  <- cowplot::plot_grid(p_h1, p_sh, ncol = 1, align = "v", rel_heights = c(1.5, 0.5))
tri_layout <- cowplot::plot_grid(p_h5, right_col, ncol = 2, align = "h", rel_widths = c(1.55, 1))

plot6c <- cowplot::plot_grid(tri_layout, legend_grob, ncol = 2, rel_widths = c(1, 0.20))
print(plot6c)

ggsave("./results/Figure6c.pdf", plot = plot6c, width = 20, height = 11)
ggsave("./results/Figure6c.png", plot = plot6c, width = 20, height = 11,dpi = 300)


######Figure6d------
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(scales)

# === 读取数据 ===
flu_H1 <- read.csv("./results/unique_flu_h1_df_clones.csv",
                   check.names = TRUE, stringsAsFactors = FALSE)

# === 清理序列 ===
clean_aa <- function(x){
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}

heavy_col <- "junction_10x_aa"
light_col <- "light_junction_10x_aa"
stopifnot(all(c(heavy_col, light_col) %in% names(flu_H1)))

df_len <- flu_H1 %>%
  transmute(
    Heavy = nchar(clean_aa(.data[[heavy_col]])),
    Light = nchar(clean_aa(.data[[light_col]]))
  ) %>%
  pivot_longer(everything(), names_to = "Chain", values_to = "Length") %>%
  filter(!is.na(Length), Length > 0)

# === 计算均值 ===
means <- df_len %>%
  group_by(Chain) %>%
  summarise(mean_len = mean(Length, na.rm = TRUE), .groups = "drop")

# === 手动设置每条链文字的垂直位置（避免重叠） ===
means <- means %>%
  mutate(
    y_pos = ifelse(Chain == "Heavy", 0.32, 0.26)   # 可根据图中情况微调
  )

# === 颜色 ===
col_map  <- c("Heavy" = "#1f77b4", "Light" = "grey40")
fill_map <- c("Heavy" = "#1f77b4", "Light" = "grey60")

# === 横轴范围自动适配并留空白 ===
x_min <- min(df_len$Length, na.rm = TRUE)
x_max <- max(df_len$Length, na.rm = TRUE)
x_pad <- max(0.1 * (x_max - x_min), 0.5)
x_limits <- c(x_min - x_pad, x_max + x_pad)

# === 绘制频率直方图（相对频率） ===
p <- ggplot(df_len, aes(x = Length, fill = Chain, color = Chain)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 1, position = "identity",
                 alpha = 0.35, linewidth = 0.3) +
  geom_vline(data = means, aes(xintercept = mean_len, color = Chain),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = means,
            aes(x = mean_len, y = y_pos,
                label = paste0("Mean: ", round(mean_len, 1)),
                color = Chain),
            hjust = -0.1, size = 6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = col_map, name = "Chain type",
                     labels = c("Heavy chain", "Light chain")) +
  scale_fill_manual(values  = fill_map, name = "Chain type",
                    labels = c("Heavy chain", "Light chain")) +
  coord_cartesian(xlim = x_limits, ylim = c(0, NA)) +
  labs(x = "CDR3 length (AA)", y = "Proportion",
       title = "flu_H1: CDR3 length histogram") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)


out_dir <- "./results/Figure6"
png_path <- file.path(out_dir, "flu_H1_CDR3_length_histogram_meanAdjusted.png")
pdf_path <- file.path(out_dir, "flu_H1_CDR3_length_histogram_meanAdjusted.pdf")

ggsave(png_path, p, width = 6, height = 5, dpi = 300)
ggsave(pdf_path, p, width = 6, height = 5)

# === cdr3 character ===
library(dplyr)
library(stringr)
library(ggseqlogo)
library(ggplot2)

# === 工具函数 ===
clean_aa <- function(x) {
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}
anchor_two_ends <- function(v, right_anchor = "[FW]") {
  if (length(v) == 0) return(character(0))
  v <- as.character(v)
  lp <- regexpr("C", v, perl = TRUE); lp[lp < 0] <- 1L
  rp <- sapply(v, function(s) {
    hits <- gregexpr(right_anchor, s, perl = TRUE)[[1]]
    if (all(hits < 0)) nchar(s) else max(hits)
  })
  left_part  <- mapply(substr, v, 1, lp, USE.NAMES = FALSE)
  mid_part   <- mapply(function(s, l, r)
    if (r - l - 1 <= 0) "" else substr(s, l + 1, r - 1),
    v, lp, rp, USE.NAMES = FALSE)
  right_part <- mapply(substr, v, rp, nchar(v), USE.NAMES = FALSE)
  max_mid <- max(nchar(mid_part))
  mapply(function(lf, mid, rt)
    paste0(lf, mid, strrep("-", max_mid - nchar(mid)), rt),
    left_part, mid_part, right_part, USE.NAMES = FALSE)
}
trim_sparse_cols <- function(seqs, min_non_gap_prop = 0.1) {
  if (length(seqs) == 0) return(seqs)
  L <- max(nchar(seqs))
  mat <- sapply(seqs, function(s) strsplit(sprintf("%-*s", L, s), "")[[1]])
  mat[mat == " "] <- "-"
  keep <- apply(mat, 1, function(col) mean(col != "-")) >= min_non_gap_prop
  if (!any(keep)) return(seqs)
  seqs_trim <- apply(mat[keep, , drop = FALSE], 2, paste0, collapse = "")
  unname(seqs_trim)
}
pad_right <- function(v) {
  if (length(v) == 0) return(character(0))
  stringr::str_pad(v, max(nchar(v)), side = "right", pad = "-")
}

stopifnot(all(c("junction_10x_aa", "clone_id") %in% names(flu_H1)))
has_light <- "light_junction_10x_aa" %in% names(flu_H1)

flu_H1_clean <- flu_H1 %>%
  transmute(
    clone_id = clone_id,
    cdr3_h   = clean_aa(junction_10x_aa),
    len_h    = nchar(cdr3_h),
    cdr3_l   = if (has_light) clean_aa(light_junction_10x_aa) else NA_character_,
    len_l    = if (has_light) nchar(cdr3_l) else NA_integer_
  ) %>%
  filter(!is.na(cdr3_h), len_h > 0)


length_stats <- flu_H1_clean %>%
  group_by(len_h) %>%
  summarise(n_clones = n_distinct(clone_id), n_seq = n(), .groups = "drop") %>%
  arrange(len_h)
print(length_stats)

clone_per_len <- flu_H1_clean %>%
  group_by(len_h, clone_id) %>%
  summarise(n_seq = n(), .groups = "drop") %>%
  arrange(len_h, desc(n_seq))
print(clone_per_len)

target_len_h <- 20
out_dir <- "./results/Figure6"

seqs_h <- flu_H1_clean %>%
  filter(len_h == target_len_h) %>%
  pull(cdr3_h) %>%
  unique()

aligned_h <- seqs_h %>%
  anchor_two_ends("[FW]") %>%
  trim_sparse_cols(0.1) %>%
  pad_right()

if (length(aligned_h) > 0) {
  p_heavy <- ggseqlogo(aligned_h, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
    labs(title = paste0("flu_H1 - Heavy CDR3 (length ", target_len_h, ")"),
         x = "Position", y = "Probability") +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
      axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
      axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
      axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
      axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
      legend.title = element_text(size = 16, face = "bold"),             # 图例标题
      legend.text = element_text(size = 14),                             # 图例文字
      strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
    )
  print(p_heavy)
  
  output_path_png <- file.path(out_dir, paste0("flu_H1_HeavyCDR3_len", target_len_h, ".png"))
  output_path_pdf <- file.path(out_dir, paste0("flu_H1_HeavyCDR3_len", target_len_h, ".pdf"))
  ggsave(output_path_png, p_heavy, width = 11, height = 4.5, dpi = 300)
  ggsave(output_path_pdf, p_heavy, width = 11, height = 4.5)
}

if (has_light) {
  seqs_l_all <- flu_H1_clean %>%
    filter(len_h == target_len_h, !is.na(cdr3_l), len_l > 0) %>%
    pull(cdr3_l) %>%
    unique()
  
  aligned_l_all <- seqs_l_all %>%
    anchor_two_ends("[FW]") %>%
    trim_sparse_cols(0.1) %>%
    pad_right()
  
  if (length(aligned_l_all) > 0) {
    p_light_all <- ggseqlogo(aligned_l_all, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
      labs(title = paste0("flu_H1 - Light CDR3 (all lengths; heavy=", target_len_h, ")"),
           x = "Position", y = "Probability") +
      theme_bw()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
        axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
        axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
        axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
        axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
        legend.title = element_text(size = 16, face = "bold"),             # 图例标题
        legend.text = element_text(size = 14),                             # 图例文字
        strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
      )
    print(p_light_all)
    
    output_path_png <- file.path(out_dir, paste0("flu_H1_LightCDR3_allLen_withHeavy", target_len_h, ".png"))
    output_path_pdf <- file.path(out_dir, paste0("flu_H1_LightCDR3_allLen_withHeavy", target_len_h, ".pdf"))
    ggsave(output_path_png, p_light_all, width = 11, height = 4.5, dpi = 300)
    ggsave(output_path_pdf, p_light_all, width = 11, height = 4.5)
  }
}
######Figure6e------
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(scales)

# === 读取数据 ===
flu_H5 <- read.csv("./results/unique_flu_H5_df_clones.csv",
                   check.names = TRUE, stringsAsFactors = FALSE)

# === 清理序列 ===
clean_aa <- function(x){
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}

heavy_col <- "junction_10x_aa"
light_col <- "light_junction_10x_aa"
stopifnot(all(c(heavy_col, light_col) %in% names(flu_H5)))

df_len <- flu_H5 %>%
  transmute(
    Heavy = nchar(clean_aa(.data[[heavy_col]])),
    Light = nchar(clean_aa(.data[[light_col]]))
  ) %>%
  pivot_longer(everything(), names_to = "Chain", values_to = "Length") %>%
  filter(!is.na(Length), Length > 0)

# === 计算均值 ===
means <- df_len %>%
  group_by(Chain) %>%
  summarise(mean_len = mean(Length, na.rm = TRUE), .groups = "drop")

# === 手动设置每条链文字的垂直位置（避免重叠） ===
means <- means %>%
  mutate(
    y_pos = ifelse(Chain == "Heavy", 0.32, 0.26)   # 可根据图中情况微调
  )

# === 颜色 ===
col_map  <- c("Heavy" = "#1f77b4", "Light" = "grey40")
fill_map <- c("Heavy" = "#1f77b4", "Light" = "grey60")

# === 横轴范围自动适配并留空白 ===
x_min <- min(df_len$Length, na.rm = TRUE)
x_max <- max(df_len$Length, na.rm = TRUE)
x_pad <- max(0.1 * (x_max - x_min), 0.5)
x_limits <- c(x_min - x_pad, x_max + x_pad)

# === 绘制频率直方图（相对频率） ===
p <- ggplot(df_len, aes(x = Length, fill = Chain, color = Chain)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 1, position = "identity",
                 alpha = 0.35, linewidth = 0.3) +
  geom_vline(data = means, aes(xintercept = mean_len, color = Chain),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = means,
            aes(x = mean_len, y = y_pos,
                label = paste0("Mean: ", round(mean_len, 1)),
                color = Chain),
            hjust = -0.1, size = 6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = col_map, name = "Chain type",
                     labels = c("Heavy chain", "Light chain")) +
  scale_fill_manual(values  = fill_map, name = "Chain type",
                    labels = c("Heavy chain", "Light chain")) +
  coord_cartesian(xlim = x_limits, ylim = c(0, NA)) +
  labs(x = "CDR3 length (AA)", y = "Proportion",
       title = "flu_H5: CDR3 length histogram") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

# === 保存输出 ===
out_dir <- "./results/Figure6"
png_path <- file.path(out_dir, "flu_H5_CDR3_length_histogram_meanAdjusted.png")
pdf_path <- file.path(out_dir, "flu_H5_CDR3_length_histogram_meanAdjusted.pdf")

ggsave(png_path, p, width = 6, height = 5, dpi = 300)
ggsave(pdf_path, p, width = 6, height = 5)


library(dplyr)
library(stringr)
library(ggseqlogo)
library(ggplot2)

# === 工具函数 ===
clean_aa <- function(x) {
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}
anchor_two_ends <- function(v, right_anchor = "[FW]") {
  if (length(v) == 0) return(character(0))
  v <- as.character(v)
  lp <- regexpr("C", v, perl = TRUE); lp[lp < 0] <- 1L
  rp <- sapply(v, function(s) {
    hits <- gregexpr(right_anchor, s, perl = TRUE)[[1]]
    if (all(hits < 0)) nchar(s) else max(hits)
  })
  left_part  <- mapply(substr, v, 1, lp, USE.NAMES = FALSE)
  mid_part   <- mapply(function(s, l, r)
    if (r - l - 1 <= 0) "" else substr(s, l + 1, r - 1),
    v, lp, rp, USE.NAMES = FALSE)
  right_part <- mapply(substr, v, rp, nchar(v), USE.NAMES = FALSE)
  max_mid <- max(nchar(mid_part))
  mapply(function(lf, mid, rt)
    paste0(lf, mid, strrep("-", max_mid - nchar(mid)), rt),
    left_part, mid_part, right_part, USE.NAMES = FALSE)
}
trim_sparse_cols <- function(seqs, min_non_gap_prop = 0.1) {
  if (length(seqs) == 0) return(seqs)
  L <- max(nchar(seqs))
  mat <- sapply(seqs, function(s) strsplit(sprintf("%-*s", L, s), "")[[1]])
  mat[mat == " "] <- "-"
  keep <- apply(mat, 1, function(col) mean(col != "-")) >= min_non_gap_prop
  if (!any(keep)) return(seqs)
  seqs_trim <- apply(mat[keep, , drop = FALSE], 2, paste0, collapse = "")
  unname(seqs_trim)
}
pad_right <- function(v) {
  if (length(v) == 0) return(character(0))
  stringr::str_pad(v, max(nchar(v)), side = "right", pad = "-")
}

# === 预处理 ===
stopifnot(all(c("junction_10x_aa", "clone_id") %in% names(flu_H5)))
has_light <- "light_junction_10x_aa" %in% names(flu_H5)

flu_H5_clean <- flu_H5 %>%
  transmute(
    clone_id = clone_id,
    cdr3_h   = clean_aa(junction_10x_aa),
    len_h    = nchar(cdr3_h),
    cdr3_l   = if (has_light) clean_aa(light_junction_10x_aa) else NA_character_,
    len_l    = if (has_light) nchar(cdr3_l) else NA_integer_
  ) %>%
  filter(!is.na(cdr3_h), len_h > 0)

# === 重链长度分布 ===
length_stats <- flu_H5_clean %>%
  group_by(len_h) %>%
  summarise(n_clones = n_distinct(clone_id), n_seq = n(), .groups = "drop") %>%
  arrange(len_h)
print(length_stats)

clone_per_len <- flu_H5_clean %>%
  group_by(len_h, clone_id) %>%
  summarise(n_seq = n(), .groups = "drop") %>%
  arrange(len_h, desc(n_seq))
print(clone_per_len)

# === 设定要分析的重链长度 ===
target_len_h <- 20
out_dir <- "./results/Figure6"

# ========== 1) 重链
seqs_h <- flu_H5_clean %>%
  filter(len_h == target_len_h) %>%
  pull(cdr3_h) %>%
  unique()

aligned_h <- seqs_h %>%
  anchor_two_ends("[FW]") %>%
  trim_sparse_cols(0.1) %>%
  pad_right()

if (length(aligned_h) > 0) {
  p_heavy <- ggseqlogo(aligned_h, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
    labs(title = paste0("flu_H5 - Heavy CDR3 (length ", target_len_h, ")"),
         x = "Position", y = "Probability") +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
      axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
      axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
      axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
      axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
      legend.title = element_text(size = 16, face = "bold"),             # 图例标题
      legend.text = element_text(size = 14),                             # 图例文字
      strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
    )
  print(p_heavy)
  
  output_path_png <- file.path(out_dir, paste0("flu_H5_HeavyCDR3_len", target_len_h, ".png"))
  output_path_pdf <- file.path(out_dir, paste0("flu_H5_HeavyCDR3_len", target_len_h, ".pdf"))
  ggsave(output_path_png, p_heavy, width = 11, height = 4.5, dpi = 300)
  ggsave(output_path_pdf, p_heavy, width = 11, height = 4.5)
}

# ========== 2) 轻链总体
if (has_light) {
  seqs_l_all <- flu_H5_clean %>%
    filter(len_h == target_len_h, !is.na(cdr3_l), len_l > 0) %>%
    pull(cdr3_l) %>%
    unique()
  
  aligned_l_all <- seqs_l_all %>%
    anchor_two_ends("[FW]") %>%
    trim_sparse_cols(0.1) %>%
    pad_right()
  
  if (length(aligned_l_all) > 0) {
    p_light_all <- ggseqlogo(aligned_l_all, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
      labs(title = paste0("flu_H5 - Light CDR3 (all lengths; heavy=", target_len_h, ")"),
           x = "Position", y = "Probability") +
      theme_bw()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
        axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
        axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
        axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
        axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
        legend.title = element_text(size = 16, face = "bold"),             # 图例标题
        legend.text = element_text(size = 14),                             # 图例文字
        strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
      )
    print(p_light_all)
    
    output_path_png <- file.path(out_dir, paste0("flu_H5_LightCDR3_allLen_withHeavy", target_len_h, ".png"))
    output_path_pdf <- file.path(out_dir, paste0("flu_H5_LightCDR3_allLen_withHeavy", target_len_h, ".pdf"))
    ggsave(output_path_png, p_light_all, width = 11, height = 4.5, dpi = 300)
    ggsave(output_path_pdf, p_light_all, width = 11, height = 4.5)
  }
}
######Figure6f------
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(scales)

# === 读取数据 ===
shared <- read.csv("./results/flu_h1_h5_shared_clones.csv",
                   check.names = TRUE, stringsAsFactors = FALSE)

# === 清理序列 ===
clean_aa <- function(x){
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}

heavy_col <- "junction_10x_aa"
light_col <- "light_junction_10x_aa"
stopifnot(all(c(heavy_col, light_col) %in% names(shared)))

df_len <- shared %>%
  transmute(
    Heavy = nchar(clean_aa(.data[[heavy_col]])),
    Light = nchar(clean_aa(.data[[light_col]]))
  ) %>%
  pivot_longer(everything(), names_to = "Chain", values_to = "Length") %>%
  filter(!is.na(Length), Length > 0)

# === 计算均值 ===
means <- df_len %>%
  group_by(Chain) %>%
  summarise(mean_len = mean(Length, na.rm = TRUE), .groups = "drop")

# === 手动设置每条链文字的垂直位置（避免重叠） ===
means <- means %>%
  mutate(
    y_pos = ifelse(Chain == "Heavy", 0.32, 0.26)   # 可根据图中情况微调
  )

# === 颜色 ===
col_map  <- c("Heavy" = "#1f77b4", "Light" = "grey40")
fill_map <- c("Heavy" = "#1f77b4", "Light" = "grey60")

# === 横轴范围自动适配并留空白 ===
x_min <- min(df_len$Length, na.rm = TRUE)
x_max <- max(df_len$Length, na.rm = TRUE)
x_pad <- max(0.1 * (x_max - x_min), 0.5)
x_limits <- c(x_min - x_pad, x_max + x_pad)

# === 绘制频率直方图（相对频率） ===
p <- ggplot(df_len, aes(x = Length, fill = Chain, color = Chain)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 1, position = "identity",
                 alpha = 0.35, linewidth = 0.3) +
  geom_vline(data = means, aes(xintercept = mean_len, color = Chain),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = means,
            aes(x = mean_len, y = y_pos,
                label = paste0("Mean: ", round(mean_len, 1)),
                color = Chain),
            hjust = -0.1, size = 6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = col_map, name = "Chain type",
                     labels = c("Heavy chain", "Light chain")) +
  scale_fill_manual(values  = fill_map, name = "Chain type",
                    labels = c("Heavy chain", "Light chain")) +
  coord_cartesian(xlim = x_limits, ylim = c(0, NA)) +
  labs(x = "CDR3 length (AA)", y = "Proportion",
       title = "shared: CDR3 length histogram") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

# === 保存输出 ===
out_dir <- "./results/Figure6"
png_path <- file.path(out_dir, "shared_CDR3_length_histogram_meanAdjusted.png")
pdf_path <- file.path(out_dir, "shared_CDR3_length_histogram_meanAdjusted.pdf")

ggsave(png_path, p, width = 6, height = 5, dpi = 300)
ggsave(pdf_path, p, width = 6, height = 5)


library(dplyr)
library(stringr)
library(ggseqlogo)
library(ggplot2)

# === 工具函数 ===
clean_aa <- function(x) {
  x <- toupper(x)
  x <- str_replace_all(x, "\\*", "")
  x <- str_replace_all(x, "[^ACDEFGHIKLMNPQRSTVWY-]", "")
  x
}
anchor_two_ends <- function(v, right_anchor = "[FW]") {
  if (length(v) == 0) return(character(0))
  v <- as.character(v)
  lp <- regexpr("C", v, perl = TRUE); lp[lp < 0] <- 1L
  rp <- sapply(v, function(s) {
    hits <- gregexpr(right_anchor, s, perl = TRUE)[[1]]
    if (all(hits < 0)) nchar(s) else max(hits)
  })
  left_part  <- mapply(substr, v, 1, lp, USE.NAMES = FALSE)
  mid_part   <- mapply(function(s, l, r)
    if (r - l - 1 <= 0) "" else substr(s, l + 1, r - 1),
    v, lp, rp, USE.NAMES = FALSE)
  right_part <- mapply(substr, v, rp, nchar(v), USE.NAMES = FALSE)
  max_mid <- max(nchar(mid_part))
  mapply(function(lf, mid, rt)
    paste0(lf, mid, strrep("-", max_mid - nchar(mid)), rt),
    left_part, mid_part, right_part, USE.NAMES = FALSE)
}
trim_sparse_cols <- function(seqs, min_non_gap_prop = 0.1) {
  if (length(seqs) == 0) return(seqs)
  L <- max(nchar(seqs))
  mat <- sapply(seqs, function(s) strsplit(sprintf("%-*s", L, s), "")[[1]])
  mat[mat == " "] <- "-"
  keep <- apply(mat, 1, function(col) mean(col != "-")) >= min_non_gap_prop
  if (!any(keep)) return(seqs)
  seqs_trim <- apply(mat[keep, , drop = FALSE], 2, paste0, collapse = "")
  unname(seqs_trim)
}
pad_right <- function(v) {
  if (length(v) == 0) return(character(0))
  stringr::str_pad(v, max(nchar(v)), side = "right", pad = "-")
}

# === 预处理 ===
stopifnot(all(c("junction_10x_aa", "clone_id") %in% names(shared)))
has_light <- "light_junction_10x_aa" %in% names(shared)

shared_clean <- shared %>%
  transmute(
    clone_id = clone_id,
    cdr3_h   = clean_aa(junction_10x_aa),
    len_h    = nchar(cdr3_h),
    cdr3_l   = if (has_light) clean_aa(light_junction_10x_aa) else NA_character_,
    len_l    = if (has_light) nchar(cdr3_l) else NA_integer_
  ) %>%
  filter(!is.na(cdr3_h), len_h > 0)

# === 重链长度分布 ===
length_stats <- shared_clean %>%
  group_by(len_h) %>%
  summarise(n_clones = n_distinct(clone_id), n_seq = n(), .groups = "drop") %>%
  arrange(len_h)
print(length_stats)

clone_per_len <- shared_clean %>%
  group_by(len_h, clone_id) %>%
  summarise(n_seq = n(), .groups = "drop") %>%
  arrange(len_h, desc(n_seq))
print(clone_per_len)

# === 设定要分析的重链长度 ===
target_len_h <- 16
out_dir <- "./results/Figure6"

# ========== 1) 重链
seqs_h <- shared_clean %>%
  filter(len_h == target_len_h) %>%
  pull(cdr3_h) %>%
  unique()

aligned_h <- seqs_h %>%
  anchor_two_ends("[FW]") %>%
  trim_sparse_cols(0.1) %>%
  pad_right()

if (length(aligned_h) > 0) {
  p_heavy <- ggseqlogo(aligned_h, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
    labs(title = paste0("shared - Heavy CDR3 (length ", target_len_h, ")"),
         x = "Position", y = "Probability") +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
      axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
      axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
      axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
      axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
      legend.title = element_text(size = 16, face = "bold"),             # 图例标题
      legend.text = element_text(size = 14),                             # 图例文字
      strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
    )
  print(p_heavy)
  
  output_path_png <- file.path(out_dir, paste0("shared_HeavyCDR3_len", target_len_h, ".png"))
  output_path_pdf <- file.path(out_dir, paste0("shared_HeavyCDR3_len", target_len_h, ".pdf"))
  ggsave(output_path_png, p_heavy, width = 11, height = 4.5, dpi = 300)
  ggsave(output_path_pdf, p_heavy, width = 11, height = 4.5)
}

# ========== 2) 轻链总体
if (has_light) {
  seqs_l_all <- shared_clean %>%
    filter(len_h == target_len_h, !is.na(cdr3_l), len_l > 0) %>%
    pull(cdr3_l) %>%
    unique()
  
  aligned_l_all <- seqs_l_all %>%
    anchor_two_ends("[FW]") %>%
    trim_sparse_cols(0.1) %>%
    pad_right()
  
  if (length(aligned_l_all) > 0) {
    p_light_all <- ggseqlogo(aligned_l_all, seq_type = "aa", method = "prob", col_scheme = "chemistry") +
      labs(title = paste0("shared - Light CDR3 (all lengths; heavy=", target_len_h, ")"),
           x = "Position", y = "Probability") +
      theme_bw()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # 主标题字体更大
        axis.title.x = element_text(size = 18, face = "bold"),             # X轴标题
        axis.title.y = element_text(size = 18, face = "bold"),             # Y轴标题
        axis.text.x = element_text(size = 14, face = "bold"),              # X轴刻度
        axis.text.y = element_text(size = 14, face = "bold"),              # Y轴刻度
        legend.title = element_text(size = 16, face = "bold"),             # 图例标题
        legend.text = element_text(size = 14),                             # 图例文字
        strip.text = element_text(size = 16, face = "bold")                # 分面标题（若存在）
      )
    print(p_light_all)
    
    output_path_png <- file.path(out_dir, paste0("shared_LightCDR3_allLen_withHeavy", target_len_h, ".png"))
    output_path_pdf <- file.path(out_dir, paste0("shared_LightCDR3_allLen_withHeavy", target_len_h, ".pdf"))
    ggsave(output_path_png, p_light_all, width = 11, height = 4.5, dpi = 300)
    ggsave(output_path_pdf, p_light_all, width = 11, height = 4.5)
  }
}



save.image(file = "./results/flu_mouse.project.RData")
