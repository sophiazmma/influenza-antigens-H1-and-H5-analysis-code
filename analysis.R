
# 01. loading library ------------------------------------------------------------
# Define vector of package names to load
pkgs <- c(
  "stringr", "survival", "glmnet", "survminer", "data.table",
  "ggpubr", "dplyr", "patchwork", "matrixStats", "readr", "tibble", "ggplot2",
  "future", "pheatmap", "msigdbr", 
  "Seurat",  "harmony", 
   "RColorBrewer"
)

# Check if packages are installed, install if missing, then load
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


library(shazam)#Package for SHM and clonal analysis
library(ggplot2)
data(ExampleDb, package="alakazam")
native_mouse2_heavy_parse <- read.table("./cell_reports/VDJ/annotation_results/heavy_native_mouse2_parse-select.tsv",sep = "\t",header = T,stringsAsFactors = F)
native_mouse2_dist_ham <- distToNearest(native_mouse2_heavy_parse, sequenceColumn="junction", #Calculate distance to nearest neighbor (Hamming distance)
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=1)

native_mouse2_output_ham <- findThreshold(native_mouse2_dist_ham$dist_nearest, method="density")
native_mouse2_output_ham@threshold


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
  
  # Use the full path provided
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
  file_path = file_paths,       # file_path
  sample_name = names(file_paths) # sample_name
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
  
  # 2.5. Merge heavy and light chain data
  merged_df <- left_join(heavy_df, light_info_to_merge, by = "cell_id")
  
  # 2.6. Construct output filename and save
  #    : [prefix].merge.bcr.tsv
  output_filename <- paste0(prefix, ".merge.bcr.tsv")
  output_filepath <- file.path(output_directory, output_filename)
  
  write_tsv(merged_df, output_filepath)
  
  print(paste("Successfully merged and saved to:", output_filepath))
}

print("--- All samples processed and saved. ---")

list.files(path = output_directory, pattern = "\\.merge\\.bcr\\.tsv$")

###shm calculation for Heavy Chain------

# --- Step 1: Set working directory ---
#  .merge.final.bcr.tsv 
work_dir <- "./cell_reports/VDJ/annotation_results"

# --- Step 2: Find all files to process ---
cat("Finding files to process...\n")

# Find all files ending with ".merge.final.bcr.tsv"
bcr_files <- list.files(path = work_dir, 
                        pattern = "\\.merge\\.final\\.bcr\\.tsv$", 
                        full.names = TRUE)

# Check if files were found
if (length(bcr_files) == 0) {
  stop("No '*.merge.final.bcr.tsv' files found in the specified directory. Please check your 'work_dir' path.")
} else {
  cat("Found the following files to process:\n")
  print(basename(bcr_files))
}

# --- Step 3: Process each file using for loop ---
cat("\nStarting mutation quantification for each file...\n")

for (file_path in bcr_files) {
  
  cat('--------------------------------------------------\n')
  cat('Processing file:', basename(file_path), '\n')
  
  # 1. Read data
  db <- tryCatch({
    read_tsv(file_path, col_types = cols(.default = "c"))
  }, error = function(e) {
    warning(paste("Could not read file:", basename(file_path), ". Skipping. Error:", e$message))
    return(NULL)
  })
  
  if (is.null(db)) {
    next
  }
  
  # 2. Check if required columns exist
  required_cols <- c("sequence_alignment", "germline_alignment_d_mask", "locus")
  if (!all(required_cols %in% colnames(db))) {
    warning(paste("File", basename(file_path), "is missing required columns. Skipping."))
    next
  }
  
  # 3. Filter for heavy chains
  heavy_chains_db <- db %>% filter(locus == "IGH")
  
  if (nrow(heavy_chains_db) > 0) {
    count_db <- observedMutations(heavy_chains_db, 
                                  sequenceColumn = "sequence_alignment", 
                                  germlineColumn = "germline_alignment_d_mask",
                                  regionDefinition=NULL,
                                  frequency = FALSE)
    
    freq_db <- observedMutations(heavy_chains_db, 
                                 sequenceColumn = "sequence_alignment", 
                                 germlineColumn = "germline_alignment_d_mask",
                                 regionDefinition=NULL,
                                 frequency = TRUE,                                 combine = TRUE)
    # 4. Prepare mutation info for merging
    mutation_counts <- count_db %>% 
      select(sequence_id, MU_COUNT_HEAVY_seq_r = mu_count_seq_r,MU_COUNT_HEAVY_seq_s = mu_count_seq_s)
    
    mutation_freqs <- freq_db %>%
      select(sequence_id, MU_FREQ_HEAVY_TOTAL = mu_freq)
    
    # 5. Merge mutation counts and frequencies back to original dataframe
    db_with_mutations <- db %>%
      left_join(mutation_counts, by = "sequence_id") %>%
      left_join(mutation_freqs, by = "sequence_id")
    
  } else {
    warning(paste("No heavy chains (locus == 'IGH') found in file:", basename(file_path)))
    db_with_mutations <- db %>%
      mutate(
        MU_COUNT_HEAVY_TOTAL = NA_integer_, 
        MU_FREQ_HEAVY_TOTAL = NA_real_
      )
  }
  
  # 6. Construct new output filename and save
  # Replace ".tsv" with ".shm.tsv"
  output_filepath <- gsub("\\.tsv$", ".shm.tsv", file_path)
  
  write_tsv(db_with_mutations, output_filepath)
  
  cat('Successfully calculated mutations and saved to:', basename(output_filepath), '\n')
}

cat('--------------------------------------------------\n')
cat('All samples processed and saved successfully!\n')

###shm calculation for Light chain------

# --- Step 0: Install and load required R packages ---

library(dplyr)
library(readr)
library(shazam)

# --- Step 1: Load, merge and prepare data ---

# 1.1 Define file paths and names
base_path <- "./cell_reports/VDJ/annotation_results"file_paths <- c(
  fluH1_mouse1 = file.path(base_path, "light_fluH1_mouse1_parse-select.tsv"),
  fluH1_mouse2 = file.path(base_path, "light_fluH1_mouse2_parse-select.tsv"),
  fluH5_mouse = file.path(base_path, "light_fluH5_mouse_parse-select.tsv"),
  naive_mouse1 = file.path(base_path, "light_naive_mouse1_parse-select.tsv"),
  naive_mouse2 = file.path(base_path, "light_naive_mouse2_parse-select.tsv")
)

# 1.2 Read files one by one
list_of_dfs <- lapply(file_paths, function(path) {
  if (file.exists(path)) {
    read_tsv(path, col_types = cols(.default = "c"))
  } else {
    warning("File not found: ", path, ". Skipped.")
    NULL
  }
})

list_of_dfs <- list_of_dfs[!sapply(list_of_dfs, is.null)]

if (length(list_of_dfs) == 0) {
  stop("Error: All specified files failed to read. Please check file names.")
}

# 1.3 Create 2 merged files
list_of_dfs$fluH1_mouse_combined <- bind_rows(list_of_dfs$fluH1_mouse1, list_of_dfs$fluH1_mouse2)
list_of_dfs$naive_mouse_combined <- bind_rows(list_of_dfs$naive_mouse1, list_of_dfs$naive_mouse2)

message("All 7 data samples loaded and ready.")
print(names(list_of_dfs))

# --- Step 2: Set output directory ---
results_directory <- "light_shm_results"
if (!dir.exists(results_directory)) {
  dir.create(results_directory)
}

# --- Step 3: Define function to calculate SHM for single sample ---
calculate_shm_for_light_chain <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- SHM: ", sample_name, " ---")
  
  # Check required columns
  required_cols <- c("sequence_id", "sequence_alignment", "germline_alignment")
  if (!all(required_cols %in% colnames(sample_data))) {
    warning("Sample '", sample_name, "' ")
    return(NULL)
  }
  
  db <- sample_data
  
  if (nrow(db) > 0) {
    message("Found ", nrow(db), " light chain sequences for calculation")
    
    
    message("  Calculating mutation counts (COUNT)...")
    count_db <- observedMutations(db, 
                                  sequenceColumn = "sequence_alignment", 
                                  germlineColumn = "germline_alignment",
                                  regionDefinition = NULL,
                                  frequency = FALSE,
                                  nproc = 1)
    
    message("  Calculating mutation frequencies (FREQUENCY)...")
    freq_db <- observedMutations(db, 
                                 sequenceColumn = "sequence_alignment", 
                                 germlineColumn = "germline_alignment",
                                 regionDefinition = NULL,
                                 frequency = TRUE,
                                 nproc = 1)
    
    # Organize count and frequency results
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
    
    # Merge count and frequency info back to original dataframe
    db_with_mutations <- db %>%
      left_join(mutation_counts, by = "sequence_id") %>%
      left_join(mutation_freqs, by = "sequence_id")
    
  } else {
    warning("Sample '", sample_name, "' Found")
    db_with_mutations <- db %>%
      mutate(
        MU_COUNT_LIGHT_TOTAL = NA_integer_, 
        MU_FREQ_LIGHT_TOTAL = NA_real_
      )
  }
  
  # Construct new output filename and save
  output_filename <- paste0(sample_name, "_light_with_shm.tsv")
  output_filepath <- file.path(output_dir, output_filename)
  
  write_tsv(db_with_mutations, output_filepath)
  
  message("Successfully calculated light chain SHM and saved to: ", output_filename)
}

# --- Step 4: Loop through all data frames ---
for (name in names(list_of_dfs)) {
  current_data <- list_of_dfs[[name]]
  tryCatch({
    calculate_shm_for_light_chain(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

message("\n\nAll samples processed! Please check '", results_directory, "' folder")

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
  
  # --- Step 1: Initial global QC ---
  clean_bcr <- raw_bcr_df %>%
    filter(toupper(as.character(productive)) == "T", 
           toupper(as.character(high_confidence)) == "TRUE",
           !is.na(umis))
  
  # If no data left after initial filtering, stop and warn
  if(nrow(clean_bcr) == 0) {
    warning("Initial filtering (productive/high_confidence)  flu_H5.bcr ")
    return(tibble()) # tibble
  }
  
  heavy_chains_clean <- clean_bcr %>%
    filter(chain == "IGH") %>%
    group_by(barcode) %>%
    slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  light_chains_clean <- clean_bcr %>%
    filter(chain %in% c("IGK", "IGL")) %>%
    group_by(barcode) %>%
    slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  final_bcr_table <- full_join(heavy_chains_clean, light_chains_clean, by = "barcode")
  
  return(final_bcr_table)
}

add_bcr_to_seurat <- function(seurat_obj, cleaned_bcr_df) {
  
  # 1. Get current Seurat metadata
  seurat_metadata <- seurat_obj@meta.data
  
  seurat_metadata$barcode_seurat <- rownames(seurat_metadata)
  
  # 3. Process BCR data: remove "_contig_X" suffix from sequence_id and rename to barcode
  cleaned_bcr_df$barcode <- sub("_contig_\\d+", "", cleaned_bcr_df$sequence_id)
  
  # 4. Move barcode column to first
  cleaned_bcr_df <- cleaned_bcr_df[, c("barcode", setdiff(names(cleaned_bcr_df), "barcode"))]
  
  # 5. Merge using left_join
  merged_metadata <- left_join(seurat_metadata, cleaned_bcr_df, by = c("barcode_seurat" = "barcode"))
  
  if (sum(is.na(merged_metadata$barcode_seurat)) > 0) {
    warning(paste(sum(is.na(merged_metadata$barcode_seurat)), "cells in Seurat object did not match any BCR data."))
  }
  
  # 7. Ensure unique column names
  colnames(merged_metadata) <- make.unique(colnames(merged_metadata))
  
  # 8. Reset row names
  rownames(merged_metadata) <- merged_metadata$barcode_seurat
  
  # 9. Assign new metadata back to Seurat object
  seurat_obj@meta.data <- merged_metadata
  
  # 10. Output message
  message("BCR data successfully merged into Seurat metadata")
  
  # 11. Return updated Seurat object
  return(seurat_obj)
}

flu_H5_mouse.obj_bcr <- add_bcr_to_seurat(flu_H5_mouse.obj, flu_H5_mouse.bcr)

#Normalization and Data scaling
flu_H5_mouse.obj_bcr <- NormalizeData(flu_H5_mouse.obj_bcr)
flu_H5_mouse.obj_bcr <- FindVariableFeatures(flu_H5_mouse.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
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


# add bcr to rna
naive_mouse2.obj_bcr <- add_bcr_to_seurat(naive_mouse2.obj, naive_mouse2.bcr)

# Normalization and Data scaling
naive_mouse2.obj_bcr <- NormalizeData(naive_mouse2.obj_bcr)
naive_mouse2.obj_bcr <- FindVariableFeatures(naive_mouse2.obj_bcr, selection.method = "vst", nfeatures = 2000)

# save rds
saveRDS(naive_mouse2.obj_bcr, file = "./processed_seurat_objects/naive_mouse2.obj_bcr.processed.rds")


## 02 Results ------------------------------------------------
####Figure
features_B<- c("Cd19","Cd79a","Ebf1","Iglc3","Fcer2a")
features_T<- c("Cd3e", "Cd3g")
features_Macrophage<- c("C1qa", "Cd68","Trf")
features_Monocyte<- c("Fn1","Ace")
features_Neutrophil<- c("S100a9","S100a8","Retnlg")
features_NKC<- c("Ncr1","Klrk1")
feature_Erythroid <-c("Hba-a1", "Hba-a2","Hbq1b")
features_Granulocyte <- c("Ccl4","Csf1")
genes_alls = c(features_B, features_T, features_Macrophage, features_Neutrophil, features_NKC,
               feature_Erythroid,features_Granulocyte)
plot1 <- DotPlot(flu_mouse.combined, features = genes_alls, assay = "RNA")+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
print(plot1)
#clustercluster

#####umap -----
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

#scRNAi_QC@meta.data
flu_mouse.annotated_clusters <- flu_mouse.combined@meta.data
write.csv(flu_mouse.annotated_clusters,'./results/Table2.csv',row.names = T)

umap_data <- data.frame(
  UMAP_1 = flu_mouse.combined@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = flu_mouse.combined@reductions$umap@cell.embeddings[, 2],
  cell_type = flu_mouse.combined$spleen_cell_subpopulations)

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

plot2 <- ggplot() +
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
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)  )
print(plot2)

#####cell counts------
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
plot3 <- ggplot(plot_data, aes(x = orig.ident, y = n, fill = spleen_cell_subpopulations)) +
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
  
  # ---  ---
  # 0.15  0.2  0.25
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

print(plot3)

b_cell_seurat_object <- subset(
  flu_mouse.combined, 
  idents = c(1,2,3,4,5,6,10,11,12,13,14,16,17,18,21,26,28,29,30,32,33,34,36)
)
b_cell_counts_from_subset <- table(b_cell_seurat_object$orig.ident)
print(b_cell_counts_from_subset)

#####heatmap ------
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

# 1)  levels 
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
grp_cols <- grp_cols[levels(Idents(flu_mouse.combined))]
genes_use <- intersect(top8_c$gene, rownames(flu_mouse.combined))
# 2)  gradientn color_manual  drop=FALSE
plot4 <- DoHeatmap(
  object       = flu_mouse.combined,
  features     = genes_use,
  size         = 5,
  draw.lines   = FALSE,
  group.by     = NULL,          #  Idents
  group.colors = grp_cols
) +
  scale_fill_gradientn(colors = c("#421863", "#437F8C", "#F4E755")) +  # 
  scale_color_manual(                                             # ← 
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
print(plot4)

#####Dotplot------

plot5 <- DotPlot(flu_mouse.combined, features = genes_alls) +
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
print(plot5)

####umap for 5 samples------
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
      size = 0.5    ) +
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
  plot6 <- wrap_plot(plot_list, ncol = 3) +
    plot_layout(guides = 'collect') 
    plot_annotation( 
      title = 'UMAP Projections per Sample',
      theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
    )
  
} else {
  plot5 <- wrap_plot(plot_list, ncol = 3)
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = 'UMAP Projections per Sample',
      theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
    )
}
print(plot6)


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

plot7<- DimPlot(B_cell_subset_flu_mouse.integrated, group.by = "orig.ident", reduction = "umap") +
  ggtitle("After B cell subset integration")
plot7b<- DimPlot(B_cell_subset_flu_mouse.integrated, group.by = "seurat_clusters", reduction = "umap", label = TRUE) +
  ggtitle("After B cell subset integration")

####Before B cell interation

B_cell_subset_flu_mouse.combined <- FindNeighbors(B_cell_subset_flu_mouse.combined, dims = 1:30)
B_cell_subset_flu_mouse.combined <- FindClusters(
  B_cell_subset_flu_mouse.combined, 
  resolution = 0.6,
  graph.name = "integrated_snn"
)
B_cell_subset_flu_mouse.combined <- RunUMAP(B_cell_subset_flu_mouse.combined, dims = 1:30)
plot7c <- DimPlot(B_cell_subset_flu_mouse.combined, group.by = "orig.ident", reduction = "umap") +
  ggtitle("Before B cell subset integration")
plot7d <- DimPlot(B_cell_subset_flu_mouse.combined, group.by = "seurat_clusters", reduction = "umap", label = TRUE) +
  ggtitle("Before B cell subset integration")
plot7 <- (plot7+plot7b)/(plot7c+plot7d)
print(plot7)

####Bcell

DefaultAssay(B_cell_subset_flu_mouse.integrated) <- "RNA"
B_cell_subset_flu_mouse.integrated <- NormalizeData(B_cell_subset_flu_mouse.integrated)
B_cell_subset_flu_mouse.integrated <- FindVariableFeatures(B_cell_subset_flu_mouse.integrated, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(B_cell_subset_flu_mouse.integrated)
B_cell_subset_flu_mouse.integrated <- ScaleData(B_cell_subset_flu_mouse.integrated, features = scale.genes)
B_cell_subset_flu_mouse.integrated <- RunPCA(B_cell_subset_flu_mouse.integrated, features = VariableFeatures(B_cell_subset_flu_mouse.integrated))
plot8a <- ElbowPlot(B_cell_subset_flu_mouse.integrated, ndims=30, reduction="pca") 
plot8b <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "pca")
plot8 <- plot8a+plot8b
print(plot8)

B_cell_subset_flu_mouse.integrated <- FindNeighbors(B_cell_subset_flu_mouse.integrated, dims = 1:20)
B_cell_subset_flu_mouse.integrated <- FindClusters(B_cell_subset_flu_mouse.integrated, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
plot9 <- clustree(B_cell_subset_flu_mouse.integrated)
print(plot9)


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
    
    scale_color_manual(
      values = bcell_cluster_colors, 
      name = "B Cell\nSubcluster",      guide = guide_legend(override.aes = list(size = 4))    ) +
    
    labs(title = sample_id) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right", # <-- 
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  plot_list_bcell[[sample_id]] <- p_sample
}

plot10 <- wrap_plot(plot_list_bcell, ncol = 3)

plot10 <- plot10 +
  plot_annotation(
    title = 'B Cell Subset UMAP Projections per Sample',
    theme = theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5))
  )
print(plot10)

####FiguerS11------
B_cell_subset_flu_mouse.integrated <- FindNeighbors(B_cell_subset_flu_mouse.integrated, dims= 1:20)
B_cell_subset_flu_mouse.integrated <- FindClusters(B_cell_subset_flu_mouse.integrated, resolution= 0.6)
table(B_cell_subset_flu_mouse.integrated@meta.data$seurat_clusters)
metadata <- B_cell_subset_flu_mouse.integrated@meta.data
B_cell_subset_flu_mouse.cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(B_cell_subset_flu_mouse.cell_cluster,'./results/Table3.csv',row.names = F)

plot11a <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "umap", label = TRUE)+ 
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15))
plot11b <- DimPlot(B_cell_subset_flu_mouse.integrated, reduction = "umap",  group.by = "orig.ident")+ 
  theme(
    legend.title = element_text(family = "Arial", size = 18), 
    legend.text = element_text(family = "Arial", size = 18),              
    plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5), 
    axis.title = element_text(family = "Arial", size = 15), 
    axis.text = element_text(family = "Arial", size = 15))
plot11<- plot11a + plot11b
print(plot11)

B_cell_subset_flu_mouse.integrated <- JoinLayers(B_cell_subset_flu_mouse.integrated, assay = "RNA")
dim(B_cell_subset_flu_mouse.integrated[["RNA"]]$data) 
B_cell_subset_flu_mouse.integrated_markers.genes <- FindAllMarkers(B_cell_subset_flu_mouse.integrated, logfc.threshold = 0.25, test.use = "wilcox", assay = "RNA",
                                      min.pct = 0.25,  only.pos = FALSE)
saveRDS(B_cell_subset_flu_mouse.integrated, "./results/B_cell_subset_flu_mouse.integrated.rds")
saveRDS(B_cell_subset_flu_mouse.integrated_markers.genes, "./results/B_cell_subset_flu_mouse.integrated_markers.genes.rds")
B_cell_subset_flu_mouse.integrated <- readRDS(
  "./results/B_cell_subset_flu_mouse.integrated.rds"
)
#marker
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

plot12 <- DotPlot(B_cell_subset_flu_mouse.integrated, features = Bcell_genes_all1)+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
print(plot12)

##### B cell heatmap ------
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
    inherit.aes = FALSE, size = 0  # 0
  ) +
  scale_color_manual(
    values = id.cols,
    breaks = names(id.cols),
    drop   = FALSE,
    name   = "Identity"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),    fill  = guide_colorbar(barwidth = 1, barheight = 10, title.position = "top")
  )

print(plot1c)

ggsave(filename="./results/Figure1.pdf", plot = plot1c, width = 9, height = 16)
ggsave(filename="./results/Figure1.png", plot = plot1c, width = 9, height = 16)

#####B cell umap ------
B_cell_subset_flu_mouse.integrated@meta.data$B_cell_subpopulations <- plyr::mapvalues(
                                                                           from = c (0,7,1,4,11,13,14,19,20,21,23,24,2,22,3,5,6,8,9,10,12,15,16,17,18),
                                                                           to = c(rep("MZ",2), 
                                                                                  rep("GC",10),
                                                                                  rep("Bmem",2),
                                                                                  rep("PB",11)
                                                                                 ),
                                                                           x = B_cell_subset_flu_mouse.integrated@meta.data$seurat_clusters)

head(B_cell_subset_flu_mouse.integrated@meta.data)
table(B_cell_subset_flu_mouse.integrated@meta.data$B_cell_subpopulations)

umap_data <- data.frame(
  UMAP_1 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = B_cell_subset_flu_mouse.integrated@reductions$umap@cell.embeddings[, 2],
  cell_type = B_cell_subset_flu_mouse.integrated$B_cell_subpopulations)

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
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)  )
print(plot1a)
ggsave("./results/Figure2.pdf", plot = plot1a, width = 9, height = 8)
ggsave("./results/Figure2.png", plot = plot1a, width = 9, height = 8, dpi = 300)

##### Cell counts and percentage, Antigenic stimulation drives robust expansion of antigen-experienced B cell populations ------------------------------------------------

meta_df <- B_cell_subset_flu_mouse.integrated@meta.data
write.csv(meta_df,'./results/Table4.csv',row.names = F)

plot_data_grouped <- meta_df %>%
  # unlist
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
  #  dplyr::count
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
    size = 6,    fontface = "bold",
    inherit.aes = FALSE 
  ) +
  
  coord_flip() +
  scale_fill_manual(values = cell_type_colors) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, .25))) +
  labs(
    title = "B Cell Proportions by Group",
    x = "Experimental Group",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 18, face = "bold"), # Y
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )
print(plot1b)
ggsave("./results/Figure3.pdf", plot = plot1b, width = 9, height = 4.5)
ggsave("./results/Figure3.png", plot = plot1b, width = 9, height = 4.5, dpi = 300)

#####Dotplot, identification of distinct B cell functional subsets -----------------------------------------------
Bcell_genes_all1 = c( "Cr2","Cd1d1",
                      "Eif4a1","Mif","Ran",#PreGC
                      "Mki67","Hmgb2", "Tuba1b",#GC
                      "Bach2", "Spib","Cd19",#Bmem
                      "Sdc1","Slpi","Jchain","Prdm1")#PB

plot1d <- DotPlot(B_cell_subset_flu_mouse.integrated, features = Bcell_genes_all1,group.by = "B_cell_subpopulations") +
   scale_color_gradientn(
    colors = c("#313695", "#FFFFBF", "#A50026"), # () - () - ()
    name = "Average\nExpression"  ) +
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
ggsave("./results/Figure4.pdf", plot = plot1d, width = 18, height = 6)
ggsave("./results/Figure4.png", plot = plot1d, width = 18, height = 6)

# ####Split UMAPs,NaiveH1N1H5N1 -----------------------------------------------

obj <- B_cell_subset_flu_mouse.integrated

obj$group_name <- dplyr::case_when(
  grepl("^naive[_-]?mouse", obj$orig.ident, ignore.case = TRUE) ~ "naive",
  grepl("^flu[_-]?h1",      obj$orig.ident, ignore.case = TRUE) ~ "flu_H1",
  grepl("^flu[_-]?h5",      obj$orig.ident, ignore.case = TRUE) ~ "flu_H5",
  TRUE ~ "other"   # NA
)

# other
obj$group_name <- factor(obj$group_name, levels = c("naive","flu_H1","flu_H5","other"))

table(obj$group_name, useNA = "ifany")

umap_data <- data.frame(
  UMAP_1 = obj@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = obj@reductions$umap@cell.embeddings[, 2],
  cell_type = obj$B_cell_subpopulations,
  group = obj$group_name # <-- 
)

cell_type_medians <- umap_data %>%
  group_by(group, cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = 'drop'
  )

allcolour <- c("MZ" = "#3A8EBA", "GC" = "#E69F00", "Bmem" = "#009E73", "PB" = "#D55E00")

# ---  3:  UMAP  ---
plot1e <- ggplot() +
  geom_point(data = umap_data[, c("UMAP_1", "UMAP_2")], aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
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
  scale_color_manual(values = allcolour, name = "Cell Type") +  facet_wrap(~ group) + # "naive", "flu_H1", "flu_H5"
  theme_void() +  theme(
    legend.position = "right",    plot.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 20, face = "bold", family = "Arial"),    panel.border = element_rect(colour = "black", fill=NA, size=1),    legend.title = element_text(size = 18, face = "bold", family = "Arial"),    legend.text = element_text(size = 18, family = "Arial"),    legend.key.size = unit(1.5, "cm"),    legend.key.height = unit(1.2, "cm"),    legend.key.width = unit(1.2, "cm")     ) +
  guides(color = guide_legend(override.aes = list(size = 8)))

print(plot1e)
ggsave("./results/Figure5.pdf", plot = plot1e, width = 18, height = 6)
ggsave("./results/Figure5.png", plot = plot1e, width = 18, height = 6)

##### Pseudotime -----------------------------------------------
# 
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

# Monocle3 CellDataSet
mat <- GetAssayData(B_cell_subset_flu_mouse.integrated, slot = "data")  
cellInfo <- B_cell_subset_flu_mouse.integrated@meta.data  
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))  

cds <- new_cell_data_set(expression_data = mat,  
                         cell_metadata = cellInfo,  
                         gene_metadata = geneInfo)  

cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(B_cell_subset_flu_mouse.integrated, "umap")

cds = cluster_cells(cds, cluster_method = 'louvain')
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len=10
))

start = c("Cluster_A")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name

flag = closest_vertex[as.character(colData(cds)$T_cell_annotation) %in% start,]
flag = as.numeric(names(which.max(table( flag ))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds)

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
  #  - 
  theme_classic(base_family = "Arial", base_size = 20) +
  theme(
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    
    # ★★★  ★★★
    axis.title.x = element_text(      size = 20, 
      face = "bold", 
      color = "black",
      margin = margin(t = 15)    ),
    axis.title.y = element_text(      size = 20,
      face = "bold", 
      color = "black",
      margin = margin(r = 15)    ),
    
    # ★★★  ★★★
    axis.text.x = element_text(      size = 20, 
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(      size = 20, 
      color = "black",
      face = "bold"
    )
  )

time = pseudotime(cds)
time[time==Inf] = 0
time = data.frame(cell=names(time), pseudotime=time)
write.table(time, file="./results/pseudotime_data.txt", sep="\t", quote=F, col.names=T, row.names=F)
saveRDS(cds, file="./results/monocle3_result.rds")
ggsave(filename="./results/Figure6.pdf", plot = plot1f, width = 10, height = 8)
ggsave(filename="./results/Figure6.png", plot = plot1f, width = 10, height = 8)

#####  DEG for PB-----------------------------------------------

###1
PB_cells <- subset(B_cell_subset_flu_mouse.integrated, 
                   subset = B_cell_subpopulations == "PB")
PB_group_names <- case_when(
  PB_cells$orig.ident %in% c("naive_mouse1", "naive_mouse2") ~ "naive", # PB_cellsnaive
  PB_cells$orig.ident %in% c("flu_H1_mouse1", "flu_H1_mouse2") ~ "flu_H1",
  PB_cells$orig.ident == "flu_H5_mouse" ~ "flu_H5"
)
PB_cells$group_name <- PB_group_names
head(PB_cells@meta.data)
table(PB_cells$group_name)
#
#  'group_name' 
deg_PB_standard <- FindMarkers(PB_cells, 
                               ident.1 = "flu_H5", 
                               ident.2 = "flu_H1", 
                               group.by = "group_name",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)

head(deg_PB_standard)

###2
# 1. ID
h1_PB_cells <- WhichCells(PB_cells, expression = group_name == "flu_H1")
h5_PB_cells <- WhichCells(PB_cells, expression = group_name == "flu_H5")

n_h1 <- length(h1_PB_cells)
n_h5 <- length(h5_PB_cells)

target_n <- min(n_h1, n_h5) # 600

# 2. 
n_iterations <- 50 # 50100
PB_deg_results_list <- list()
for (i in 1:n_iterations) {
  set.seed(i)
  h5_sampled_cells <- sample(h5_PB_cells, size = target_n)
  combined_cells <- c(h1_PB_cells, h5_sampled_cells)
  temp_subset <- subset(PB_cells, cells = combined_cells)
  
  # FindMarkers
  # logfc.threshold  min.pct 0
  PB_deg_temp <- FindMarkers(temp_subset,
                          ident.1 = "flu_H5",
                          ident.2 = "flu_H1",
                          group.by = "group_name",
                          logfc.threshold = 0, # FC
                          min.pct = 0,                          verbose = FALSE)
  
  PB_deg_temp$gene <- rownames(PB_deg_temp)
  
  PB_deg_results_list[[i]] <- PB_deg_temp
  
  print(paste("Completed iteration:", i))
}

### ****

library(dplyr)
library(tibble)

# 1. 
PB_all_iterations_df <- bind_rows(PB_deg_results_list)

# 2.
# avg_log2FC  p_val_adj
PB_aggregated_results <- PB_all_iterations_df %>%
  group_by(gene) %>%
  summarise(
    mean_avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    mean_p_val_adj = mean(p_val_adj, na.rm = TRUE),
    n_significant = sum(p_val_adj < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# 3. 
#  -log10 of the mean adjusted p-value
PB_volcano_data <- PB_aggregated_results %>%
  mutate(
    neg_log10_padj = -log10(mean_p_val_adj)
  )

head(PB_volcano_data)

library(ggplot2)
library(ggrepel)

# 1. 
log2FC_threshold <- 1  # log2FC > 0.5  < -0.5
padj_threshold <- 0.05   #  -log10(0.05)  1.3

# 2.
PB_volcano_data <- PB_volcano_data %>%
  mutate(
    significance = case_when(
      mean_avg_log2FC > log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Upregulated in H5",
      mean_avg_log2FC < -log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Downregulated in H5",
      TRUE ~ "Not Significant"
    )
  )

table(PB_volcano_data$significance)

# 3. 
genes_to_label <- PB_volcano_data %>%
  filter(abs(mean_avg_log2FC) > log2FC_threshold & mean_p_val_adj < padj_threshold) %>%
  arrange(mean_p_val_adj) %>%
  head(30)

# 4. 
plot1g <- ggplot(PB_volcano_data, aes(x = mean_avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  
  scale_color_manual(values = c("Upregulated in H5" = "#D55E00", 
                                "Downregulated in H5" = "#3A8EBA", 
                                "Not Significant" = "grey80")) +
  
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey50") +
  
  # ggrepel
  geom_text_repel(data = genes_to_label, 
                  aes(label = gene),
                  size = 4,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  max.overlaps = Inf) +  labs(
    title = "Volcano Plot of PB cells: H5 vs H1",
    x = "Mean Log2 Fold Change",
    y = "-log10(Mean Adjusted P-value)",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, family = "Arial"),
    axis.title.x = element_text(size = 18, family = "Arial"),  # x 
    axis.title.y = element_text(size = 18, family = "Arial"),  # y 
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Arial
    legend.text = element_text(size = 16, family = "Arial"),  # Arial
    legend.margin = margin(t = 8),  
    text = element_text(size = 18, family = "Arial")  # Arial
  )

print(plot1g)
ggsave(filename="./results/Figure7.pdf", plot = plot1g, width = 8.5, height = 7)
ggsave(filename="./results/Figure7.png", plot = plot1g, width = 8.5, height = 7)

GC_cells <- subset(B_cell_subset_flu_mouse.integrated, 
                   subset = B_cell_subpopulations == "GC")
GC_group_names <- case_when(
  GC_cells$orig.ident %in% c("naive_mouse1", "naive_mouse2") ~ "naive", # PB_cellsnaive
  GC_cells$orig.ident %in% c("flu_H1_mouse1", "flu_H1_mouse2") ~ "flu_H1",
  GC_cells$orig.ident == "flu_H5_mouse" ~ "flu_H5"
)
GC_cells$group_name <- GC_group_names
head(GC_cells@meta.data)
table(GC_cells$group_name)
#
#  'group_name' 
deg_GC_standard <- FindMarkers(GC_cells, 
                               ident.1 = "flu_H5", 
                               ident.2 = "flu_H1", 
                               group.by = "group_name",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)

head(deg_GC_standard)

###2
# 1. ID
h1_GC_cells <- WhichCells(GC_cells, expression = group_name == "flu_H1")
h5_GC_cells <- WhichCells(GC_cells, expression = group_name == "flu_H5")

n_h1 <- length(h1_GC_cells)
n_h5 <- length(h5_GC_cells)

target_n <- min(n_h1, n_h5) # 600

# 2. 
n_iterations <- 50 # 50100
GC_deg_results_list <- list()
for (i in 1:n_iterations) {
  set.seed(i)
  h5_sampled_cells <- sample(h5_GC_cells, size = target_n)
  combined_cells <- c(h1_GC_cells, h5_sampled_cells)
  temp_subset <- subset(GC_cells, cells = combined_cells)
  
  # FindMarkers
  # logfc.threshold  min.pct 0
  GC_deg_temp <- FindMarkers(temp_subset,
                             ident.1 = "flu_H5",
                             ident.2 = "flu_H1",
                             group.by = "group_name",
                             logfc.threshold = 0, # FC
                             min.pct = 0,                             verbose = FALSE)
  
  GC_deg_temp$gene <- rownames(GC_deg_temp)
  
  GC_deg_results_list[[i]] <- GC_deg_temp
  
  print(paste("Completed iteration:", i))
}
### ****

library(dplyr)
library(tibble)

# 1. 
GC_all_iterations_df <- bind_rows(GC_deg_results_list)

# 2.
# avg_log2FC  p_val_adj
GC_aggregated_results <- GC_all_iterations_df %>%
  group_by(gene) %>%
  summarise(
    mean_avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
    mean_p_val_adj = mean(p_val_adj, na.rm = TRUE),
    n_significant = sum(p_val_adj < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# 3. 
#  -log10 of the mean adjusted p-value
GC_volcano_data <- GC_aggregated_results %>%
  mutate(
    neg_log10_padj = -log10(mean_p_val_adj)
  )

head(GC_volcano_data)

library(ggplot2)
library(ggrepel)

# 1. 
log2FC_threshold <- 1  # log2FC > 0.5  < -0.5
padj_threshold <- 0.05   #  -log10(0.05)  1.3

# 2.
GC_volcano_data <- GC_volcano_data %>%
  mutate(
    significance = case_when(
      mean_avg_log2FC > log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Upregulated in H5",
      mean_avg_log2FC < -log2FC_threshold & mean_p_val_adj < padj_threshold ~ "Downregulated in H5",
      TRUE ~ "Not Significant"
    )
  )

table(GC_volcano_data$significance)

# 3. 
GC_genes_to_label <- GC_volcano_data %>%
  filter(abs(mean_avg_log2FC) > log2FC_threshold & mean_p_val_adj < padj_threshold) %>%
  arrange(mean_p_val_adj) %>%
  head(30)

# 4. 
plot1h <- ggplot(GC_volcano_data, aes(x = mean_avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 1.5) +
  
  scale_color_manual(values = c(
    "Upregulated in H5" = "#D55E00", 
    "Downregulated in H5" = "#3A8EBA", 
    "Not Significant" = "grey80"
  )) +
  
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_threshold),
             linetype = "dashed", color = "grey50") +
  
  #  ggrepel ✅  x / y 
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
    axis.title.x = element_text(size = 18, family = "Arial"),  # x 
    axis.title.y = element_text(size = 18, family = "Arial"),  # y 
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold", family = "Arial"),  # Arial
    legend.text = element_text(size = 16, family = "Arial"),  # Arial
    legend.margin = margin(t = 8),  
    text = element_text(size = 18, family = "Arial")  # Arial
  )

print(plot1h)
ggsave(filename="./results/Figure8.pdf", plot = plot1h, width = 10, height = 8)
ggsave(filename="./results/Figure8.png", plot = plot1h, width = 10, height = 8.5,dpi = 300)


obj <- B_cell_subset_flu_mouse.integrated  

# 1)  condition
if (!"condition" %in% colnames(obj@meta.data)) {
  obj$condition <- NA_character_
  obj$condition[grepl("^naive_", obj$orig.ident)]   <- "naive"
  obj$condition[grepl("^flu_H1_", obj$orig.ident)]  <- "flu_H1"
  obj$condition[grepl("^flu_H5",  obj$orig.ident)]  <- "flu_H5"
  obj$condition <- factor(obj$condition, levels = c("naive","flu_H1","flu_H5"))
}

# 2) /
md <- obj@meta.data
getcol <- function(nm) if (nm %in% colnames(md)) md[[nm]] else rep(NA, nrow(md))

vh <- getcol("v_call_10x");      vl <- getcol("light_v_call_10x")
vh2<- getcol("junction_aa");  vl2<- getcol("light_junction_10x_aa")prod <- getcol("productive");     inframe <- getcol("vj_in_frame")
clon <- getcol("clone_id")

has_heavy <- (!is.na(vh)  & grepl("^IGH", vh))  | (!is.na(vh2)  & grepl("^IGH", vh2))
has_light <- (!is.na(vl)  & grepl("^IG[KL]", vl))| (!is.na(vl2) & grepl("^IG[KL]", vl2))

obj$BCR_paired <- has_heavy & has_light
obj$BCR_paired_strict <- obj$BCR_paired &
  (is.na(prod)   | prod) &
  (is.na(inframe)| inframe) &
  (is.na(clon)   | !is.na(clon))

table(obj$BCR_paired, useNA = "ifany")
table(obj$BCR_paired_strict, useNA = "ifany")

#  RNA–BCR 
obj_bcr <- subset(obj, subset = BCR_paired)# obj_bcr <- subset(obj, subset = BCR_paired_strict)  # 

# UMAP
if (!"umap" %in% SeuratObject::Reductions(obj_bcr)) {
  obj_bcr <- RunUMAP(obj_bcr, dims = 1:30)}

# ① 
DimPlot(obj_bcr, reduction = "umap", group.by = "condition", pt.size = 0.25)

# ②
DimPlot(obj_bcr, reduction = "umap", split.by = "condition", group.by = "B_cell_subpopulations", ncol = 3, pt.size = 0.25)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cell_type_colors <- c(
  "MZ" = "#3A8EBA", 
  "GC" = "#E69F00", 
  "PB" = "#D55E00", 
  "Bmem" = "#009E73"
)
#  Seurat  obj_bcr
conditions_ordered <- c("naive", "flu_H1", "flu_H5")
obj_bcr$condition <- factor(obj_bcr$condition, levels = conditions_ordered)
summary_df <- obj_bcr@meta.data %>%
  group_by(condition, B_cell_subpopulations) %>%
  summarise(count = n(), .groups = 'drop_last') %>%
  #  condition 
  mutate(proportion = count / sum(count)) %>%
  ungroup()

print(summary_df)

plot_list <- list()

for (cond in conditions_ordered) {
  
  # --- a.  condition  UMAP  ---
  obj_subset <- subset(obj_bcr, subset = condition == cond)
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
    theme(plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial")) +    NoLegend()
  
  # --- b.  condition  ---
  summary_subset <- summary_df %>% filter(condition == cond)
  
  barplot <- ggplot(summary_subset, aes(x = " ", y = proportion, fill = B_cell_subpopulations)) +
    geom_col(width = 0.3) +    
    geom_text(aes(label = scales::percent(proportion, accuracy = 1)),
              position = position_stack(vjust = 0.2), 
              color = "black", 
              fontface = "bold",
              size = 4,              family = "Arial"  #  Arial
    ) +
    
    scale_fill_manual(values = cell_type_colors, name = "B Cell") +    theme_void() +
    theme(legend.position = "right",          legend.text = element_text(family = "Arial", size = 16),
          legend.title = element_text(family = "Arial", size = 16, face = "bold"))
  
  # --- c.  UMAP  ---
  # UMAP  3  1
  combined_panel <- umap_plot + barplot + plot_layout(widths = c(3, 1), heights = c(1, 0.3))+
    theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
  
  plot_list[[cond]] <- combined_panel
}
# 3
plot3a <- wrap_plot(plot_list, ncol = 3)+plot_layout(guides = "collect")

print(plot3a)
ggsave("./results/Figure9.pdf", plot = plot3a, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure9.png", plot = plot3a, width = 15, height = 5, dpi = 300, device = "png")


obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(treatment_group = case_when(
    grepl("naive", orig.ident, ignore.case = TRUE) ~ "Naive",
    grepl("H1", orig.ident, ignore.case = TRUE)    ~ "Flu_H1",
    grepl("H5", orig.ident, ignore.case = TRUE)    ~ "Flu_H5",
    TRUE                                          ~ "Other"
  )) %>%
  mutate(treatment_group = factor(treatment_group, levels = c("Naive", "Flu_H1", "Flu_H5")))

cat(" ( Naive, Flu_H1, Flu_H5):\n")
print(levels(obj_bcr$treatment_group))

head(obj_bcr@meta.data$c_call)
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(isotype = case_when(
    c_call == "IGHM" ~ "IgM",
    c_call == "IGHD" ~ "IgD",
    c_call == "IGHG1" ~ "IgG1",
    c_call == "IGHG2B" ~ "IgG2B",    c_call == "IGHG2C" ~ "IgG2C",
    c_call == "IGHG3" ~ "IgG3",
    c_call == "IGHA" ~ "IgA",   #  IGHA
    c_call == "IGHE" ~ "IgE",
    TRUE ~ "Unknown" #  ( NA)  "Unknown"
  ))
print(table(obj_bcr$isotype))
isotype_colors <- c(
  "IgM" = "#5A5A5A",    # /
  "IgD" = "#0072B2",    # 
  "IgG1" = "#F0E442",   # 
  "IgG2B" = "#E69F00",   # 
  "IgG2C" = "#D55E00",   # 
  "IgG3" = "#CC79A7",   # 
  "IgA" = "#d62728",   # 
  "IgG3" = "#009E73",    # 
  "Unknown" = "#D3D3D3" # 
)
#  DimPlot 
plot3b <- DimPlot(
  obj_bcr,
  reduction = "umap",
  group.by = "isotype",           #  Seurat  Isotype 
  split.by = "treatment_group",   # Seurat
  cols = isotype_colors,  pt.size = 1,
  ncol = 3                        #  (3)
) +
  # ---  ggplot2  ---
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold") #  (e.g., "Naive")
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 5),
    ncol = 1,
    title = NULL
  ))+
  labs(title = "") 

print(plot3b)
ggsave("./results/Figure10.pdf", plot = plot3b, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure10.png", plot = plot3b, width = 15, height = 5, dpi = 300, device = "png")


naive_mouse.bcr$barcode <- sub("_contig_\\d+", "", naive_mouse.bcr$sequence_id)

library(dplyr)
library(stringr)

## 0)  treatment_group
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(treatment_group = case_when(
    str_detect(orig.ident, regex("^naive",  ignore_case = TRUE)) ~ "Naive",
    str_detect(orig.ident, regex("^flu_H1", ignore_case = TRUE)) ~ "Flu_H1",
    str_detect(orig.ident, regex("^flu_H5", ignore_case = TRUE)) ~ "Flu_H5",
    TRUE ~ "Other"
  ))

# # 1) meta  join key barcode_seurat -1
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  mutate(barcode_key = as.character(barcode_seurat))

# # 2) BCRnaive join key _contig_ -1
# #     barcode  contig max mean/median
bcr_naive_keyed <- naive_mouse.bcr %>%
  mutate(barcode_key = barcode) %>%
  group_by(barcode_key) %>%
  summarise(MU_FREQ_naive = max(MU_FREQ_HEAVY_TOTAL, na.rm = TRUE), .groups = "drop")

# # 3)  naive
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

## 4)  naive  NA
cat("Non-NA SHM counts by group:\n")
print(obj_bcr@meta.data %>%
        group_by(treatment_group) %>%
        summarise(non_NA = sum(!is.na(MU_FREQ_HEAVY_TOTAL)), .groups = "drop"))

## 5)  naive 
obj_bcr@meta.data %>%
  filter(treatment_group == "Naive") %>%
  select(orig.ident, barcode_seurat, MU_FREQ_HEAVY_TOTAL) %>%
  slice_head(n = 15)

## 6)  SHM NA
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
obj_bcr@meta.data <- obj_bcr@meta.data %>%
  #  1a: 
  mutate(treatment_group = case_when(
    grepl("^naive", orig.ident, ignore.case = TRUE) ~ "Naive",
    grepl("^flu_H1", orig.ident, ignore.case = TRUE) ~ "Flu_H1",
    grepl("^flu_H5", orig.ident, ignore.case = TRUE) ~ "Flu_H5",
    TRUE ~ "Other"
  )) %>%
  mutate(treatment_group = factor(treatment_group, levels = c("Naive", "Flu_H1", "Flu_H5"))) %>%
  
  #  1b:  SHM  NA  0
  mutate(MU_FREQ_HEAVY_TOTAL = replace_na(MU_FREQ_HEAVY_TOTAL, 0))

# "shm_freq_meta"
obj_bcr$shm_freq_meta <- obj_bcr$MU_FREQ_HEAVY_TOTAL

summary(obj_bcr$shm_freq_meta)
# 1.  UMAP 
umap_coords <- as.data.frame(Embeddings(obj_bcr, reduction = "umap"))

# 2. 
metadata_to_plot <- obj_bcr@meta.data %>%
  select(treatment_group, shm_freq_meta)
# 3.  data.frame
plotting_df <- cbind(umap_coords, metadata_to_plot)

head(plotting_df)
plotting_df <- plotting_df %>%
  arrange(shm_freq_meta)

plot3c <- ggplot(
  plotting_df, 
  aes(x = umap_1, y = umap_2, color = shm_freq_meta) #  x, y 
) +
  geom_point(size = 1.0) +
  
  # facet_wrap
  facet_wrap(~ treatment_group, ncol = 3) +
  
  scale_color_gradientn(
    colors = c("grey85", "#FEE0D2", "#DE2D26", "#A50F15"),
    name = "SHM Freq."
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),    strip.background = element_blank(),    legend.position = "right",    
    #  Arial
    legend.title = element_text(size = 16, family = "Arial", face = "bold"),  # Arial
    legend.text = element_text(size = 16, family = "Arial")  # Arial
  )

print(plot3c)
ggsave("./results/Figure11.pdf", plot = plot3c, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure11.png", plot = plot3c, width = 15, height = 5, dpi = 300, device = "png")


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

# meta.data .x/.y
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
    strip.background = element_blank(),   #  facet 
    legend.position = "right",
    axis.title  = element_text(size = 16, family = "Arial", face = "bold"),    axis.text   = element_text(size = 16, family = "Arial")  )
print(plot3d)

table_summary <- plot_df %>%
  filter(tag != "Other") %>%
  dplyr::count(group, tag, clone_id, subpop, name = "n_cells") %>%
  arrange(group, tag, desc(n_cells))

print(table_summary)
ggsave("./results/Figure12.pdf", plot = plot3d, width = 15, height = 5, device = "pdf")
ggsave("./results/Figure12.png", plot = plot3d, width = 15, height = 5, dpi = 300, device = "png")


flu_H1_mouse1.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH1_mouse1.merge.final.bcr.shm.tsv")
flu_H1_mouse2.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH1_mouse2.merge.final.bcr.shm.tsv")
flu_H5_mouse.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/fluH5_mouse.merge.final.bcr.shm.tsv")
naive_mouse1.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/naive_mouse1.merge.final.bcr.shm.tsv")
naive_mouse2.bcr <- read_tsv("./cell_reports/VDJ/annotation_results/naive_mouse2.merge.final.bcr.shm.tsv")

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

results_directory <- "./results/Chord_Diagrams_Final_Corrected"

# folder
if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE)
  message("folder: ", results_directory)
}

# ---  2:  ---
# R
list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

# ---  3:  ---
generate_chord_diagram_final_corrected <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- : ", sample_name, " ---")
  message("----------------------------------------------------")
  
  # 3.1 
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
    message("Warning: Sample '", sample_name, "' ")
    return(NULL)
  }
  
  # 3.2 
  data_matrix <- xtabs(value ~ from + to, data = pairing_freq)
  heavy_genes <- sort(unique(pairing_freq$from))
  light_genes <- sort(unique(pairing_freq$to))
  gene_order <- c(heavy_genes, light_genes)
  all_genes <- unique(gene_order)
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  grid_colors <- setNames(colors[1:length(all_genes)], all_genes)
  
  # 3.3 
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
    
    # --- !!!  !!! ---
    # ylim
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
  
  # 3.4 
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Chord_Diagram.pdf"))
  message("Generating PDF file: ", pdf_file_name)
  pdf(pdf_file_name, width = 12, height = 12)
  create_plot_commands()
  dev.off()
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Chord_Diagram.png"))
  message("Generating PNG file: ", png_file_name)
  png(png_file_name, width = 1800, height = 1800, res = 150)
  create_plot_commands()
  dev.off()
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

# ---  4:  ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_chord_diagram_final_corrected(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

message("\n\nAll samples processed! Please check '", results_directory, "' folder")

###top 30

# ---  0: R ---
library(readr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(circlize) # circlize

# ---  1:  ---
results_directory <- "./results/Chord_Diagrams_Top30_Heavy"

if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE)
  message("folder: ", results_directory)
}

# ---  2:  ---

list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

generate_chord_diagram_top30_heavy <- function(sample_data, sample_name, output_dir) {
  
  message("\n----------------------------------------------------")
  message("--- : ", sample_name, " ---")
  message("----------------------------------------------------")
  
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
    message("Warning: Sample '", sample_name, "' ")
    return(NULL)
  }
  
  # Top 30
  heavy_gene_usage <- pairing_freq %>%
    group_by(from) %>%
    summarise(total_value = sum(value), .groups = 'drop') %>%
    arrange(desc(total_value))
  
  top_heavy_genes <- head(heavy_gene_usage$from, 20)
  
  message(sprintf(" '%s'  %d  Top %d ",
                  sample_name, nrow(heavy_gene_usage), length(top_heavy_genes)))
  
  pairing_freq_filtered <- pairing_freq %>%
    filter(from %in% top_heavy_genes)
  
  if (nrow(pairing_freq_filtered) == 0) {
    message("Warning: Sample '", sample_name, "' Top 30")
    return(NULL)
  }
  
  data_matrix <- xtabs(value ~ from + to, data = pairing_freq_filtered)
  heavy_genes <- sort(unique(pairing_freq_filtered$from))
  light_genes <- sort(unique(pairing_freq_filtered$to))
  gene_order <- c(heavy_genes, light_genes)
  all_genes <- unique(gene_order)
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  grid_colors <- setNames(rep_len(colors, length.out = length(all_genes)), all_genes)
  
  create_plot_commands <- function() {
    circos.clear()
    circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2), start.degree = 90)
    
    chordDiagram(data_matrix, order = gene_order, grid.col = grid_colors,
                 annotationTrack = "grid", 
                 preAllocateTracks = list(
                   list(track.height = 0.05), # 1 ()
                   list(track.height = 0.15)  # 3 ()
                 ))
    
    circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {})
    
    circos.track(track.index = 2, bg.border = NA, panel.fun = function(x, y) {
      circos.text(get.cell.meta.data("xcenter"), 
                  get.cell.meta.data("ylim")[1] + mm_y(4),
                  get.cell.meta.data("sector.index"), 
                  facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), 
                  cex = 0.6)    })
  }
  
  save_plot_with_titles <- function(file_path, device, width, height, res = NA) {
    message(": ", file_path)
    
    if (device == "pdf") {
      pdf(file_path, width = width, height = height)
    } else {
      png(file_path, width = width, height = height, res = res)
    }
    
    par(mar = c(1, 1, 4, 1))    create_plot_commands()
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
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_chord_diagram_top30_heavy(
      sample_data = current_data, 
      sample_name = name, 
      output_dir = results_directory
    )
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}


# 2.1 5
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"))

# 2.2 Create 2 merged files
flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 2.3 7
list_of_samples <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

message("7")

# ---  3:  ---
results_directory <- "./results/Treemap"

tryCatch({
  if (!dir.exists(results_directory)) {
    dir.create(results_directory, recursive = TRUE)
    message("folder: ", results_directory)
  } else {
    message("folder: ", results_directory)
  }
}, error = function(e) {
  # recursive=TRUE
  message("folder")
  message("Original error message: ", e$message)
})

# ---  4:  ---
generate_publication_treemap <- function(sample_data, sample_name, output_dir) {
  message("\n--- Start processing sample: ", sample_name, " ---")
  
  processed_H_fre <- sample_data %>%
    tidyr::unnest(v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = v_call_10x) %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    count(heavy_v, name = "total_value") %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH")) %>%
    arrange(desc(total_value))
  
  if (nrow(processed_H_fre) == 0) {
    message("Warning: Sample '", sample_name, "' has no valid heavy chain data, skipping plot.")
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
  message("Generating PDF file: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 8, device = cairo_pdf) # cairo_pdf
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.png"))
  message("Generating PNG file: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 8, dpi = 300, device = "png")
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

# ---  5:  ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_publication_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

####VL density
# ---  3:  ---
results_directory <- "./results/Treemap_Light"

tryCatch({
  if (!dir.exists(results_directory)) {
    dir.create(results_directory, recursive = TRUE)
    message("folder: ", results_directory)
  } else {
    message("folder: ", results_directory)
  }
}, error = function(e) {
  # recursive=TRUE
  message("folder")
  message("Original error message: ", e$message)
})

# ---  4:  ---
generate_publication_treemap <- function(sample_data, sample_name, output_dir) {
  message("\n--- Start processing sample: ", sample_name, " ---")
  
  processed_H_fre <- sample_data %>%
    tidyr::unnest(light_v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = light_v_call_10x) %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    count(heavy_v, name = "total_value") %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH")) %>%
    arrange(desc(total_value))
  
  if (nrow(processed_H_fre) == 0) {
    message("Warning: Sample '", sample_name, "' has no valid heavy chain data, skipping plot.")
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
  message("Generating PDF file: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 8, device = cairo_pdf) # cairo_pdf
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Treemap.png"))
  message("Generating PNG file: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 8, dpi = 300, device = "png")
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

# ---  5:  ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_publication_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

####top30

# --- Step 1: Load, merge and prepare data ---
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
message("7")

# --- Step 2: Set output directory ---
# folder
results_directory <- "./results/Treemap_Heavy_Top30"
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)
message("Output directory ready: ", results_directory)

generate_heavy_chain_top30_treemap <- function(sample_data, sample_name, output_dir) {
  
  message("\n--- Start processing sample: ", sample_name, " ---")
  
  # 3.1 Top 30
  processed_H_fre <- sample_data %>%
    tidyr::unnest(v_call_10x, keep_empty = TRUE) %>%
    transmute(heavy_v = v_call_10x) %>%
    filter(!is.na(heavy_v) & heavy_v != "") %>%
    mutate(heavy_v = str_remove(heavy_v, "\\*.*")) %>%
    count(heavy_v, name = "total_value") %>%
    arrange(desc(total_value)) %>%
    # --- !!! : 30 !!! ---
    slice_head(n = 20) %>%
    mutate(label = str_replace(heavy_v, "IGHV", "VH"))
  
  if (nrow(processed_H_fre) == 0) {
    message("Warning: Sample '", sample_name, "' has no valid heavy chain data, skipping plot.")
    return(NULL)
  }
  
  # 3.2 
  num_colors <- nrow(processed_H_fre)
  color_vector <- RColorBrewer::brewer.pal.info %>%
    filter(category == "qual") %>%
    row.names() %>%
    lapply(function(pal) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal)) %>%
    unlist() %>%
    unique()
  color_palette <- setNames(rep_len(color_vector, length.out = num_colors), processed_H_fre$label)
  
  # 3.3 
  treemap_plot <- ggplot(processed_H_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6
    ) +
    scale_fill_manual(values = color_palette) +
    # Top 30
    labs(title = paste("Top 30 Heavy Chain V Gene Usage:", sample_name)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 3.4 
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Heavy_Top30_Treemap.pdf"))
  message("Generating PDF file: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 4, device = "pdf")
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Heavy_Top30_Treemap.png"))
  message("Generating PNG file: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 4, dpi = 300, device = "png")
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

# ---  4:  ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_heavy_chain_top30_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

# --- Light top 30 ---
# ---  2:  ---
results_directory <- "./results/Treemap_Light_Top30"
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)
message("Output directory ready: ", results_directory)

generate_light_chain_top30_treemap <- function(sample_data, sample_name, output_dir) {
  
  message("\n--- Start processing sample: ", sample_name, " (Light Chain) ---")
  
  # 3.1 Top 30
  processed_L_fre <- sample_data %>%
    # --- !!!  1:  light_v_call_10x  !!! ---
    tidyr::unnest(light_v_call_10x, keep_empty = TRUE) %>%
    transmute(light_v = light_v_call_10x) %>%
    filter(!is.na(light_v) & light_v != "") %>%
    mutate(light_v = str_remove(light_v, "\\*.*")) %>%
    count(light_v, name = "total_value") %>%
    arrange(desc(total_value)) %>%
    slice_head(n = 20) %>%
    mutate(label = case_when(
      grepl("IGKV", light_v) ~ str_replace(light_v, "IGKV", "VK"),
      grepl("IGLV", light_v) ~ str_replace(light_v, "IGLV", "VL"),
      TRUE ~ light_v    ))
  
  if (nrow(processed_L_fre) == 0) {
    message("Warning: Sample '", sample_name, "' ")
    return(NULL)
  }
  
  # 3.2 
  num_colors <- nrow(processed_L_fre)
  color_vector <- RColorBrewer::brewer.pal.info %>%
    filter(category == "qual") %>%
    row.names() %>%
    lapply(function(pal) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal)) %>%
    unlist() %>%
    unique()
  color_palette <- setNames(rep_len(color_vector, length.out = num_colors), processed_L_fre$label)
  
  # 3.3 
  treemap_plot <- ggplot(processed_L_fre, 
                         aes(area = total_value, fill = label, label = label)) +
    geom_treemap(color = "white", size = 1.5) +
    geom_treemap_text(
      fontface = "bold", color = "white", place = "centre",
      grow = TRUE, min.size = 6
    ) +
    scale_fill_manual(values = color_palette) +
    # --- !!!  3:  !!! ---
    labs(title = paste("Top 30 Light Chain V Gene Usage:", sample_name)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 3.4 
  pdf_file_name <- file.path(output_dir, paste0(sample_name, "_Light_Top30_Treemap.pdf"))
  message("Generating PDF file: ", pdf_file_name)
  ggsave(pdf_file_name, plot = treemap_plot, width = 12, height = 5, device = "pdf")
  
  png_file_name <- file.path(output_dir, paste0(sample_name, "_Light_Top30_Treemap.png"))
  message("Generating PNG file: ", png_file_name)
  ggsave(png_file_name, plot = treemap_plot, width = 12, height = 5, dpi = 300, device = "png")
  
  message("--- Sample '", sample_name, "' processing complete ---")
}

# ---  4:  ---
for (name in names(list_of_samples)) {
  current_data <- list_of_samples[[name]]
  tryCatch({
    generate_light_chain_top30_treemap(sample_data = current_data, sample_name = name, output_dir = results_directory)
  }, error = function(e) {
    message("Error processing sample '", name, "' error: ", e$message)
  })
}

base_path <- "./cell_reports/VDJ/annotation_results/"
setwd(base_path)

output_dir <- "./results/VDJ_Sankey_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

sample_list <- list(
  flu_H1_mouse1 = flu_H1_mouse1.bcr,
  flu_H1_mouse2 = flu_H1_mouse2.bcr,
  flu_H5_mouse = flu_H5_mouse.bcr,
  naive_mouse1 = naive_mouse1.bcr,
  naive_mouse2 = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined = naive_mouse.bcr
)

prepare_sankey_data <- function(data) {
  sankey_data <- data %>%
    select(v_call_10x, j_call_10x, light_v_call_10x, light_j_call_10x) %>%
    mutate(
      VH = gsub("\\*.*", "", v_call_10x),
      JH = gsub("\\*.*", "", j_call_10x),
      VL = gsub("\\*.*", "", light_v_call_10x),
      JL = gsub("\\*.*", "", light_j_call_10x)
    ) %>%
    filter(!is.na(VH) & !is.na(JH) & !is.na(VL) & !is.na(JL)) %>%
    # ID
    mutate(clone_id = paste0("clone_", row_number()))
  
  return(sankey_data)
}

generate_sankey_plot <- function(sample_data, sample_name, output_dir) {
  sankey_data <- prepare_sankey_data(sample_data)
  
  if (nrow(sankey_data) == 0) {
    message("Sample '", sample_name, "' has no valid V-J pairing data, skipping plot.")
    return(NULL)
  }
  
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
  
  # VJ
  v_genes <- unique(sankey_long$gene[sankey_long$gene_type == "V"])
  j_genes <- unique(sankey_long$gene[sankey_long$gene_type == "J"])
  
  v_palette <- colorRampPalette(brewer.pal(8, "Set1"))(length(v_genes))
  j_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(j_genes))
  
  gene_colors <- setNames(c(v_palette, j_palette), c(v_genes, j_genes))
  
  sankey_plot <- ggplot(sankey_long,
                        aes(x = axis, stratum = gene, alluvium = clone_id,
                            y = 1, label = gene, fill = gene)) +
    geom_flow(alpha = 0.7, color = "white", width = 1/3) +
    geom_stratum(width = 1/3, color = "black", linewidth = 0.2) +
    
    scale_fill_manual(values = gene_colors) +
    
    geom_text(stat = "stratum", size = 2.5) +
    
    scale_x_discrete(
      limits = c("VH", "JH", "VL", "JL"),
      labels = c("VH" = "Heavy V", "JH" = "Heavy J", 
                 "VL" = "Light V", "JL" = "Light J"),
      expand = c(0.15, 0.05)
    ) +
    
    labs(
      title = paste("BCR V-J Gene Pairing:", sample_name),
      x = "Gene Segment",
      y = "Frequency",
      fill = "Gene"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
      axis.text.x = element_text(hjust = 1, size = 15),
      axis.text.y = element_blank(),  # y
      axis.title.y = element_blank(),  # y
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)    )
  
  pdf_file <- file.path(output_dir, paste0(sample_name, "_sankey_plot.pdf"))
  png_file <- file.path(output_dir, paste0(sample_name, "_sankey_plot.png"))
  
  ggsave(pdf_file, plot = sankey_plot, width = 7, height = 14, device = "pdf")
  ggsave(png_file, plot = sankey_plot, width = 7, height = 14, dpi = 300, device = "png", bg = "white")
  
  message("Sankey plot generated for sample '", sample_name, "'")
  
  return(sankey_plot)
}

for (sample_name in names(sample_list)) {
  message("Processing sample: ", sample_name)
  tryCatch({
    generate_sankey_plot(sample_list[[sample_name]], sample_name, output_dir)
  }, error = function(e) {
    message("Error processing sample '", sample_name, "' : ", e$message)
  })
}


# --- Step 1: Load, merge and prepare data ---

# 1.1 
base_path <- "./cell_reports/VDJ/annotation_results/"
output_dir <- "./results/VDJ_Sankey_Top30_VH_Plots_Large"dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1.2 
message("Loading data...")
flu_H1_mouse1.bcr <- readr::read_tsv(file.path(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- readr::read_tsv(file.path(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr  <- readr::read_tsv(file.path(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr  <- readr::read_tsv(file.path(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr  <- readr::read_tsv(file.path(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# 1.3 
flu_H1_mouse.bcr <- dplyr::bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr  <- dplyr::bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

# 1.4 
sample_list <- list(
  flu_H1_mouse1         = flu_H1_mouse1.bcr,
  flu_H1_mouse2         = flu_H1_mouse2.bcr,
  flu_H5_mouse          = flu_H5_mouse.bcr,
  naive_mouse1          = naive_mouse1.bcr,
  naive_mouse2          = naive_mouse2.bcr,
  flu_H1_mouse_combined = flu_H1_mouse.bcr,
  naive_mouse_combined  = naive_mouse.bcr
)
message("")

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
  
  message(sprintf(" %d VH Top %d ", 
                  nrow(vh_usage), length(top_vh_genes)))
  
  sankey_data <- processed_data %>%
    filter(VH %in% top_vh_genes) %>%
    mutate(clone_id = paste0("clone_", row_number()))
  
  return(sankey_data)
}

generate_sankey_plot <- function(
    sample_data, sample_name, output_dir,
    # --- !!!  1:  !!! ---
    output_width = 10,        #  (: )
    base_height = 1.5,          #  (: )
    height_per_gene = 0.15,
    # --- !!!  2:  !!! ---
    title_size = 22,    axis_text_size = 18,      #  (Heavy V, etc.) 
    gene_label_size = 4) {
  
  if (nrow(sample_data) == 0) {
    message("Sample '", sample_name, "' V-J")
    return(NULL)
  }
  
  sankey_long <- sample_data %>%
    select(clone_id, VH, JH, VL, JL) %>%
    pivot_longer(cols = -clone_id, names_to = "axis", values_to = "gene") %>%
    mutate(axis = factor(axis, levels = c("VH", "JH", "VL", "JL")))
  
  all_genes <- unique(sankey_long$gene)
  num_colors <- length(all_genes)
  color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
  gene_colors <- setNames(color_palette, all_genes)
  
  sankey_plot <- ggplot(data = sankey_long,
                        aes(x = axis, stratum = gene, alluvium = clone_id,
                            y = 1, fill = gene, label = gene)) +
    geom_flow(alpha = 0.6, color = "white", width = 0.4) +
    geom_stratum(width = 0.4, color = "grey30") +
    geom_text(stat = "stratum", size = gene_label_size, color = "black") + # <--  gene_label_size
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
      plot.title = element_text(hjust = 0.2, size = title_size, face = "bold"), 
      axis.text.x = element_text(
        size = axis_text_size, 
        face = "bold",
        margin = margin(t = -30, unit = "pt") # -10
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  
  num_strata <- length(unique(sankey_long$gene))
  plot_height <- max(base_height, num_strata * height_per_gene) # <--  base_height  height_per_gene
  
  pdf_file <- file.path(output_dir, paste0(sample_name, "_Top30VH_Sankey.pdf"))
  png_file <- file.path(output_dir, paste0(sample_name, "_Top30VH_Sankey.png"))
  
  ggsave(pdf_file, plot = sankey_plot, width = output_width, height = plot_height, device = "pdf", limitsize = FALSE) # <--  output_width
  ggsave(png_file, plot = sankey_plot, width = output_width, height = plot_height, dpi = 300, device = "png", bg = "white", limitsize = FALSE) # <--  output_width
  
  message("Sankey plot generated for sample '", sample_name, "'")
}

for (name in names(sample_list)) {
  message("\n====================================================")
  message("--- : ", name, " ---")
  message("====================================================")
  
  tryCatch({
    data_for_plotting <- prepare_sankey_data(sample_list[[name]], top_n = 20)
    
    generate_sankey_plot(data_for_plotting, name, output_dir)
    
  }, error = function(e) {
    message("!!!!!! Error processing sample '", name, "': ", e$message, " !!!!!!")
  })
}

### Figure: ------------------------------------------------
#####Isotype ------------------------------------------------

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
  dplyr::count(group, isotype, name = "count") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(proportion))

all_isotypes_found <- unique(isotype_proportions$isotype)
isotype_levels <- sort(all_isotypes_found) 
isotype_proportions$isotype <- factor(isotype_proportions$isotype, levels = rev(isotype_levels))

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
ggsave("./results/Figure13.pdf", plot = plot4a, width = 6, height = 5, device = "pdf")
ggsave("./results/Figure13.png", plot = plot4a, width = 6, height = 5, dpi = 300, device = "png")

# 1.1  
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
  dplyr::count(group, isotype, name = "count") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(proportion))

all_isotypes_found <- unique(isotype_proportions$isotype)
isotype_levels <- sort(all_isotypes_found) 
isotype_proportions$isotype <- factor(isotype_proportions$isotype, levels = rev(isotype_levels))

group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")
isotype_proportions$group <- factor(isotype_proportions$group, levels = group_levels)
print(isotype_proportions) 

# ---  3: 100% ---
unique_isotypes <- levels(isotype_proportions$isotype)
color_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique_isotypes)),
  unique_isotypes
)

plot13 <- ggplot(isotype_proportions, aes(x = group, y = proportion, fill = isotype)) +
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

print(plot13)

#####Figure4b Shannon's Entropy -------------------

library(vegan)

# 5
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr  <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"),  col_types = cols(.default = "c"))
naive_mouse1.bcr  <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr  <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

# 2
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
    panel.background = element_rect(fill = "transparent", colour = NA),    plot.background = element_rect(fill = "transparent", colour = NA),    panel.grid = element_blank(),    axis.line = element_line(color = "black"),    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    axis.text = element_text(size = 16, face = "bold",family = "Arial"),
    axis.title = element_text(size = 18, face = "bold", family = "Arial"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Arial")
  )
print(plot4b)
ggsave("./results/Figure14.pdf", plot = plot4b, width = 6, height = 5)
ggsave("./results/Figure14.png", plot = plot4b, width = 6, height = 5)

#####CDR3 clones share in H1 and H5 -------------------
AA_COLUMN <- "junction_10x_aa"          # CDR3 AA 
SUBSET_COL     <- "B_cell_subpopulations"
COND_LEVELS    <- c("flu_H1","flu_H5")       #  H1  H5 naive
TOP_N          <- 10                          # CDR3 Other

meta_df1 <- obj_bcr@meta.data %>%
  tibble::rownames_to_column(var = "barcode") %>%
  tibble::as_tibble()

dat <- meta_df1 %>%
  dplyr::filter(!is.na(.data[[CDR3_AA_COLUMN]]), .data[[CDR3_AA_COLUMN]] != "") %>%
  dplyr::transmute(
    condition = factor(.data$condition),    subset    = .data[[SUBSET_COL]],    cdr3      = toupper(.data[[CDR3_AA_COLUMN]])  ) %>%
  dplyr::filter(condition %in% COND_LEVELS) %>%          #  H1  H5
  dplyr::mutate(condition = factor(condition, levels = COND_LEVELS))

shared_h1h5 <- dat %>%
  dplyr::distinct(condition, cdr3) %>%
  dplyr::add_count(cdr3, name = "n_groups") %>%
  dplyr::filter(n_groups == length(COND_LEVELS)) %>%  dplyr::pull(cdr3) %>%
  unique()

df_flow <- dat %>%
  dplyr::filter(cdr3 %in% shared_h1h5) %>%
  dplyr::count(condition, cdr3, name = "count") %>%      #  CDR3 
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
  # CDR3  0
  tidyr::complete(condition, cdr3, fill = list(count = 0, percentage = 0)) %>%
  dplyr::mutate(
    condition = factor(condition, levels = COND_LEVELS),
    cdr3 = forcats::fct_reorder(cdr3, percentage, .fun = max, .desc = TRUE)  )
pal_example <- c(
  "#9ACD66", "#8FCAD0", "#F39A89", "#A7B7D9", "#CFE6A0",
  "#F4D57E", "#D9B07E", "#AFBE77", "#D6C7E8", "#B390C6"
)

# CDR3 TOP_N≦10 n
n_lv <- nlevels(df_flow_top$cdr3)
pal_use <- if (n_lv <= length(pal_example)) {
  pal_example[seq_len(n_lv)]
} else {
  # >10
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
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(plot4c)

ggsave("./results/Figure15.pdf", plot = plot4c, width = 8, height = 5)
ggsave("./results/Figure15.png", plot = plot4c, width = 8, height = 5)

#####HCDR3 Length ------------------------------------------------
# --- Step 0: Install and load required R packages ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

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
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

# HCDR3
hcdr3_lengths <- all_data %>%
  filter(locus == "IGH") %>%
  select(group, cdr3_aa) %>%
  filter(!is.na(cdr3_aa) & cdr3_aa != "") %>%
  # HCDR3
  mutate(hcdr3_length = str_length(cdr3_aa)) %>%
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

sample_sizes <- hcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3 length calculation complete.")
print(head(hcdr3_lengths))
print(sample_sizes)

stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = hcdr3_lengths,
  method = "wilcox.test"
)
message("Statistical test results:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\nSignificant group differences (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("No significant group differences found at p.adj < 0.05.")
}

my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

color_palette <- c("Naive" = "# 0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") #

plot4d <- ggplot(hcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  scale_fill_manual(values = color_palette) +
  
  # P
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", #  *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 5, fontface = "bold", color = "black") +
  
  labs(
    title = "HCDR3 Length Distribution",
    y = "HCDR3 Length (aa)",
    x = "" # X
  ) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 30, 10)) +
  
  theme_classic() +
  theme(
    legend.position = "none",    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot4d)
ggsave("./results/Figure16.pdf", plot = plot4d, width = 6, height = 5)
ggsave("./results/Figure16.png", plot = plot4d, width = 6, height = 5,dpi = 300)

#####5 samples for HCDR3 Length-------
# ---  1: 5 ---

files_to_process <- list(
  "Flu_H1_Mouse1" = "fluH1_mouse1.merge.final.bcr.shm.tsv",
  "Flu_H1_Mouse2" = "fluH1_mouse2.merge.final.bcr.shm.tsv",
  "Flu_H5_Mouse" = "fluH5_mouse.merge.final.bcr.shm.tsv",
  "Naive_Mouse1" = "naive_mouse1.merge.final.bcr.shm.tsv",
  "Naive_Mouse2" = "naive_mouse2.merge.final.bcr.shm.tsv"
)

#  purrr::map_dfr 
all_data <- map_dfr(files_to_process, ~ read_tsv(.x, col_types = cols(.default = "c")), .id = "group")

# ---  2:  - HCDR3 ---

# X
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

hcdr3_lengths <- all_data %>%
  filter(locus == "IGH") %>%
  select(group, cdr3_aa) %>%
  filter(!is.na(cdr3_aa) & cdr3_aa != "") %>%
  mutate(hcdr3_length = str_length(cdr3_aa)) %>%
  mutate(group = factor(group, levels = group_levels))

sample_sizes <- hcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3 length calculation complete.")

# ---  2.5:  ---
message("\n----------------------------------------------------")
message("--- Performing HCDR3 length significance test across 5 samples ---")
message("----------------------------------------------------")
stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = hcdr3_lengths,
  method = "wilcox.test"
)
message("Statistical test results:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\nSignificant group differences (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("No significant group differences found at p.adj < 0.05.")
}
message("----------------------------------------------------\n")

# ---  3:  ggplot2  ggpubr  ---

message("Starting violin plot...")

my_comparisons <- list( 
  c("Naive_Mouse1", "Naive_Mouse2"),
  c("Naive_Mouse2", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Flu_H5_Mouse")
)

color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plot14 <- ggplot(hcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = color_palette) +
  
  # P
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, size = 5) +
  
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 6, fontface = "bold", color = "black") +
  
  labs(
    title = "HCDR3 Length Distribution (5 Samples)",
    y = "HCDR3 Length (aa)",
    x = ""
  ) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 30, 10)) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    # X
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot14)

#####LCDR3 Length ------------------------------------------------

# --- Step 0: Install and load required R packages ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(purrr)

# ---  1:  ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)
list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

# ---  2:  - HCDR3 ---

# LCDR3
lcdr3_lengths <- all_data %>%
  select(group, light_cdr3_aa) %>%
  filter(!is.na(light_cdr3_aa) & light_cdr3_aa != "") %>%
  # HCDR3
  mutate(hcdr3_length = str_length(light_cdr3_aa)) %>%
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

sample_sizes <- lcdr3_lengths %>%
  dplyr::count(group) %>%
  dplyr::rename(label = n)

message("HCDR3 length calculation complete.")
print(head(lcdr3_lengths))
print(sample_sizes)

message("\n----------------------------------------------------")
message("--- Performing group-wise HCDR3 length significance test ---")
message("----------------------------------------------------")
stat_test_results <- compare_means(
  formula = hcdr3_length ~ group,
  data = lcdr3_lengths,
  method = "wilcox.test"
)
message("Statistical test results:")
print(stat_test_results)
significant_comparisons <- stat_test_results %>% filter(p.adj < 0.05)
message("\nSignificant group differences (p.adj < 0.05):")
if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("No significant group differences found at p.adj < 0.05.")
}
message("----------------------------------------------------\n")

# ---  3:  ggplot2  ggpubr  ---

my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

color_palette <- c("Naive" = "# 0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") #

plot4e <- ggplot(lcdr3_lengths, aes(x = group, y = hcdr3_length, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  scale_fill_manual(values = color_palette) +
  
  # P
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", #  *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  geom_text(data = sample_sizes, aes(x = group, y = 2, label = label), 
            size = 5, fontface = "bold", color = "black") +
  
  labs(
    title = "LCDR3 Length Distribution",
    y = "LCDR3 Length (aa)",
    x = "" # X
  ) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 30, 10)) +
  
  theme_classic() +
  theme(
    legend.position = "none",    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(plot4e)
ggsave("./results/Figure17.pdf", plot = plot4e, width = 6, height = 5)
ggsave("./results/Figure17.png", plot = plot4e, width = 6, height = 5,dpi = 300)

#####H_SHM ----------------------
# --- Step 0: Install and load required R packages ---

packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr", "scales")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr) 
library(purrr)  
library(scales) 

# ---  1:  ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr <- bind_rows(flu_H1_mouse1.bcr, flu_H1_mouse2.bcr)
naive_mouse.bcr <- bind_rows(naive_mouse1.bcr, naive_mouse2.bcr)

list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr,
  "Naive" = naive_mouse.bcr,
  "Flu_H5" = flu_H5_mouse.bcr
)

all_data <- bind_rows(list_of_samples, .id = "group")

# ---  2:  - SHM ---
# Check if required columns exist
if (!"MU_FREQ_HEAVY_TOTAL" %in% names(all_data)) {
  stop("Error: 'MU_FREQ_HEAVY_TOTAL' column not found in merged data. Please ensure input files are correct and SHM calculated.")
}

shm_frequencies <- all_data %>%
  select(group, MU_FREQ_HEAVY_TOTAL) %>%
  filter(!is.na(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = c("Naive", "Flu_H1", "Flu_H5")))

print(head(shm_frequencies))

# ---  3:  ggplot2  ggpubr  ---

my_comparisons <- list( c("Naive", "Flu_H1"), c("Naive", "Flu_H5"), c("Flu_H1", "Flu_H5") )

color_palette <- c("Naive" = "# 0072B2", "Flu_H1" = "#E69F00", "Flu_H5" = "#D55E00") #

plot4f <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  
  scale_fill_manual(values = color_palette) +
  
  # P
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", #  *, **, ***
                     bracket.size = 0.6,
                     size = 6) +
  
  labs(
    title = "Heavy Chain Mutation Frequency",
    y = "Heavy chain mutation frequency (%)",
    x = "" # X
  ) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 30, 10)) +
  
  theme_classic() +
  theme(
    legend.position = "none",    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

# R
print(plot4f)
ggsave("./results/Figure18.pdf", plot = plot4f, width = 6, height = 5)
ggsave("./results/Figure18.png", plot = plot4f, width = 6, height = 5,dpi = 300)

# --- Step 0: Install and load required R packages ---
packages_to_install <- c("dplyr", "readr", "stringr", "ggplot2", "ggpubr", "purrr", "scales")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(dplyr); library(readr); library(stringr); library(ggplot2); library(ggpubr); library(purrr); library(scales)

# ---  1:  ---
base_path <- "./cell_reports/VDJ/annotation_results/"
flu_H1_mouse1.bcr <- read_tsv(paste0(base_path, "fluH1_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr <- read_tsv(paste0(base_path, "fluH1_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr <- read_tsv(paste0(base_path, "fluH5_mouse.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr <- read_tsv(paste0(base_path, "naive_mouse1.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr <- read_tsv(paste0(base_path, "naive_mouse2.merge.final.bcr.shm.tsv"), col_types = cols(.default = "c"))

list_of_samples <- list(
  "Flu_H1_Mouse1" = flu_H1_mouse1.bcr,
  "Flu_H1_Mouse2" = flu_H1_mouse2.bcr,
  "Flu_H5_Mouse" = flu_H5_mouse.bcr,
  "Naive_Mouse1" = naive_mouse1.bcr,
  "Naive_Mouse2" = naive_mouse2.bcr
)

# bind_rows  'group'
all_data <- bind_rows(list_of_samples, .id = "group")

message("5")

# ---  3:  - SHM ---
if (!"MU_FREQ_HEAVY_TOTAL" %in% names(all_data)) {
  stop("Error: 'MU_FREQ_HEAVY_TOTAL' column not found in data.")
}

# X
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

shm_frequencies <- all_data %>%
  # locus
  filter(locus == "IGH") %>%
  select(group, MU_FREQ_HEAVY_TOTAL) %>%
  filter(!is.na(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_HEAVY_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM frequency data processing complete.")

message("Starting violin plot...")

my_comparisons <- list( 
  c("Naive_Mouse1", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Naive_Mouse1", "Flu_H5_Mouse"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Naive_Mouse2")
)

color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plot15 <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
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
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 30, 10)) + # Y
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

print(plot15)

#####L_SHM -----------------------

base_path <- "./cell_reports/VDJ/annotation_results/light_shm_results/"
flu_H1_mouse1.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr.light <- read_tsv(paste0(base_path, "fluH5_mouse_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr.light <- read_tsv(paste0(base_path, "naive_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr.light <- read_tsv(paste0(base_path, "naive_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))

flu_H1_mouse.bcr.light <- bind_rows(flu_H1_mouse1.bcr.light, flu_H1_mouse2.bcr.light)
naive_mouse.bcr.light <- bind_rows(naive_mouse1.bcr.light, naive_mouse2.bcr.light)

list_of_samples <- list(
  "Flu_H1" = flu_H1_mouse.bcr.light,
  "Flu_H5" = flu_H5_mouse.bcr.light,
  "Naive" = naive_mouse.bcr.light
)

# 1.3  bind_rows  'group'
all_data <- bind_rows(list_of_samples, .id = "group")

# ---  2:  - SHM ---
if (!"MU_FREQ_LIGHT_TOTAL" %in% names(all_data)) {
  stop("Error: 'MU_FREQ_LIGHT_TOTAL' column not found in data.")
}

group_levels <- c("Naive", "Flu_H1", "Flu_H5")

shm_frequencies <- all_data %>%
  select(group, MU_FREQ_LIGHT_TOTAL) %>%
  filter(!is.na(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM frequency data processing complete.")

message("Starting violin plot...")

my_comparisons <- list( 
  c("Naive", "Flu_H1"), 
  c("Naive", "Flu_H5"),
  c("Flu_H1", "Flu_H5")
)

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
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) + # Y
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
ggsave("./results/Figure19.pdf", plot = plot4g, width = 6, height = 5)
ggsave("./results/Figure19.png", plot = plot4g, width = 6, height = 5,dpi = 300)

base_path <- "./cell_reports/VDJ/annotation_results/light_shm_results/"
flu_H1_mouse1.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H1_mouse2.bcr.light <- read_tsv(paste0(base_path, "fluH1_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))
flu_H5_mouse.bcr.light <- read_tsv(paste0(base_path, "fluH5_mouse_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse1.bcr.light <- read_tsv(paste0(base_path, "naive_mouse1_light_with_shm.tsv"), col_types = cols(.default = "c"))
naive_mouse2.bcr.light <- read_tsv(paste0(base_path, "naive_mouse2_light_with_shm.tsv"), col_types = cols(.default = "c"))

list_of_samples <- list(
  "Flu_H1_Mouse1" = flu_H1_mouse1.bcr.light,
  "Flu_H1_Mouse2" = flu_H1_mouse2.bcr.light,
  "Flu_H5_Mouse" = flu_H5_mouse.bcr.light,
  "Naive_Mouse1" = naive_mouse1.bcr.light,
  "Naive_Mouse2" = naive_mouse2.bcr.light
)

# bind_rows  'group'
all_data <- bind_rows(list_of_samples, .id = "group")

message("5")

# ---  3:  - SHM ---
if (!"MU_FREQ_LIGHT_TOTAL" %in% names(all_data)) {
  stop("Error: 'MU_FREQ_LIGHT_TOTAL' column not found in data.")
}

# X
group_levels <- c("Naive_Mouse1", "Naive_Mouse2", "Flu_H1_Mouse1", "Flu_H1_Mouse2", "Flu_H5_Mouse")

shm_frequencies <- all_data %>%
  # locus
  select(group, MU_FREQ_LIGHT_TOTAL) %>%
  filter(!is.na(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq = as.numeric(MU_FREQ_LIGHT_TOTAL)) %>%
  mutate(shm_freq_percent = shm_freq * 100) %>%
  mutate(group = factor(group, levels = group_levels))

message("SHM frequency data processing complete.")

message("Starting violin plot...")

my_comparisons <- list( 
  c("Naive_Mouse1", "Flu_H1_Mouse1"), 
  c("Flu_H1_Mouse1", "Flu_H1_Mouse2"),
  c("Naive_Mouse1", "Flu_H5_Mouse"),
  c("Flu_H1_Mouse1", "Flu_H5_Mouse"),
  c("Naive_Mouse1", "Naive_Mouse2")
)

color_palette <- c(
  "Naive_Mouse1" = "#86BBD8", "Naive_Mouse2" = "#33658A",
  "Flu_H1_Mouse1" = "#F6AE2D", "Flu_H1_Mouse2" = "#F26419",
  "Flu_H5_Mouse" = "#D55E00"
)

plot16 <- ggplot(shm_frequencies, aes(x = group, y = shm_freq_percent, fill = group)) +
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
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 30, 10)) + # Y
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

print(plot16)

## 0)  meta.data
obj <- B_cell_subset_flu_mouse.integrated

## 1)  meta.data  SHM 
# >0.02 = high(0,0.02] = low==0 = none
obj$MU_FREQ_HEAVY_TOTAL <- suppressWarnings(as.numeric(obj$MU_FREQ_HEAVY_TOTAL))

obj$SHM_group <- dplyr::case_when(
  obj$MU_FREQ_HEAVY_TOTAL > 0.02 ~ "high",
  obj$MU_FREQ_HEAVY_TOTAL > 0    & obj$MU_FREQ_HEAVY_TOTAL <= 0.02 ~ "low",
  obj$MU_FREQ_HEAVY_TOTAL == 0   ~ "none",
  TRUE ~ NA_character_
) %>% factor(levels = c("high","low","none"))

table(obj$SHM_group)
h1_GC <- subset(
  obj,
  subset = grepl("^flu_H1", orig.ident)&
    B_cell_subpopulations %in% c("GC")
)

## 3)  high/low  DEG
h1_GC <- subset(h1_GC, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h1_GC$SHM_group)

Idents(h1_GC) <- "SHM_group"

# RNA assay 

## 4) high vs low
deg_high_vs_low <- FindMarkers(
  h1_GC,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  #  flu_H1_* 
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 
LFC_cut  <- 0.25    # |log2FC| 
LABEL_N  <- 25       # Up/Down 

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), #  -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

#  LABEL_N 
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)#  FDR  |log2FC| 
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
    min.segment.length = 0, max.overlaps = Inf  ) +
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
ggsave("./results/Figure20.pdf", plot = plot4h, width = 8, height = 7)
ggsave("./results/Figure20.png", plot = plot4h, width = 8, height = 7,dpi = 300)


h1_PB <- subset(
  obj,
  subset = grepl("^flu_H1", orig.ident)&
    B_cell_subpopulations %in% c("PB")
)

## 3)  high/low  DEG
h1_PB <- subset(h1_PB, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h1_PB$SHM_group)

Idents(h1_PB) <- "SHM_group"

# RNA assay 

## 4) high vs low
deg_high_vs_low <- FindMarkers(
  h1_PB,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  #  flu_H1_* 
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 
LFC_cut  <- 0.25    # |log2FC| 
LABEL_N  <- 25       # Up/Down 

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), #  -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

#  LABEL_N 
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)#  FDR  |log2FC| 
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
ggsave("./results/Figure21.pdf", plot = plot4i, width = 8, height = 7)
ggsave("./results/Figure21.png", plot = plot4i, width = 8, height = 7,dpi = 300)


h5_GC <- subset(
  obj,
  subset = grepl("^flu_H5", orig.ident)&
    B_cell_subpopulations %in% c("GC")
)

## 3)  high/low  DEG
h5_GC <- subset(h5_GC, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h5_GC$SHM_group)

Idents(h5_GC) <- "SHM_group"

## 4) high vs low
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

FDR_cut  <- 0.05    # FDR 
LFC_cut  <- 0.25    # |log2FC| 
LABEL_N  <- 25       # Up/Down 

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), #  -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

#  LABEL_N 
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N)
)
set.seed(1)#  FDR  |log2FC| 
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
ggsave("./results/Figure22.pdf", plot = plot4j, width = 8, height = 7)
ggsave("./results/Figure22.png", plot = plot4j, width = 8, height = 7,dpi = 300)


h5_PB <- subset(
  obj,
  subset = grepl("^flu_H5", orig.ident)&
    B_cell_subpopulations %in% c("PB")
)

## 3)  high/low  DEG
h5_PB <- subset(h5_PB, subset = !is.na(SHM_group) & SHM_group %in% c("high","low"))
table(h5_PB$SHM_group)

Idents(h5_PB) <- "SHM_group"

## 4) high vs low
deg_high_vs_low <- FindMarkers(
  h5_PB,
  ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0.25, min.pct = 0.10,
  test.use = "wilcox",
  latent.vars = c("nCount_RNA", "percent.mt", "orig.ident")  #  flu_h5_* 
) %>%
  arrange(p_val_adj) %>%
  tibble::rownames_to_column("gene")

head(deg_high_vs_low, 20)
if (!"gene" %in% colnames(deg_high_vs_low)) {
  deg_high_vs_low <- tibble::rownames_to_column(deg_high_vs_low, "gene")
}

FDR_cut  <- 0.05    # FDR 
LFC_cut  <- 0.25    # |log2FC| 
LABEL_N  <- 25       # Up/Down 

df <- deg_high_vs_low %>%
  mutate(
    FDR      = p_val,                 # Seurat 
    log2FC   = avg_log2FC,
    negLog10 = -log10(pmax(FDR, 1e-300)), #  -log10(0)
    sig = case_when(
      FDR < FDR_cut & log2FC >=  LFC_cut ~ "Up",
      FDR < FDR_cut & log2FC <= - LFC_cut ~ "Down",
      TRUE                               ~ "NS"
    )
  )

n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
n_down <- sum(df$sig == "Down", na.rm = TRUE)

#  LABEL_N 
labels_df <- bind_rows(
  df %>% filter(sig == "Up")   %>% arrange(FDR) %>% head(LABEL_N),
  df %>% filter(sig == "Down") %>% arrange(FDR) %>% head(LABEL_N))

set.seed(1)
#  FDR  |log2FC| 
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
    max.overlaps = Inf                     # ← 
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
    plot.margin = margin(10, 24, 10, 10)  ) +
  coord_cartesian(clip = "off") +         # ← 
  annotate("text",
           x = min(df$log2FC, na.rm=TRUE), y = max(df$negLog10, na.rm=TRUE),
           hjust = 0, vjust = -0.5, label = paste0(n_down, " DOWN"), fontface = 2) +
  annotate("text",
           x = max(df$log2FC, na.rm=TRUE), y = max(df$negLog10, na.rm=TRUE),
           hjust = 1, vjust = -0.5, label = paste0(n_up, " UP"), fontface = 2)

print(plot4k)
ggsave("./results/Figure23.pdf", plot = plot4k, width = 8, height = 7)
ggsave("./results/Figure23.png", plot = plot4k, width = 8, height = 7,dpi = 300)