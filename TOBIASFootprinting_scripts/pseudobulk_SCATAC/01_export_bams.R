###############################################################################
# R script to merge bam reads for each cell type per disease condition
# Yuyan Cheng
# 05-10-2023
###############################################################################

###############################################################################
# Import packages and functions
###############################################################################

# R packages
liblist <- c("tidyverse", "readxl", "Seurat")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = FALSE)))

# to prepare barcode files for subsetting
prep_barcodes <- function(...) {
  cr <- list(...)
  dir.create(cr$subset_dir, showWarnings = FALSE)
  
  
  if(!file.exists(cr$barcode_path)){
    cell_ids_strip <- cr$data$barcode # this is barcode id
    bam_name <- cr$data$cluster # this is 
    
    barcode <- data.frame('cell_ids_strip' = cell_ids_strip, 'bam_name' = paste(cr$library_id, bam_name, sep = "_"))
    write_delim(barcode, cr$barcode_path, col_names=F, delim = '\t')}
  
}

# to obtain sorted bam files from cell ranger output directories
get_sample_tb <- function( cellranger_count_dirs) { 
  #sample_meta <- read_xlsx(meta_xlsx)
  sample_dirs <- list.dirs(cellranger_count_dirs, recursive = FALSE)
  sample_tb <- tibble(sample_dirs) %>%
    mutate(
      library_id = basename(sample_dirs),
      fragments = file.path(sample_dirs, "outs", "fragments.tsv.gz"),
      bam = file.path(sample_dirs, "outs", "possorted_bam.bam"),
      cells = file.path(sample_dirs, "outs", "singlecell.csv"),
    ) %>%
    filter(file.exists(fragments), file.exists(cells))
  return(sample_tb)
}

# to subset bam files from each sample by cell type
subset_bam_by_sample <- function(sample_meta_tb, sample_tb) { 
  args_tb <- left_join(sample_meta_tb, sample_tb, by = c("library_id")) %>%
    mutate(
      subset_cmd = str_glue("{SUBSET_BIN} filterbarcodes --bam {bam} --cells {barcode_path} --outdir {subset_dir} --nproc 4")
    )
  
  print(args_tb)
  
  walk(args_tb$subset_cmd, function(x) {
    print(x)
    if (system(x) != 0) { stop() }
  })
}

# to join the above subset bams by disease condition (contorl or asd) 
join_bams_by_cluster <- function(cluster_meta_tb) {
  
  args_tb <- cluster_meta_tb %>%
    mutate(cluster_aggr = pmap(., function(...) {
      cr <- list(...)
      this_meta_tb <- cr$data %>%
        group_by(library_id) %>%
        group_nest %>%
        mutate(subset_bam = file.path(out_dir, library_id, str_glue("{library_id}_{cr$cluster}.bam")))
      
      return(str_glue(
        "{SAMTOOLS_BIN}",
        "cat",
        "-o {cr$bam_path}",
        "--threads 4",
        "{paste(this_meta_tb$subset_bam, collapse = ' ')}",
        .sep = " "
      ))
    }))
  print(args_tb)
  walk(args_tb$cluster_aggr, function(x) {
    print(x)
    
    if (system(x) != 0) { stop() }
  })
}


cleanup_subset_dir <- function(cluster_meta_tb, sample_meta_tb) {
  cluster_meta_tb <- cluster_meta_tb %>%
    filter(file.exists(bam_path))
  unlink(sample_meta_tb$subset_dir, recursive = TRUE)
}

###############################################################################
#  Set up the environment for analysis
###############################################################################

base_dir <- normalizePath('/u/home/y/ycheng41/nobackup-dhg/Brie_ATAC')
out_dir <- file.path(base_dir, "1_bam_export", "merged_bam")
CELLRANGER_COUNT_DIRS <- "/u/project/geschwind/shxiao/ASD/cellRanger_atac/align_atac_cellR"

# import meta file derived from Signac object 
# contains each cell's barcode, disease condition, cell type, and library id. 
SAMPLE_META <- file.path(base_dir, "/0_resource/Brie_snATAC_meta.csv")

# use sinto and samtools for bam subsetting and merging
SUBSET_BIN <- system("which sinto", intern = TRUE)
SAMTOOLS_BIN <- system("which samtools", intern = TRUE)

###############################################################################
#  main function to merge reads from bam output from cellranger 
###############################################################################

main <- function() {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- read.csv(SAMPLE_META, row.names = 1)
  
  #create output directories for each sample (library id)
  sample_meta_tb <- meta %>%
    group_by(library_id) %>% 
    rename(cluster = obj_cellType) %>%
    group_nest %>% 
    mutate(
      subset_dir = file.path(out_dir, str_glue("{library_id}")),
      barcode_path = file.path(out_dir, str_glue("{library_id}/{library_id}_barcodes.csv"))
      ) 
  
  #prepare barcode files for each sample, with cell type each barcode belongs to
  barcode_tb <- pmap(sample_meta_tb, prep_barcodes)
  
  # obtain file paths of cellranger processed bams
  sample_tb <- get_sample_tb(CELLRANGER_COUNT_DIRS)
  
  # subset bam files given the barcode files
  subset_bam_by_sample(sample_meta_tb, sample_tb)
  
  # join bam files by cell type
  cluster_meta_tb <- meta %>%
    group_by(obj_cellType) %>% 
    rename(cluster = obj_cellType) %>%
    group_nest %>% 
    mutate(
      bam_path = file.path(out_dir, str_glue("{cluster}.bam"))
    ) 
  
  join_bams_by_cluster(cluster_meta_tb)
  
  #cleanup_subset_dir(cluster_tb)
}
if (!interactive()) main()
