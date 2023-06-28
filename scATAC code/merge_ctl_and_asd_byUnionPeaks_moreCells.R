# merge the ctl and asd scATAC datasets by the union of peaks, using more cells from the cell ranger "singlecells.csv" files
# reference: script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", lines 29-44, 80-154

library(Signac)
library(Seurat)

set.seed(1234)

# 1) create Seurat objects
# extract metadata (per the Signac "pbmc" vignette)
# ctl metadata
md.ctl <- read.csv(
  file = "/u/project/geschwind/bwamsley/ATAC/CTL_AGR_10/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)    # copy the first column as "row.names"
dim(md.ctl)
# asd metadata
md.asd <- read.csv(
  file = "/u/project/geschwind/bwamsley/ATAC/ASD_AGR_12/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)
dim(md.asd)
# perform an initial filtering of low count cells (per the Signac vignette "merging datasets")
md.ctl <- md.ctl[md.ctl$passed_filters > 500, ]
dim(md.ctl)
md.asd <- md.asd[md.asd$passed_filters > 500, ]
dim(md.asd)

# for the ctl dataset - read in the feature matrix
ctl_unionP_counts <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/featureMatrix_unionPeaks_ctl_moreCells.rds")
dim(ctl_unionP_counts)    
# read in fragments object
frags.ctl <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/fragmentObject_unionPeaks_ctl_moreCells.rds")
frags.ctl
# create a chromatin assay
ctl_unionP_assay <- CreateChromatinAssay(
  counts = ctl_unionP_counts,
  fragments = frags.ctl,
  genome = 'hg38',
  min.cells = 1)    # follow the parameters in the "mouse brain" vignette, sep = c("-", "-").
ctl_unionP_assay
# create a seurat object
ctl_unionP_seu <- CreateSeuratObject(
  counts = ctl_unionP_assay,
  assay = "unionPeaks",
  project = "ATAC",
  meta.data = md.ctl)
ctl_unionP_seu
# save the ctl object
saveRDS(ctl_unionP_seu, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/seuObj_unionPeaks_ctl_moreCells.rds")

# for the asd dataset
asd_unionP_counts <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/featureMatrix_unionPeaks_asd_moreCells.rds")
dim(asd_unionP_counts)
frags.asd <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/fragmentObject_unionPeaks_asd_moreCells.rds")
frags.asd
asd_unionP_assay <- CreateChromatinAssay(
  counts = asd_unionP_counts,
  fragments = frags.asd,
  genome = 'hg38',
  min.cells = 1)
asd_unionP_assay
asd_unionP_seu <- CreateSeuratObject(
  counts = asd_unionP_assay,
  assay = "unionPeaks",
  project = "ATAC",
  meta.data = md.asd)
asd_unionP_seu
# save the asd object
saveRDS(asd_unionP_seu, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/prep_objects/seuObj_unionPeaks_asd_moreCells.rds")

# 2) merge objects
# check meta.data of the two object
colnames(ctl_unionP_seu@meta.data)
head(ctl_unionP_seu@meta.data, 3)
colnames(asd_unionP_seu@meta.data)
head(asd_unionP_seu@meta.data, 3)

# add new columns in metadata to identify object of origin
ctl_unionP_seu$object <- 'ctl'
asd_unionP_seu$object <- 'asd'

# merge datasets, adding a cell ID to make sure cell names are unique 
combined <- merge(
  x = ctl_unionP_seu,
  y = asd_unionP_seu,
  add.cell.ids = c("ctl", "asd")
)
combined[["unionPeaks"]]
# save the merged object
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells.rds")
# end of script