# create feature matrices of fragment counts in the union of ctl and asd scATAC peaks
# use more cells directly from the cell ranger "singlecell.csv" files
# reference: script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", lines 14-15, 17, 25-26, 29-44, 64-77
# anaconda3-r/4.0.3, h_data=64G,h_vmem=256G -pe shared 4 (4 cores).

library(Signac)
library(Seurat)
library(future)

# set "multicore" processing
plan("multicore", workers =4)
options(future.globals.maxSize = 50000 * 1024^2)    # for 50 Gb RAM

set.seed(1234)

## read in "combined.peaks" (same as that used in script "integrate_ctl_asd_scATAC_onUnionPeaks.R")
combined.peaks <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/preparative_objects/unionPeaks_ctl_asd_atac.rds")
combined.peaks

## load metadata 
md.ctl <- read.csv(
  file = "/u/project/geschwind/bwamsley/ATAC/CTL_AGR_10/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)
dim(md.ctl)

md.asd <- read.csv(
  file = "/u/project/geschwind/bwamsley/ATAC/ASD_AGR_12/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)
dim(md.asd)

# perform an initial filtering of low count cells
md.ctl <- md.ctl[md.ctl$passed_filters > 500, ]
dim(md.ctl)
md.asd <- md.asd[md.asd$passed_filters > 500, ]
dim(md.asd)

## read in fragment objects
frags.ctl <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/fragmentObject_unionPeaks_ctl_moreCells.rds")
frags.ctl
    # the object was saved in line 54 in the reference script)
frags.asd <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/fragmentObject_unionPeaks_asd_moreCells.rds")
frags.asd
    # the object was saved in line 62 in the reference script)

## quantify the combined.peaks in each dataset
ctl_unionP_counts <- FeatureMatrix(
  fragments = frags.ctl,
  features = combined.peaks,
  cells = rownames(md.ctl))
# save the feature matrix
saveRDS(ctl_unionP_counts, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/featureMatrix_unionPeaks_ctl_moreCells.rds")

asd_unionP_counts <- FeatureMatrix(
  fragments = frags.asd,
  features = combined.peaks,
  cells = rownames(md.asd))
saveRDS(asd_unionP_counts, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/featureMatrix_unionPeaks_asd_moreCells.rds")
# end of script
