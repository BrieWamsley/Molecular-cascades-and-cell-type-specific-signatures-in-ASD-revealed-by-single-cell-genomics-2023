# subset the merged scATAC object of ctl_asd_unionPeaks started with cells in the "singlecells.csv" files
# without using the TSS enrichment score criterion.
# Then proceed with normalization, remension reduction, harmony integration, and clustering.
# reference: adapted from script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", lines 270-382
# run in a job with h_rt=48hours, h_data=250G

# load libraries
library(Signac)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

## 1) subset the merged_ctl_asd_object
# read in the merged object
combined <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells_QC.rds")
    # saved in script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", line 255.
combined

# subset without using TSS enrichment score as a criterion
combined <- subset(
  x = combined,
  subset = unionPeak_region_fragments > 3000 &
    unionPeak_region_fragments < 250000 &
    pct_reads_in_unionPeaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4
)
combined
# save the subset object
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/seuObj_merged_ctl_asd_unionPeaks_moreCells_subsetNoTSS.rds")

# QC plot after subset
Idents(combined) <- combined$object
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/QC/merged_ctl_asd_unionPeaks_moreCells_subsetNoTSS_QCplots.pdf", width = 15, height = 4)
VlnPlot(
  object = combined,
  features = c('pct_reads_in_unionPeaks', 'unionPeak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
dev.off()

## 2) normalization and dimension reduction of the merged object
set.seed(1234)    

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)    # per Signac vignette "merge datasets"
# linear dimension reduction
combined <- RunSVD(combined)
# plot the correlation between LSI components and sequencing depth
pdf("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/plots/DepthCorPlot_merged_ctl_asd_unionPeaks_moreCells_subsetNoTSS.pdf", width = 8, height = 5)
DepthCor(combined)
dev.off()
# non-linear dimension reduction. use 49 lsi.
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
# dimplot to show merging of datasets
pdf("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/plots/DimPlot_merged_ctl_asd_unionPeaks_subsetNoTSS_moreCells_byObject.pdf")
DimPlot(combined, group.by = 'object', pt.size = 0.1) + ggplot2::ggtitle("Merged")
dev.off()
# save the subset and umapped merged object
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/seuObj_merged_ctl_asd_unionPeaks_moreCells_subsetNoTSS_umap.rds")

## 3) harmony integration
set.seed(1234)

# check cell counts in each dataset
table(combined$object)

# integration
hm.integrated <- RunHarmony(
  object = combined,
  group.by.vars = 'object',
  reduction = 'lsi',
  assay.use = 'unionPeaks',
  project.dim = FALSE
)

# check for correlation between the harmony-corrected lsi reduced components and sequencing depth
pdf("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/plots/DepthCorPlot_hmIntegrated_ctl_asd_unionPeaks_moreCells_subsetNoTSS.pdf", width = 8, height = 5)
DepthCor(hm.integrated, reduction = 'harmony')
dev.off()

# re-compute UMAP using corrected LSI embeddings (reduction='harmony')
# remove the first LSI component, which normally reflects technical variation. use 49 lsi.
hm.integrated <- RunUMAP(hm.integrated, dims = 2:50, reduction = 'harmony')

# find neighbors. use 49 lsi.
set.seed(1234)
hm.integrated <- FindNeighbors(
  object = hm.integrated,
  reduction = 'harmony',
  dims = 2:50
)
# cluster with a resolution of 1.0 (- the command has an internal random seed.)
hm.integrated <- FindClusters(
  object = hm.integrated,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE
)

# dimplot after harmony_integration, color by 'Idents'/clusters
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/integrated_objects/plots/dimplot_hmIntegrated_ctl_asd_unionPeaks_moreCells_subsetNoTSS_res1.pdf")
DimPlot(object = hm.integrated, label = TRUE) + NoLegend()
dev.off()
# dimplot after harmony_integration, group by 'object'
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/integrated_objects/plots/dimplot_hmIntegrated_ctl_asd_unionPeaks_moreCells_subsetNoTSS_byObject.pdf")
DimPlot(hm.integrated, group.by = 'object', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
dev.off()

# save the harmony_integrated object after "RunUMAP" reduction and reclustering
saveRDS(hm.integrated, file = "~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/noTSS_subset/integrated_objects/seuObj_hmIntegrated_ctl_asd_unionPeaks_moreCells_subsetNoTSS_umap_clustered.rds")
# end of script.