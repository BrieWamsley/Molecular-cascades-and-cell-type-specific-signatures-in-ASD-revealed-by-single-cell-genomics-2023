# run QC metrics (nucleosome signal and TSS enrichment) on the merged atac-seq object "ctl_asd_unionPeaks_moreCells"
# reference: script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", lines 156-231

# load libraries
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)    
library(ggplot2)
library(patchwork)

# read in the merged seurat object
combined <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells.rds")
combined
combined[["unionPeaks"]]
head(combined@meta.data)

# 1) add gene annotation to the 'unionPeaks' ChromatinAssay
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(combined[["unionPeaks"]]) <- annotations
# check the "unionPeaks" ChromatinAssay
combined[["unionPeaks"]]

# 2) compute QC metrics - nucleosome signal
# check metadata first
colnames(combined@meta.data)
# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)
# plot fragment length histogram (signal cutoff = 4)
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# save the nucleosome signal plot
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/QC/merged_ctl_asd_unionPeaks_moreCells_histoFragments.pdf", width = 8, height = 5)
FragmentHistogram(object = combined, group.by = 'nucleosome_group')
dev.off()
# count the number of cells with NS>4
table(combined$nucleosome_group)
# save the object with nucleosome signal calculated
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells_withNS.rds")

# 3) compute QC metrics - TSS enrichment score
# calculate TSS enrichment score
combined <- TSSEnrichment(object = combined, fast = FALSE)    
# plot the TSS enrichment scores (cutoff = 2)
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/QC/merged_ctl_asd_unionPeaks_moreCells_TSSplot.pdf", width = 8, height = 5)
TSSPlot(combined, group.by = 'high.tss') + NoLegend()
dev.off()
# save the object with TSS score calculated
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells_withNS_TSS.rds")
# count cells with high and low tss scores (cutoff = 2)
table(combined$high.tss)
table(combined$object, combined$high.tss)
# end of script
