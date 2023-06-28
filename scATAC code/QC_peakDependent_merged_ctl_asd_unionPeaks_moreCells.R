# calculate peak-dependent QC metrics of the "merged_ctl_asd_unionPeaks_moreCells_withNS_TSS" object
# reference: script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", lines 242-268

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

# read in the merged object
combined <- readRDS("/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells_withNS_TSS.rds")
# add fraction of reads in unionPeaks
combined <- FRiP(object = combined, assay = 'unionPeaks', total.fragments = 'passed_filters')
combined$pct_reads_in_unionPeaks <- combined$FRiP * 100
# number of fragments in unionPeaks
combined$unionPeak_region_fragments <- combined$FRiP * combined$passed_filters
# blacklist ratio
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$unionPeak_region_fragments
# save the updated object
saveRDS(combined, file = "/u/home/s/shxiao/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/merged_ctl_asd_unionPeaks_moreCells_QC.rds")

# plot QC metrics by object (ctl and asd)
Idents(combined) <- combined$object
pdf("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/QC/merged_ctl_asd_unionPeaks_moreCells_QCplots.pdf", width = 15, height = 4)
VlnPlot(
  object = combined,
  features = c('pct_reads_in_unionPeaks', 'unionPeak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
dev.off()
# end of script