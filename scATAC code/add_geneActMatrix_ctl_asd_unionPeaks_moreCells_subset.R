# create a gene activity matrix in the harmony-integrated and subset ctl_asd_unionPeaks_morecells_scATAC object
# The merging and integration were on the union of ctl and asd peaks.
# Subset was based on pct_fragment_in_peaks, nucleosomal signal, and TSS score.
# reference: adapted from "add_geneActMatrix_ctl_asd_ctlPeaks40k_subset.R"
#            <https://satijalab.org/signac/articles/pbmc_vignette.html>
# anaconda3-r/4.0, h_data=250G, h_rt=24 hours

# load packages
library(Signac)
library(Seurat)

set.seed(1234)

# load in the integrated ctl_asd_unionPeaks scATAC object that has been subset and clustered with 49 lsi components
ctl_asd <- readRDS("~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/integrated_objects/seuObj_hmIntegrated_ctl_asd_unionPeaks_moreCells_subset_umap_clustered.rds")
    # the object was saved in script "integrate_ctl_asd_scATAC_moreCells_onUnionPeaks.R", line 383.
ctl_asd

# check default assay, which should be "unionPeaks"
DefaultAssay(ctl_asd)
# check the annotation in the "unionPeaks" assay
ctl_asd[["unionPeaks"]]

# check the path of the fragment files in the "unionPeaks" assay
# per Signac vignette "Data structures and object interaction"
GetFragmentData(object = Fragments(ctl_asd)[[1]], slot = "path")
GetFragmentData(object = Fragments(ctl_asd)[[2]], slot = "path")

# compute gene activity, i.e. counts per cell in gene body and promoter region.
gene_act <- GeneActivity(ctl_asd)    
    # Use the default setting. extract gene coordinates and extend them to include the 2 kb upstream region.
    # The gene activity matrix is calculated from the reads per cell in fragments files (see function description), 
    # thus I think it is independent to the peaks in chromatin assays, but changes with different cells.

# add the gene activity matrix to the Seurat object as a new assay
ctl_asd[['unionPeaks_RNA']] <- CreateAssayObject(counts = gene_act)
saveRDS(object = ctl_asd, file = "~/project-geschwind/ASD/Signac/merge_objects_byUnionPeaks/with_moreCells/integrated_objects/seuObj_hmintegrated_ctl_asd_unionPeaks_moreCells_subset_clustered_geneAct.rds")
# end of script
