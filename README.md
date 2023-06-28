# Molecular-cascades-and-cell-type-specific-signatures-in-ASD-revealed-by-single-cell-genomics-2023
Human postmortem single-cell integrative genomic analysis Wamsley et al Science 2023

# scRNA scripts: 

 Script to run single sample analysis to obtain doublets and basic QC:
  1)run_pegasas_singleSample.sh
  2)singlePegasusAnalysis.py
  3)combine_meta.R

 Integrated analysis using pegasus pipeline:
  1)Pegasus 1.0
  2)Pegasus_scASD_integration.py

 Differntial cell composition analysis using limma voom:
  1)Diff_Cell_Type_Compo_voom lm.R 

Pseudobulk differential gene expression analysis by cluster:
 1)Pseudo_bulk_analysis_DUP15q_balanced.R 
 2)Pseudo_bulk_analysis.R 
 3)Corr_table.R 
 
# scATAC scripts:

scATAC processing:
 1)featureMatrix_unionPeaks_ctl_and_asd_moreCells.R
 2)merge_ctl_and_asd_byUnionPeaks_moreCells.R
 3)QC_ns_tss_merged_ctl_asd_unionPeaks_moreCells.R
 4)QC_peakDependent_merged_ctl_asd_unionPeaks_moreCells.R
 5)subset_noTSS_inntegrate_ctl_asd_unionPeaks_moreCells.R
 6)add_geneActMatrix_ctl_asd_unionPeaks_moreCells_subeset.R

Footprinting:
 1)scATAcpseudo-bulk: export and filter bams files: 01_export_bams.R & 01_export_bams.R
 2)Footprinting: 03_TOBIAS_footprint.sh, 03_TOBIAS_BINDetect.sh, 03_TOBIAS_ATACCorrect.sh

 
