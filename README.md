# Molecular-cascades-and-cell-type-specific-signatures-in-ASD-revealed-by-single-cell-genomics-2023
Human postmortem single-cell integrative genomic analysis Wamsley et al Science 2023

Files in this scRNA folder are: 

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
