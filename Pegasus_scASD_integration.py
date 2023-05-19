import pegasus as pg
data_subset1=pg.read_input("sample_rmDB.main.zarr.zip") # QC has been run for this data already
data_subset1.obs['Channel'].value_counts()
data_subset1.obs['doublet_call'].value_counts()

# Basic QC to look at nGene, mito, then filter cells 
pg.qc_metrics(data_subset1,min_genes=250,percent_mito=10)
pg.qcviolin(data_subset1,plot_type='gene',dpi=200,n_violin_per_panel=74,panel_size=(12,4))
pg.qcviolin(data_subset1,plot_type='mito',dpi=200,n_violin_per_panel=74,panel_size=(12,4))
pg.filter_data(data_subset1)
pg.identify_robust_genes(data_subset1)

# dimension reduction by pca using low number of components
pg.pca(data_subset1, robust=True,n_components=10)
pc_regressed=pg.pc_regress_out(data_subset1, attrs=['n_genes', 'percent_mito'])
harmony_key = pg.run_harmony(data_subset1,rep=pc_regressed)


# clustering and umap visualization
pg.neighbors(data_subset1,rep=harmony_key)
pg.louvain(data_subset1,rep=harmony_key,resolution=0.4)
pg.umap(data_subset1,rep=harmony_key)
pg.scatter(data_subset1, attrs=['louvain_labels', 'Diagnosis','BrainBrank','Clinical','COD1','COD2','sequencingBatch','Chemistry10X'], basis='umap')
pg.scatter(data_subset1, attrs=['louvain_labels'], basis='umap',legend_loc='on data',dpi=200)
pg.scatter(data_subset1, attrs=['percent_mito', 'n_genes'], basis='umap')


# Get marker genes for each cluster
pg.de_analysis(data_subset1, cluster='louvain_labels')
celltype_dict = pg.infer_cell_types(data_subset1, markers = 'human_brain',output_file='annot.txt')

# Revise cluster names
anno_dict=dict({'1':'EXN','2':'EXN','3':'ODC','4':'AST','5':'INN','6':'OPC','7':'MG','8':'END','9':'EXN','10':'NA'})
pg.annotate(data_subset1, name='anno', based_on='louvain_labels', anno_dict=anno_dict)
data_subset1.obs['anno'].value_counts()
pg.scatter(data_subset1, attrs='anno', basis='umap',legend_loc='on data',dpi=200)
pg.write_output(data_subset1,"sample.main.lowRes.rmDB.zarr.zip",file_type='zarr.zip')


# Second round to generate high resolution clusters
data_subset1=pg.read_input("sample.main.lowRes.rmDB.zarr.zip")

# remove low quanlity cells
from pegasusio import MultimodalData
data_clean=MultimodalData(data_subset1[~data_subset1.obs['anno'].isin(['NA']),].copy( ))
pg.scatter(data_clean, attrs='anno', basis='umap',legend_loc='on data',dpi=100)
pg.write_output(data_clean,"sample.main.clean.zarr.zip",file_type='zarr.zip')

# data_clean=pg.read_input("sample.main.clean.zarr.zip")
# repeat the same process with high n_components as first round without low quality cells
pg.qc_metrics(data_clean,min_genes=250,percent_mito=10)
pg.filter_data(data_clean)
pg.identify_robust_genes(data_clean)
pg.highly_variable_features(data_clean, consider_batch=True)
pg.pca(data_clean, robust=True,n_components=65)
pc_regressed=pg.pc_regress_out(data_clean, attrs=['n_genes', 'percent_mito'])
harmony_key = pg.run_harmony(data_clean,rep=pc_regressed)
pg.neighbors(data_clean,rep=harmony_key)
pg.louvain(data_clean,rep=harmony_key,resolution=2.0)
pg.umap(data_clean,rep=harmony_key)
pg.scatter(data_clean, attrs=['louvain_labels'], basis='umap',legend_loc='on data',dpi=200)
pg.scatter(data_clean, attrs=['percent_mito', 'n_genes'], basis='umap')
pg.de_analysis(data_clean, cluster='louvain_labels')
marker_dict = pg.markers(data_clean)
pg.write_results_to_excel(marker_dict, "result.de.clean.dim65.res2.0.xlsx")
celltype_dict = pg.infer_cell_types(data_clean, markers = 'human_brain',output_file='annot_clean.dim65.0.res2.0.txt')
pg.write_output(data_clean,"sample.main.clean.dim65.res2.0.zarr")
pg.write_output(data_clean,"sample.main.clean.dim65res2.0.h5ad",file_type='h5ad')
pg.scatter(data_clean, attrs=['anno'], basis='umap',legend_loc='on data',dpi=200)

