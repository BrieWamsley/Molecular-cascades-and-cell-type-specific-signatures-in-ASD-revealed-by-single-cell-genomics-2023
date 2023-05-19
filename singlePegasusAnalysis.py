import pegasus as pg
import os
from matplotlib import pyplot as plt

sname=os.environ["i"]
print(sname)

data=pg.read_input("".join(['../rawData/filtered_feature_bc_matrix_',sname,'/matrix.mtx.gz']),genome="GRCh38")

# Get basic QC first
pg.qc_metrics(data,min_genes=250,percent_mito=20)
pg.qcviolin(data,plot_type='gene',dpi=60)
pg.qcviolin(data,plot_type='mito',dpi=60)

# get total cell number
cell_num=data.obs['Channel'].value_counts();print(cell_num)
print(type(cell_num))
# calculate double rate
drate=(0.0008*cell_num+0.0527)/100;print(drate)

# Filter low quality cells and use robustly expressed genes for dimensional reduction
pg.filter_data(data)
pg.identify_robust_genes(data)
pg.log_norm(data)
data.obs['Channel'].value_counts()
pg.highly_variable_features(data, consider_batch=False)
pg.pca(data, robust=True)
pc_regressed=pg.pc_regress_out(data, attrs=['n_genes', 'percent_mito'])
pg.neighbors(data)
pg.louvain(data)
pg.umap(data)
pg.scatter(data, attrs=['louvain_labels'],dpi=100)

# get marker genes for each cluster
pg.de_analysis(data, cluster='louvain_labels')
marker_dict = pg.markers(data)
celltype_dict = pg.infer_cell_types(data, markers = 'human_brain',output_file='annot_M9H3_72BW.txt')
celltype_dict
pg.run_scrublet(data,expected_doublet_rate=drate[0])
pg.scatter(data, attrs='louvain_labels', basis='umap',legend_loc='on data',dpi=100)
plt.savefig('Umap.pdf')
pg.scatter(data, attrs=['scrublet_score'], basis='umap',dpi=100)
plt.savefig('Umap_scrubletScore.pdf')
import pandas as pd
print(data.obs.groupby('louvain_labels')['scrublet_score'].mean())
print(data.uns['scrublet_stats'])
#data.obs['louvain_labels'].value_counts()
F=data.obs.groupby('louvain_labels')['scrublet_score'].mean()
F.to_csv("".join([sname,'_cluster_dScore.csv']))
import csv
with open("".join([sname,'_scrublet_stats.csv']), 'w') as f:
    for key in data.uns['scrublet_stats'].keys():
        f.write("%s,%s\n"%(key,data.uns['scrublet_stats'][key]))
# save scrublet data 
metaData=data.obs
A=metaData.index.to_list()
metaData.index=["".join([sname,"-"]) + sub for sub in A] 
metaData.to_csv("".join(['meta_',sname,'.csv']))
