import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('data/my_data.h5ad')
# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('data/samples/mm-mg-blz-rep1.loom', cache=True)
ldata2 = scv.read('data/samples/mm-mg-blz-rep2.loom', cache=True)
ldata3 = scv.read('data/samples/mm-mg-blz-rep3.loom', cache=True)
ldata4 = scv.read('data/samples/mm-mg-blz-rep4.loom', cache=True)


# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_8' for bc in barcodes]
ldata4.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()

# concatenate four loom files
ldata = ldata1.concatenate([ldata2, ldata3,ldata4])

# make variable names unique in adata
adata.var_names_make_unique()

# check common var names and obs names in adata and ldata
common_obs = adata.obs_names.intersection(ldata.obs_names)
common_vars = adata.var_names.intersection(ldata.var_names)
print(len(common_obs), len(common_vars))

# clean all names in adata and ldata
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)


# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)


# plot umap to check
sc.pl.umap(adata, color='celltypes', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

# check proportions of spilec and unspliced counts in the dataset
scv.pl.proportions(adata, groupby='celltypes')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# visualize velocities
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltypes', save='embedding_grid.pdf', title='', scale=0.25)

# see how conditions change using velocity vectors
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltypes', 'condition'], save='embedding_stream.pdf', title='')

# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['Tmem119'], color='celltypes')

# downstream analysis
scv.tl.rank_velocity_genes(adata, groupby='celltypes', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

# check velocity confidence 
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

scv.pl.velocity_graph(adata, threshold=.1, color='celltypes')
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

# pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

# paga graph construction
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='celltypes')

df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)

 




