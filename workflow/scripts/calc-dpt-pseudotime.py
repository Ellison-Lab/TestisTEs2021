import scanpy as sc
import numpy as np
import random
import matplotlib.pyplot as plt

random.seed(snakemake.params['seed'])
np.random.seed(snakemake.params['seed'])

ad = sc.read(snakemake.input['ad'])

germ_key = snakemake.params['germ_key']

germline = ad.obs.clusters[ad.obs.clusters.str.contains(germ_key)].unique().astype(str).tolist()

root_key = snakemake.params['root_key']

ad = ad[ad.obs.clusters.str.contains(germ_key),]

sc.pp.neighbors(ad, n_pcs=15)

sc.tl.diffmap(ad, n_comps=5)

#sc.pp.neighbors(ad, use_rep='X_diffmap')

#sc.tl.paga(ad, groups='clusters')

#sc.pl.paga(ad, threshold=0.1, layout='eq_tree')

ad.uns['iroot'] = np.flatnonzero(ad.obs['clusters'].str.contains(root_key))[0]

sc.tl.dpt(ad, n_dcs=5)

ax = sc.pl.umap(ad, color=['clusters','dpt_pseudotime'], legend_loc='on data', gene_symbols='gene_symbol', use_raw=True, show=False)

# DEBUG = plt.savefig("test.png")
plt.savefig(snakemake.output["dps_umap"])

res = ad.obs.loc[:,'dpt_pseudotime']

res.index.rename("index",inplace=True)

res.to_csv(snakemake.output['csv'])
