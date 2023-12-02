import numpy as np

import scvelo as scv

adata = scv.datasets.pancreas(file_path="data/pancreas.h5ad")
obs_names = adata.obs_names.tolist()
obs_names[0] = f"blah-{obs_names[0]}"
obs_names[1] = f"blah-{obs_names[0]}"
adata.obs_names = obs_names

bdata = adata.copy()
bdata.obs_names = np.arange(0, bdata.n_obs).astype(str)

adata = scv.utils.merge(adata, bdata)
