import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams


def compute_dynamics(adata, basis, key='true'):
    from ..tools.dynamical_model_utils import unspliced, spliced, vectorize
    idx = np.where(adata.var_names == basis)[0][0] if isinstance(basis, str) else basis

    alpha = adata.var[key + '_alpha'][idx]
    beta = adata.var[key + '_beta'][idx]
    gamma = adata.var[key + '_gamma'][idx]
    t_ = adata.var[key + '_t_'][idx]

    t = adata.obs[key + '_t'].values if key is 'true' else adata.layers[key + '_t'][:, idx]
    tau, alpha, u0, s0 = vectorize(np.sort(t), t_, alpha, beta, gamma)

    ut = unspliced(tau, u0, alpha, beta)
    st = spliced(tau, s0, u0, alpha, beta, gamma)

    vt = spliced(tau + np.max(tau)/100, s0, u0, alpha, beta, gamma) - st
    return alpha, ut, st, vt


def show_full_dynamics(adata, basis, key='true'):
    color = 'grey' if key is 'true' else 'purple'
    linewidth = .5 if key is 'true' else 1

    _, ut, st, _ = compute_dynamics(adata, basis, key)
    pl.plot(st, ut, color=color, linewidth=linewidth)

    idx = np.where(adata.var_names == basis)[0][0]
    beta, gamma = adata.var[key + '_beta'][idx], adata.var[key + '_gamma'][idx]
    xnew = np.linspace(0, st.max())
    pl.plot(xnew, gamma / beta * xnew, color=color, linestyle='--', linewidth=linewidth)


def simulation(adata, var_names='all', ):
    from ..tools.utils import make_dense

    idx_sorted = np.argsort(adata.obs.true_t)
    data = adata[idx_sorted]

    t = data.obs.true_t
    var_names = adata.var_names if var_names is 'all' else [name for name in var_names if var_names in adata.var_names]

    figsize = rcParams['figure.figsize']
    ncols = len(var_names)
    for i, gs in enumerate(pl.GridSpec(1, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1]), dpi=100))):
        idx = np.where(data.var_names == var_names[i])[0][0]

        alpha, ut, st, _ = compute_dynamics(adata, idx)

        u = make_dense(data.layers['unspliced'][:, idx])
        s = make_dense(data.layers['spliced'][:, idx])

        pl.subplot(gs)

        pl.plot(t, alpha, label='alpha', linestyle='--', color='grey', linewidth=1)

        pl.plot(t, ut, label='unspliced', color='darkblue')
        pl.plot(t, st, label='spliced', color='red')

        pl.plot(t, u, color='darkblue', linewidth=.3)
        pl.plot(t, s, color='red', linewidth=.3)

        pl.xlim(0); pl.ylim(0); pl.legend()
        pl.xlabel('t'); pl.ylabel('counts')
