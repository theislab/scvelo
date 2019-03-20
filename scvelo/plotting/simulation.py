from ..tools.dynamical_model_utils import unspliced, spliced, vectorize, tau_u
from .utils import make_dense

import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams


def compute_dynamics(adata, basis, key='true', extrapolate=None, sort=True):
    idx = np.where(adata.var_names == basis)[0][0] if isinstance(basis, str) else basis

    alpha = adata.var[key + '_alpha'][idx] if key + '_alpha' in adata.var.keys() else 1
    beta = adata.var[key + '_beta'][idx] if key + '_beta' in adata.var.keys() else 1
    gamma = adata.var[key + '_gamma'][idx]
    t_ = adata.var[key + '_t_'][idx]
    scaling = adata.var[key + '_scaling'][idx] if key + '_scaling' in adata.var.keys() else 1

    if extrapolate:
        u0_ = unspliced(t_, 0, alpha, beta)
        tmax = t_ + tau_u(u0_ * 1e-4, u0_, 0, beta)
        t = np.concatenate([np.linspace(0, t_, num=500), t_ + np.linspace(0, tmax, num=500)])
    else:
        t = adata.obs[key + '_t'].values if key is 'true' else adata.layers[key + '_t'][:, idx]
    tau, alpha, u0, s0 = vectorize(np.sort(t) if sort else t, t_, alpha, beta, gamma)

    ut = unspliced(tau, u0, alpha, beta) * scaling
    st = spliced(tau, s0, u0, alpha, beta, gamma)

    vt = spliced(tau + np.max(tau)/100, s0, u0, alpha, beta, gamma) - st
    return alpha, ut, st, vt


def show_full_dynamics(adata, basis, key='true', use_raw=False, linewidth=1):
    color = 'grey' if key is 'true' else 'purple'
    linewidth = .5 * linewidth if key is 'true' else linewidth

    _, ut, st, _ = compute_dynamics(adata, basis, key, extrapolate=True)
    pl.plot(st, ut, color=color, linewidth=linewidth)

    if key is not 'true':
        _, ut, st, _ = compute_dynamics(adata, basis, key, extrapolate=False, sort=False)
        pl.scatter(st, ut, color=color, s=1)
        if len(st) < 500:
            skey, ukey = ('spliced', 'unspliced') if use_raw or 'Ms' not in adata.layers.keys() else ('Ms', 'Mu')
            s, u = make_dense(adata[:, basis].layers[skey]), make_dense(adata[:, basis].layers[ukey])
            pl.plot(np.array([s, st]), np.array([u, ut]), color='grey', linewidth=.2 * linewidth)

    idx = np.where(adata.var_names == basis)[0][0]
    beta, gamma = adata.var[key + '_beta'][idx], adata.var[key + '_gamma'][idx]
    scaling = adata.var[key + '_scaling'][idx] if key + '_scaling' in adata.var.keys() else 1

    xnew = np.linspace(0, st.max())
    label = 'learned dynamics' if key is 'fit' else 'true dynamics'
    pl.plot(xnew, gamma / beta * scaling * xnew, color=color, linestyle='--', linewidth=linewidth, label=label)
    return label


def simulation(adata, var_names='all'):
    from ..tools.utils import make_dense
    var_names = adata.var_names if var_names is 'all' else [name for name in var_names if var_names in adata.var_names]

    figsize = rcParams['figure.figsize']
    ncols = len(var_names)
    for i, gs in enumerate(pl.GridSpec(1, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1]), dpi=100))):
        idx = np.where(adata.var_names == var_names[i])[0][0]

        alpha, ut, st, _ = compute_dynamics(adata, idx)

        t = adata.obs['true_t'] if 'true_t' in adata.obs.keys() else make_dense(adata.layers['fit_t'][:, idx])
        idx_sorted = np.argsort(t)
        t = t[idx_sorted]

        u = make_dense(adata.layers['unspliced'][:, idx])[idx_sorted]
        s = make_dense(adata.layers['spliced'][:, idx])[idx_sorted]

        pl.subplot(gs)

        pl.plot(t, alpha, label='alpha', linestyle='--', color='grey', linewidth=1)

        pl.plot(t, ut, label='unspliced', color='darkblue')
        pl.plot(t, st, label='spliced', color='red')

        pl.plot(t, u, color='darkblue', linewidth=.3)
        pl.plot(t, s, color='red', linewidth=.3)

        pl.xlim(0); pl.ylim(0); pl.legend()
        pl.xlabel('t'); pl.ylabel('counts')
