from ..tools.dynamical_model_utils import unspliced, mRNA, vectorize, tau_inv, get_vars
from .utils import make_dense

import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams


def get_dynamics(adata, key='fit', extrapolate=False, sorted=False, t=None):
    alpha, beta, gamma, scaling, t_ = get_vars(adata, key=key)
    if extrapolate:
        u0_ = unspliced(t_, 0, alpha, beta)
        tmax = t_ + tau_inv(u0_ * 1e-4, u0=u0_, alpha=0, beta=beta)
        t = np.concatenate([np.linspace(0, t_, num=500), t_ + np.linspace(0, tmax, num=500)])
    elif t is None or t is True:
        t = adata.obs[key + '_t'].values if key is 'true' else adata.layers[key + '_t']

    tau, alpha, u0, s0 = vectorize(np.sort(t) if sorted else t, t_, alpha, beta, gamma)
    ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)

    return alpha, ut, st


def compute_dynamics(adata, basis, key='true', extrapolate=None, sort=True, t_=None, t=None):
    idx = np.where(adata.var_names == basis)[0][0] if isinstance(basis, str) else basis
    key = 'fit' if key + '_gamma' not in adata.var_keys() else key
    alpha, beta, gamma, scaling, t_ = get_vars(adata[:, basis], key=key)

    if 'fit_u0' in adata.var.keys():
        u0_offset, s0_offset = adata.var['fit_u0'][idx], adata.var['fit_s0'][idx]
    else:
        u0_offset, s0_offset = 0, 0

    if t is None or isinstance(t, bool) or len(t) < adata.n_obs:
        t = adata.obs[key + '_t'].values if key is 'true' else adata.layers[key + '_t'][:, idx]

    if extrapolate:
        u0_ = unspliced(t_, 0, alpha, beta)
        tmax = np.max(t) if True else tau_inv(u0_ * 1e-4, u0=u0_, alpha=0, beta=beta)
        t = np.concatenate([np.linspace(0, t_, num=500), np.linspace(t_, tmax, num=500)])

    tau, alpha, u0, s0 = vectorize(np.sort(t) if sort else t, t_, alpha, beta, gamma)

    ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)
    ut, st = ut * scaling + u0_offset, st + s0_offset
    return alpha, ut, st


def show_full_dynamics(adata, basis, key='true', use_raw=False, linewidth=1, show_assignments=None, ax=None):
    if ax is None: ax = pl.gca()
    color = 'grey' if key is 'true' else 'purple'
    linewidth = .5 * linewidth if key is 'true' else linewidth
    label = 'learned dynamics' if key is 'fit' else 'true dynamics'
    line = None

    if key is not 'true':
        _, ut, st = compute_dynamics(adata, basis, key, extrapolate=False, sort=False, t=show_assignments)
        if show_assignments is not 'only':
            ax.scatter(st, ut, color=color, s=1)
        if show_assignments is not None and show_assignments is not False:
            skey, ukey = ('spliced', 'unspliced') if use_raw or 'Ms' not in adata.layers.keys() else ('Ms', 'Mu')
            s, u = make_dense(adata[:, basis].layers[skey]).flatten(), make_dense(adata[:, basis].layers[ukey]).flatten()
            ax.plot(np.array([s, st]), np.array([u, ut]), color='grey', linewidth=.1 * linewidth)

    if show_assignments is not 'only':
        _, ut, st = compute_dynamics(adata, basis, key, extrapolate=True, t=show_assignments)
        line, = ax.plot(st, ut, color=color, linewidth=linewidth, label=label)

        idx = np.where(adata.var_names == basis)[0][0]
        beta, gamma = adata.var[key + '_beta'][idx], adata.var[key + '_gamma'][idx]
        xnew = np.linspace(np.min(st), np.max(st))
        ynew = gamma / beta * (xnew - np.min(xnew)) + np.min(ut)
        ax.plot(xnew, ynew, color=color, linestyle='--', linewidth=linewidth)
    return line, label


def simulation(adata, var_names='all', legend_loc='upper right', legend_fontsize=20, linewidth=None, dpi=None,
               xkey='true_t', ykey=['unspliced', 'spliced', 'alpha'], colors=['darkblue', 'darkgreen', 'grey'], **kwargs):
    from ..tools.utils import make_dense
    from .scatter import scatter
    var_names = adata.var_names if var_names is 'all' else [name for name in var_names if name in adata.var_names]

    figsize = rcParams['figure.figsize']
    ncols = len(var_names)
    for i, gs in enumerate(pl.GridSpec(1, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1]), dpi=dpi))):
        idx = np.where(adata.var_names == var_names[i])[0][0]

        alpha, ut, st = compute_dynamics(adata, idx)

        t = adata.obs[xkey] if xkey in adata.obs.keys() else make_dense(adata.layers['fit_t'][:, idx])
        idx_sorted = np.argsort(t)
        t = t[idx_sorted]

        ax = pl.subplot(gs)
        _kwargs = {'alpha': .3, 'title': '', 'xlabel': 'time', 'ylabel': 'counts'}
        _kwargs.update(kwargs)
        linewidth = 1 if linewidth is None else linewidth

        ykey = [ykey] if isinstance(ykey, str) else ykey
        for j, key in enumerate(ykey):
            if key in adata.layers:
                y = make_dense(adata.layers[key][:, idx])[idx_sorted]
                ax = scatter(x=t, y=y, color=colors[j], ax=ax, show=False, **_kwargs)

            if key is 'unspliced':
                ax.plot(t, ut, label='unspliced', color=colors[j], linewidth=linewidth)
            elif key is 'spliced':
                ax.plot(t, st, label='spliced', color=colors[j], linewidth=linewidth)
            elif key is 'alpha':
                ax.plot(t, alpha, label='alpha', linestyle='--', color=colors[j], linewidth=linewidth)

        pl.xlim(0)
        pl.ylim(0)
        if legend_loc is not 'none' and i == ncols-1:
            pl.legend(loc=legend_loc, fontsize=legend_fontsize)

