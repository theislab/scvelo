import numpy as np

import matplotlib.pyplot as pl
from matplotlib import rcParams

from scvelo.core import SplicingDynamics
from scvelo.tools._em_model_utils import get_vars, tau_inv, unspliced, vectorize
from .utils import make_dense


# TODO: Add docstrings
def get_dynamics(adata, key="fit", extrapolate=False, sorted=False, t=None):
    """TODO."""
    alpha, beta, gamma, scaling, t_ = get_vars(adata, key=key)
    if extrapolate:
        u0_ = unspliced(t_, 0, alpha, beta)
        tmax = t_ + tau_inv(u0_ * 1e-4, u0=u0_, alpha=0, beta=beta)
        t = np.concatenate(
            [np.linspace(0, t_, num=500), t_ + np.linspace(0, tmax, num=500)]
        )
    elif t is None or t is True:
        t = adata.obs[f"{key}_t"].values if key == "true" else adata.layers[f"{key}_t"]

    tau, alpha, u0, s0 = vectorize(np.sort(t) if sorted else t, t_, alpha, beta, gamma)
    ut, st = SplicingDynamics(
        alpha=alpha, beta=beta, gamma=gamma, initial_state=[u0, s0]
    ).get_solution(tau)
    return alpha, ut, st


# TODO: Add docstrings
def compute_dynamics(
    adata, basis, key="true", extrapolate=None, sort=True, t_=None, t=None
):
    """TODO."""
    idx = adata.var_names.get_loc(basis) if isinstance(basis, str) else basis
    key = "fit" if f"{key}_gamma" not in adata.var_keys() else key
    alpha, beta, gamma, scaling, t_ = get_vars(adata[:, basis], key=key)

    if "fit_u0" in adata.var.keys():
        u0_offset, s0_offset = adata.var["fit_u0"][idx], adata.var["fit_s0"][idx]
    else:
        u0_offset, s0_offset = 0, 0

    if t is None or isinstance(t, bool) or len(t) < adata.n_obs:
        t = (
            adata.obs[f"{key}_t"].values
            if key == "true"
            else adata.layers[f"{key}_t"][:, idx]
        )

    if extrapolate:
        u0_ = unspliced(t_, 0, alpha, beta)
        tmax = np.max(t) if True else tau_inv(u0_ * 1e-4, u0=u0_, alpha=0, beta=beta)
        t = np.concatenate(
            [np.linspace(0, t_, num=500), np.linspace(t_, tmax, num=500)]
        )

    tau, alpha, u0, s0 = vectorize(np.sort(t) if sort else t, t_, alpha, beta, gamma)

    ut, st = SplicingDynamics(
        alpha=alpha, beta=beta, gamma=gamma, initial_state=[u0, s0]
    ).get_solution(tau, stacked=False)
    ut, st = ut * scaling + u0_offset, st + s0_offset
    return alpha, ut, st


# TODO: Add docstrings
def show_full_dynamics(
    adata,
    basis,
    key="true",
    use_raw=False,
    linewidth=1,
    linecolor=None,
    show_assignments=None,
    ax=None,
):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    color = linecolor if linecolor else "grey" if key == "true" else "purple"
    linewidth = 0.5 * linewidth if key == "true" else linewidth
    label = "learned dynamics" if key == "fit" else "true dynamics"
    line = None

    if key != "true":
        _, ut, st = compute_dynamics(
            adata, basis, key, extrapolate=False, sort=False, t=show_assignments
        )
        if not isinstance(show_assignments, str) or show_assignments != "only":
            ax.scatter(st, ut, color=color, s=1)
        if show_assignments is not None and show_assignments is not False:
            skey, ukey = (
                ("spliced", "unspliced")
                if use_raw or "Ms" not in adata.layers.keys()
                else ("Ms", "Mu")
            )
            s, u = (
                make_dense(adata[:, basis].layers[skey]).flatten(),
                make_dense(adata[:, basis].layers[ukey]).flatten(),
            )
            ax.plot(
                np.array([s, st]),
                np.array([u, ut]),
                color="grey",
                linewidth=0.1 * linewidth,
            )

    if not isinstance(show_assignments, str) or show_assignments != "only":
        _, ut, st = compute_dynamics(
            adata, basis, key, extrapolate=True, t=show_assignments
        )
        (line,) = ax.plot(st, ut, color=color, linewidth=linewidth, label=label)

        idx = adata.var_names.get_loc(basis)
        beta, gamma = adata.var[f"{key}_beta"][idx], adata.var[f"{key}_gamma"][idx]
        xnew = np.linspace(np.min(st), np.max(st))
        ynew = gamma / beta * (xnew - np.min(xnew)) + np.min(ut)
        ax.plot(xnew, ynew, color=color, linestyle="--", linewidth=linewidth)
    return line, label


# TODO: Add docstrings
def simulation(
    adata,
    var_names="all",
    legend_loc="upper right",
    legend_fontsize=20,
    linewidth=None,
    dpi=None,
    xkey="true_t",
    ykey=None,
    colors=None,
    **kwargs,
):
    """TODO."""
    from scvelo.tools.utils import make_dense
    from .scatter import scatter

    if ykey is None:
        ykey = ["unspliced", "spliced", "alpha"]
    if colors is None:
        colors = ["darkblue", "darkgreen", "grey"]
    var_names = (
        adata.var_names
        if isinstance(var_names, str) and var_names == "all"
        else [name for name in var_names if name in adata.var_names]
    )

    figsize = rcParams["figure.figsize"]
    ncols = len(var_names)
    for i, gs in enumerate(
        pl.GridSpec(
            1, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1]), dpi=dpi)
        )
    ):
        idx = adata.var_names.get_loc(var_names[i])
        alpha, ut, st = compute_dynamics(adata, idx)
        t = (
            adata.obs[xkey]
            if xkey in adata.obs.keys()
            else make_dense(adata.layers["fit_t"][:, idx])
        )
        idx_sorted = np.argsort(t)
        t = t[idx_sorted]

        ax = pl.subplot(gs)
        _kwargs = {"alpha": 0.3, "title": "", "xlabel": "time", "ylabel": "counts"}
        _kwargs.update(kwargs)
        linewidth = 1 if linewidth is None else linewidth

        ykey = [ykey] if isinstance(ykey, str) else ykey
        for j, key in enumerate(ykey):
            if key in adata.layers:
                y = make_dense(adata.layers[key][:, idx])[idx_sorted]
                ax = scatter(x=t, y=y, color=colors[j], ax=ax, show=False, **_kwargs)

            if key == "unspliced":
                ax.plot(t, ut, label="unspliced", color=colors[j], linewidth=linewidth)
            elif key == "spliced":
                ax.plot(t, st, label="spliced", color=colors[j], linewidth=linewidth)
            elif key == "alpha":
                largs = {"linewidth": linewidth, "linestyle": "--"}
                ax.plot(t, alpha, label="alpha", color=colors[j], **largs)

        pl.xlim(0)
        pl.ylim(0)
        if legend_loc != "none" and i == ncols - 1:
            pl.legend(loc=legend_loc, fontsize=legend_fontsize)
