# IN DEVELOPMENT

from .. import settings
from .. import logging as logg
from ..preprocessing.moments import get_connectivities
from .utils import make_unique_list, test_bimodality
from .dynamical_model_utils import BaseDynamics, mRNA, linreg, convolve, tau_inv, compute_dt, unspliced, spliced

import warnings
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams
from scipy.optimize import minimize, leastsq
from scipy.stats import t


class DynamicsRecovery(BaseDynamics):
    def __init__(self, adata=None, gene=None, load_pars=None, **kwargs):
        super(DynamicsRecovery, self).__init__(adata, gene, **kwargs)
        if load_pars and 'fit_alpha' in adata.var.keys():
            self.load_pars(adata, gene)
        else:
            self.initialize()

    def initialize(self):
        # set weights
        u, s, w, perc = self.u, self.s, self.weights, 98
        u_w, s_w,  = u[w], s[w]

        # initialize scaling
        self.std_u, self.std_s = np.std(u_w), np.std(s_w)
        scaling = self.std_u / self.std_s if isinstance(self.fit_scaling, bool) else self.fit_scaling
        u, u_w = u / scaling, u_w / scaling

        # initialize beta and gamma from extreme quantiles of s
        #print(self.gene)
        weights_s = s_w >= np.percentile(s_w, perc, axis=0)
        weights_u = u_w >= np.percentile(u_w, perc, axis=0)
        beta, gamma = 1, linreg(convolve(u_w, weights_s), convolve(s_w, weights_s))

        u_inf, s_inf = u_w[weights_u | weights_s].mean(), s_w[weights_s].mean()
        u0_, s0_ = u_inf, s_inf
        alpha = u_inf * beta  # np.mean([s_inf * gamma, u_inf * beta])  # np.mean([s0_ * gamma, u0_ * beta])

        # initialize switching from u quantiles and alpha from s quantiles
        tstat_u, pval_u, means_u = test_bimodality(u_w, kde=True)
        tstat_s, pval_s, means_s = test_bimodality(s_w, kde=True)
        self.pval_steady = max(pval_u, pval_s)
        self.u_steady = means_u[1]
        self.s_steady = means_s[1]

        if self.pval_steady < 1e-3:
            u_inf = np.mean([u_inf, self.u_steady])
            #s_inf = self.s_steady # np.mean([s_inf, self.s_steady])
            alpha = gamma * s_inf
            beta = alpha / u_inf

            #weights_u = u_w >= np.percentile(u_w, perc, axis=0)
            #u0_, s0_ = u_w[weights_u].mean(), s_w[weights_u].mean()
            u0_, s0_ = u_inf, s_inf

        # alpha, beta, gamma = np.array([alpha, beta, gamma]) * scaling
        t_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)

        # update object with initialized vars
        alpha_ = 0
        self.alpha, self.beta, self.gamma, self.alpha_, self.scaling = alpha, beta, gamma, alpha_, scaling
        self.u0_, self.s0_, self.t_ = u0_, s0_, t_
        self.pars = np.array([alpha, beta, gamma, self.t_, self.scaling])[:, None]

        # initialize time point assignment
        self.t, self.tau, self.o = self.get_time_assignment()
        self.loss = [self.get_loss()]

        self.initialize_scaling(sight=.5)
        self.initialize_scaling(sight=.1)

        self.steady_state_ratio = self.gamma / self.beta

        if self.high_pars_resolution:
            self.fit_t_and_alpha()
            self.fit_scaling_()
            self.fit_rates()
            self.fit_t_()
            self.fit_t_and_rates()

        # Overwrite callbacks
        if not self.high_pars_resolution:
            self.cb_fit_t_and_alpha = None
            self.cb_fit_scaling_ = None
            self.cb_fit_rates = None
            self.cb_fit_t_ = None
            self.cb_fit_t_and_rates = None
            self.cb_fit_rates_all = None

    def initialize_scaling(self, sight=.5):  # fit scaling and update if improved
        z_vals = self.scaling + np.linspace(-1, 1, num=5) * self.scaling * sight
        for z in z_vals:
            u0_ = self.u0_ * self.scaling / z
            beta = self.beta / self.scaling * z
            self.update(scaling=z, beta=beta, u0_=u0_)

    def fit(self, assignment_mode=None):
        if self.max_iter > 0:

            # pre-train with explicit time assignment
            self.fit_t_and_alpha()
            self.fit_scaling_()
            self.fit_rates()
            self.fit_t_()
            self.fit_t_and_rates()

            # train with optimal time assignment (oth. projection)
            self.assignment_mode = assignment_mode
            self.update(adjust_t_=False)
            self.fit_t_and_rates(refit_time=False)

        # self.update(adjust_t_=False)
        # self.t, self.tau, self.o = self.get_time_assignment()
        self.update()
        self.tau, self.tau_ = self.get_divergence(mode='tau')
        self.likelihood = self.get_likelihood(refit_time=False)

    def fit_t_and_alpha(self, **kwargs):
        alpha_vals = self.alpha + np.linspace(-1, 1, num=5) * self.alpha / 10
        for alpha in alpha_vals: self.update(alpha=alpha)

        def mse(x):
            return self.get_mse(t_=x[0], alpha=x[1], **kwargs)
        res = minimize(mse, np.array([self.t_, self.alpha]), callback=self.cb_fit_t_and_alpha, **self.simplex_kwargs)# method='Nelder-Mead')
        self.update(t_=res.x[0], alpha=res.x[1])

    def fit_rates(self, **kwargs):
        def mse(x):
            return self.get_mse(alpha=x[0], gamma=x[1], **kwargs)
        res = minimize(mse, np.array([self.alpha, self.gamma]), tol=1e-2, callback=self.cb_fit_rates, **self.simplex_kwargs)
        self.update(alpha=res.x[0], gamma=res.x[1])

    def fit_t_(self, **kwargs):
        def mse(x):
            return self.get_mse(t_=x[0], **kwargs)
        res = minimize(mse, self.t_, callback=self.cb_fit_t_, **self.simplex_kwargs)
        self.update(t_=res.x[0])

    def fit_rates_all(self, **kwargs):
        def mse(x):
            return self.get_mse(alpha=x[0], beta=x[1], gamma=x[2], **kwargs)
        res = minimize(mse, np.array([self.alpha, self.beta, self.gamma]), tol=1e-2, callback=self.cb_fit_rates_all, **self.simplex_kwargs)
        self.update(alpha=res.x[0], beta=res.x[1], gamma=res.x[2])

    def fit_t_and_rates(self, **kwargs):
        def mse(x):
            return self.get_mse(t_=x[0], alpha=x[1], beta=x[2], gamma=x[3], **kwargs)
        res = minimize(mse, np.array([self.t_, self.alpha, self.beta, self.gamma]), tol=1e-2,
                       callback=self.cb_fit_t_and_rates, **self.simplex_kwargs)
        self.update(t_=res.x[0], alpha=res.x[1], beta=res.x[2], gamma=res.x[3])

    def fit_scaling_(self, **kwargs):
        def mse(x):
            return self.get_mse(t_=x[0], beta=x[1], scaling=x[2], **kwargs)
        res = minimize(mse, np.array([self.t_, self.beta, self.scaling]), callback=self.cb_fit_scaling_, **self.simplex_kwargs)
        self.update(t_=res.x[0], beta=res.x[1], scaling=res.x[2])

    # Callback functions for the Optimizer
    def cb_fit_t_and_alpha(self, x):
        self.update(t_=x[0], alpha=x[1])

    def cb_fit_scaling_(self, x):
        self.update(t_=x[0], beta=x[1], scaling=x[2])

    def cb_fit_rates(self, x):
        self.update(alpha=x[0], gamma=x[1])

    def cb_fit_t_(self, x):
        self.update(t_=x[0])

    def cb_fit_t_and_rates(self, x):
        self.update(t_=x[0], alpha=x[1], beta=x[2], gamma=x[3])

    def cb_fit_rates_all(self, x):
        self.update(alpha=x[0], beta=x[1], gamma=x[2])

    def update(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None, adjust_t_=True):
        loss_prev = self.loss[-1] if len(self.loss) > 0 else 1e6

        alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
        t, tau, o = self.get_time_assignment(alpha, beta, gamma, scaling, t_, u0_, s0_, t)
        loss = self.get_loss(t, t_, alpha, beta, gamma, scaling)
        perform_update = loss < loss_prev

        on = self.o == 1
        if adjust_t_ and np.any(on):
            if not perform_update:
                alpha, beta, gamma, scaling, t_ = self.get_vars()
                t, tau, o = self.get_time_assignment()
                loss = self.get_loss()

            alt_t_ = t[on].max()
            if 0 < alt_t_ < t_:
                # alt_u0_, alt_s0_ = mRNA(alt_t_, 0, 0, alpha, beta, gamma)
                alt_t_ += np.max(t) / len(t) * np.sum(t == t_) # np.sum((self.u / self.scaling >= alt_u0_) | (self.s >= alt_s0_))
                alt_t, alt_tau, alt_o = self.get_time_assignment(alpha, beta, gamma, scaling, alt_t_)
                alt_loss = self.get_loss(alt_t, alt_t_, alpha, beta, gamma, scaling)
                ut_cur = unspliced(t_, 0, alpha, beta)
                ut_alt = unspliced(alt_t_, 0, alpha, beta)

                if alt_loss * .99 <= np.min([loss, loss_prev]) or ut_cur * .99 < ut_alt:
                    t, tau, o, t_, loss, perform_update = alt_t, alt_tau, alt_o, alt_t_, alt_loss, True

            if False:
                steady_states = t == t_
                if perform_update and np.any(steady_states):
                    t_ += t.max() / len(t) * np.sum(steady_states)
                    t, tau, o = self.get_time_assignment(alpha, beta, gamma, scaling, t_)
                    loss = self.get_loss(t, t_, alpha, beta, gamma, scaling)

        if perform_update:
            if u0_ is not None: self.u0_ = u0_
            if s0_ is not None: self.s0_ = s0_
            self.t, self.tau, self.o = t, tau, o
            self.alpha, self.beta, self.gamma, self.scaling, self.t_ = alpha, beta, gamma, scaling, t_
            self.pars = np.c_[self.pars, np.array([alpha, beta, gamma, t_, scaling])[:, None]]
            self.loss.append(loss)

        return perform_update


def read_pars(adata, pars_names=['alpha', 'beta', 'gamma', 't_', 'scaling', 'std_u', 'std_s', 'likelihood', 'u0', 's0'], key='fit'):
    pars = []
    for name in pars_names:
        pkey = key + '_' + name
        par = adata.var[pkey].values if pkey in adata.var.keys() else np.zeros(adata.n_vars) * np.nan
        pars.append(par)
    return pars


def write_pars(adata, pars, pars_names=['alpha', 'beta', 'gamma', 't_', 'scaling', 'std_u', 'std_s', 'likelihood', 'u0', 's0'], add_key='fit'):
    for i, name in enumerate(pars_names):
        adata.var[add_key + '_' + name] = pars[i]


def recover_dynamics(data, var_names='velocity_genes', max_iter=10, assignment_mode='projection', t_max=None,
                     fit_scaling=True, fit_time=True, fit_steady_states=True, fit_connected_states=True, use_raw=False,
                     load_pars=None, return_model=True, plot_results=False, add_key='fit', copy=False, **kwargs):
    """Estimates velocities in a gene-specific manner

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.

    Returns
    -------
    Returns or updates `adata`
    """
    adata = data.copy() if copy else data
    logg.info('recovering dynamics', r=True)

    if isinstance(var_names, str) and var_names not in adata.var_names:
        var_names = adata.var_names[adata.var[var_names] == True] if 'genes' in var_names and var_names in adata.var.keys() \
            else adata.var_names if 'all' in var_names or 'genes' in var_names else var_names
    var_names = make_unique_list(var_names, allow_array=True)
    var_names = [name for name in var_names if name in adata.var_names]
    if len(var_names) == 0:
        raise ValueError('Variable name not found in var keys.')

    if fit_connected_states is True:
        fit_connected_states = get_connectivities(adata)

    alpha, beta, gamma, t_, scaling, std_u, std_s, likelihood, u0, s0 = read_pars(adata)
    idx, L, P = [], [], []
    T = adata.layers['fit_t'] if 'fit_t' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    Tau = adata.layers['fit_tau'] if 'fit_tau' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    Tau_ = adata.layers['fit_tau_'] if 'fit_tau_' in adata.layers.keys() else np.zeros(adata.shape) * np.nan

    progress = logg.ProgressReporter(len(var_names))
    for i, gene in enumerate(var_names):
        dm = DynamicsRecovery(adata, gene, use_raw=use_raw, load_pars=load_pars, max_iter=max_iter, fit_time=fit_time,
                              fit_scaling=fit_scaling, fit_steady_states=fit_steady_states,
                              fit_connected_states=fit_connected_states, **kwargs)
        dm.fit(assignment_mode=assignment_mode)

        ix = np.where(adata.var_names == gene)[0][0]
        idx.append(ix)

        T[:, ix], Tau[:, ix], Tau_[:, ix] = dm.t, dm.tau, dm.tau_
        alpha[ix], beta[ix], gamma[ix], t_[ix], scaling[ix] = dm.pars[:, -1]
        u0[ix], s0[ix] = dm.u0, dm.s0
        std_u[ix], std_s[ix], likelihood[ix] = dm.std_u, dm.std_s, dm.likelihood
        L.append(dm.loss)
        if plot_results and i < 4:
            P.append(np.array(dm.pars))

        progress.update()
    progress.finish()

    write_pars(adata, [alpha, beta, gamma, t_, scaling, std_u, std_s, likelihood, u0, s0])
    adata.layers['fit_t'] = T
    adata.layers['fit_tau'] = Tau
    adata.layers['fit_tau_'] = Tau_

    cur_len = adata.varm['loss'].shape[1] if 'loss' in adata.varm.keys() else 2
    max_len = max(np.max([len(l) for l in L]), cur_len)
    loss = np.ones((adata.n_vars, max_len)) * np.nan

    if 'loss' in adata.varm.keys():
        loss[:, :cur_len] = adata.varm['loss']

    loss[idx] = np.vstack([np.concatenate([l, np.ones(max_len-len(l)) * np.nan]) for l in L])
    adata.varm['loss'] = loss

    if t_max is not False:
        dm, m = align_dynamics(adata, t_max=t_max, dm=dm, idx=idx)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added \n' 
              '    \'' + add_key + '_pars' + '\', fitted parameters for splicing dynamics (adata.var)')

    if plot_results:  # Plot Parameter Stats
        n_rows, n_cols = len(var_names[:4]), 6
        figsize = [2 * n_cols, 1.5 * n_rows]  # rcParams['figure.figsize']
        fontsize = rcParams['font.size']
        fig, axes = pl.subplots(nrows=n_rows, ncols=6, figsize=figsize)
        pl.subplots_adjust(wspace=0.7, hspace=0.5)
        for i, gene in enumerate(var_names[:4]):
            if t_max is not False: P[i] *= np.array([1 / m[i], 1 / m[i], 1 / m[i], m[i], 1])[:, None]
            ax = axes[i] if n_rows > 1 else axes
            for j, pij in enumerate(P[i]):
                ax[j].plot(pij)
            ax[len(P[i])].plot(L[i])
            if i == 0:
                for j, name in enumerate(['alpha', 'beta', 'gamma', 't_', 'scaling', 'loss']):
                    ax[j].set_title(name, fontsize=fontsize)
    return dm if return_model else adata if copy else None


def align_dynamics(data, t_max=None, dm=None, idx=None, copy=False):
    adata = data.copy() if copy else data
    alpha, beta, gamma, t_, m = read_pars(adata, pars_names=['alpha', 'beta', 'gamma', 't_', 'alignment_scaling'])
    T = adata.layers['fit_t'] if 'fit_t' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    Tau = adata.layers['fit_tau'] if 'fit_tau' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    Tau_ = adata.layers['fit_tau_'] if 'fit_tau_' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    idx = ~ np.isnan(np.sum(T, axis=0)) if idx is None else idx

    m = np.ones(adata.n_vars)
    if t_max is not False:
        dt = compute_dt(T[:, idx])

        dt_min = t_[idx] / dt.shape[0] / 5
        idx_bool = dt > dt_min
        dt *= idx_bool / idx_bool

        dt_mean = np.nanmean(dt, axis=0)
        dt_mean += dt_mean == 0

        t_max = 100 if t_max is None else t_max
        m[idx] = t_max / (dt_mean * len(dt))

    alpha, beta, gamma, T, t_, Tau, Tau_ = alpha / m, beta / m, gamma / m, T * m, t_ * m, Tau * m, Tau_ * m

    write_pars(adata, [alpha, beta, gamma, t_, m], pars_names=['alpha', 'beta', 'gamma', 't_', 'alignment_scaling'])
    adata.layers['fit_t'] = T
    adata.layers['fit_tau'] = Tau
    adata.layers['fit_tau_'] = Tau_

    m = m[idx]
    if dm is not None:
        dm.alpha, dm.beta, dm.gamma, dm.pars[:3] = np.array([dm.alpha, dm.beta, dm.gamma, dm.pars[:3]]) / m[-1]
        dm.t, dm.tau, dm.t_, dm.pars[4] = np.array([dm.t, dm.tau, dm.t_, dm.pars[4]]) * m[-1]

    return adata if copy else (dm, m)


def recover_latent_time(data, copy=False):
    adata = data.copy() if copy else data

    logg.info('computing shared time', r=True)
    from .dynamical_model_utils import root_time, compute_shared_time
    from .terminal_states import terminal_states
    from ..utils import get_connectivities

    if 'iroot' not in adata.uns.keys():
        if 'root_cells' not in adata.obs.keys(): terminal_states(adata)
        adata.uns['iroot'] = get_connectivities(adata, mode='distances').dot(adata.obs['root_cells']).argmax()
    t, t_ = root_time(adata.layers['fit_t'], root=adata.uns['iroot'])
    adata.obs['latent_time'] = compute_shared_time(t)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added \n'
              '    \'latent_time\', shared time (adata.obs)')
    return adata if copy else None
