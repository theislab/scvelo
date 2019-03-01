from .. import settings
from .. import logging as logg
from .utils import make_dense
from .dynamical_model_utils import unspliced, spliced, vectorize, derivatives, find_swichting_time, assign_timepoints, fit_alpha, fit_scaling, linreg, convolve

import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.sparse import issparse


class DynamicsRecovery:
    def __init__(self, adata=None, gene=None, u=None, s=None, use_raw=False, reinitialize=True, deduct_dropouts=True):
        _layers = adata[:, gene].layers
        self.use_raw = use_raw = use_raw or 'Ms' not in _layers.keys()

        # extract actual data
        if u is None or s is None:
            u = _layers['unspliced'] if use_raw else _layers['Mu']
            s = _layers['spliced'] if use_raw else _layers['Ms']

        self.s = s = make_dense(s)
        self.u = u = make_dense(u)

        # set weights (exclude dropouts and extreme outliers)
        self.weights = np.ravel(s < np.percentile(s, 99)) & np.ravel(u < np.percentile(u, 99))
        if deduct_dropouts:
            self.weights &= np.ravel(make_dense(_layers['spliced']) > 0) & np.ravel(make_dense(_layers['unspliced']) > 0)

        # initialize vars
        zeros, zeros3 = np.zeros(adata.n_obs), np.zeros((3, 1))
        self.u0, self.s0, self.u0_, self.s0_, self.t_, self.scaling = None, None, None, None, None, None
        self.ut, self.st, self.t, self.tau, self.o = zeros, zeros, zeros, zeros, zeros

        self.alpha, self.beta, self.gamma, self.alpha_, self.pars = None, None, None, None, None
        self.dpars, self.m_dpars, self.v_dpars, self.loss = zeros3, zeros3, zeros3, []

        # self.u_b, self.s_b = None, None
        if reinitialize or 'fit_alpha' not in adata.var.keys():
            self.initialize()
        else:
            self.load_pars(adata, gene)

    def initialize(self):
        self.scaling = self.u.sum(0) / self.s.sum(0) * 1.3
        u, s = self.u / self.scaling, self.s

        # initialize beta and gamma from extreme quantiles of s
        perc = 95
        weights_s = s >= np.percentile(s, perc, axis=0)
        u_w, s_w = convolve(u, weights_s), convolve(s, weights_s)

        beta, gamma = 1, linreg(u_w, s_w)

        # initialize alpha and switching points from extreme quantiles of u
        weights_u = u >= np.percentile(u, perc, axis=0)
        u_w, s_w = u[weights_u], s[weights_u]

        alpha, u0_, s0_ = u_w.mean(), u_w.mean(), s_w.mean()
        alpha_, u0, s0, = 0, 0, 0

        t, tau, o = assign_timepoints(u, s, alpha, beta, gamma, u0_=u0_, s0_=s0_)

        # update object with initialized vars
        self.alpha, self.beta, self.gamma, self.alpha_ = alpha, beta, gamma, alpha_
        self.u0, self.s0, self.u0_, self.s0_ = u0, s0, u0_, s0_
        self.t, self.tau, self.o, self.t_ = t, tau, o, np.max(tau * o)
        self.pars = np.array([alpha, beta, gamma, self.t_, self.scaling])[:, None]

        self.loss = [self.get_loss()]
        self.update_state_dependent()
        self.update_scaling()

    def load_pars(self, adata, gene):
        idx = np.where(adata.var_names == gene)[0][0] if isinstance(gene, str) else gene
        self.alpha = adata.var['fit_alpha'][idx]
        self.beta = adata.var['fit_beta'][idx]
        self.gamma = adata.var['fit_gamma'][idx]
        self.scaling = adata.var['fit_scaling'][idx]
        self.t_ = adata.var['fit_t_'][idx]
        self.t = adata.layers['fit_t'][:, idx]

        self.u0, self.s0, self.alpha_ = 0, 0, 0
        self.u0_ = unspliced(self.t_, self.u0, self.alpha, self.beta)
        self.s0_ = spliced(self.t_, self.u0, self.s0, self.alpha, self.beta, self.gamma)
        self.update_state_dependent()

    def fit(self, n_iter=100, r=None, method=None, clip_loss=None):
        improved, idx_update = True, np.clip(int(n_iter / 10), 1, None)

        for i in range(n_iter):
            self.update_vars(r=r, method=method, clip_loss=clip_loss)
            if improved or (i % idx_update == 1) or i == n_iter - 1:
                improved = self.update_state_dependent()
                # improved = improved or self.update_scaling()

    def update_state_dependent(self):
        u, s, w = self.u / self.scaling, self.s, self.weights
        u_w, s_w = (u, s) if w is None else (u[w], s[w])
        alpha, beta, gamma, scaling = self.alpha, self.beta, self.gamma, self.scaling

        improved_tau, improved_alpha, improved_scaling = False, False, False
        # find optimal switching (generalized lin.reg) & assign timepoints/states (explicit)
        tau_w, o_w = (self.tau, self.o) if w is None else (self.tau[w], self.o[w])
        t0_ = find_swichting_time(u_w, s_w, tau_w, o_w, alpha, beta, gamma)

        t0_vals = t0_ + np.linspace(-1, 1, num=5) * t0_ / 10
        for t0_ in t0_vals:
            improved_tau = improved_tau or self.update_loss(t_=t0_, reassign_time=True)  # update if improved

        # fit alpha (generalized lin.reg)
        tau_w, o_w = (self.tau, self.o) if w is None else (self.tau[w], self.o[w])
        alpha = fit_alpha(u_w, s_w, tau_w, o_w, beta, gamma)

        alpha_vals = alpha + np.linspace(-1, 1, num=5) * alpha / 10
        for alpha in alpha_vals:
            improved_alpha = improved_alpha or self.update_loss(alpha=alpha, reassign_time=True)  # update if improved

        # fit scaling (generalized lin.reg)
        t_w, tau_w, o_w = (self.t, self.tau, self.o) if w is None else (self.t[w], self.tau[w], self.o[w])
        alpha, t0_ = self.alpha, self.t_
        scaling = fit_scaling(u_w, t_w, t0_, alpha, beta)
        # t0_ = find_swichting_time(u_w / scaling, s_w, tau_w, o_w, alpha, beta, gamma)
        # t, tau, o = assign_timepoints(u / scaling, s, alpha, beta, gamma, t0_)
        improved_scaling = self.update_loss(scaling=scaling * self.scaling, reassign_time=True)  # update if improved

        return improved_tau or improved_alpha or improved_scaling

    def update_scaling(self):
        # fit scaling and update if improved
        u, s, t, tau, o = self.u / self.scaling, self.s, self.t, self.tau, self.o
        alpha, beta, gamma, alpha_ = self.alpha, self.beta, self.gamma, self.alpha_
        u0, s0, u0_, s0_, t_ = self.u0, self.s0, self.u0_, self.s0_, 0 if self.t_ is None else self.t_

        # fit alpha and scaling and update if improved
        alpha = fit_alpha(u, s, tau, o, beta, gamma)
        t0_ = find_swichting_time(u, s, tau, o, alpha, beta, gamma)
        t, tau, o = assign_timepoints(u, s, alpha, beta, gamma, t0_)
        improved_alpha = self.update_loss(t, t0_, alpha=alpha)

        # fit scaling and update if improved
        scaling = fit_scaling(u, t, t_, alpha, beta) * self.scaling
        t0_ = find_swichting_time(u, s, tau, o, alpha, beta, gamma)
        t, tau, o = assign_timepoints(u, s, alpha, beta, gamma, t0_)
        improved_scaling = self.update_loss(t, t0_, scaling=scaling)

        return improved_alpha or improved_scaling

    def update_vars(self, r=None, method=None, clip_loss=None):
        if r is None:
            r = 1e-3 if method is 'adam' else 1e-6
        if clip_loss is None:
            clip_loss = False if method is 'adam' else True
        # if self.weights is None:
        #    self.uniform_weighting(n_regions=5, perc=95)
        t, t_, alpha, beta, gamma, scaling = self.t, self.t_, self.alpha, self.beta, self.gamma, self.scaling
        dalpha, dbeta, dgamma, dalpha_, dtau, dt_ = derivatives(self.u, self.s, t, t_, alpha, beta, gamma, scaling)

        if method is 'adam':
            b1, b2, eps = 0.9, 0.999, 1e-8

            # update 1st and 2nd order gradient moments
            dpars = np.array([dalpha, dbeta, dgamma])
            m_dpars = b1 * self.m_dpars[:, -1] + (1 - b1) * dpars
            v_dpars = b2 * self.v_dpars[:, -1] + (1 - b2) * dpars**2

            self.dpars = np.c_[self.dpars, dpars]
            self.m_dpars = np.c_[self.m_dpars, m_dpars]
            self.v_dpars = np.c_[self.v_dpars, v_dpars]

            # correct for bias
            t = len(self.m_dpars[0])
            m_dpars /= (1 - b1 ** t)
            v_dpars /= (1 - b2 ** t)

            # Adam parameter update
            alpha -= r * m_dpars[0] / (np.sqrt(v_dpars[0]) + eps)
            beta -= r * m_dpars[1] / (np.sqrt(v_dpars[1]) + eps)
            gamma -= r * m_dpars[2] / (np.sqrt(v_dpars[2]) + eps)

        else:
            alpha -= r * dalpha
            beta -= r * dbeta
            gamma -= r * dgamma
            # tau -= r * dtau
            # t_ -= r * dt_
            # t_ = np.max(self.tau * self.o)
            # t = tau * self.o + (tau + t_) * (1 - self.o)

        improved_vars = self.update_loss(alpha=alpha, beta=beta, gamma=gamma, clip_loss=clip_loss)

    def get_vals(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False):
        alpha = self.alpha if alpha is None else alpha
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma
        scaling = self.scaling if scaling is None else scaling
        t = self.t if t is None else t
        u, s, tau, o, w = self.u / scaling, self.s, self.tau, self.o, self.weights

        if reassign_time:
            u_w, s_w, tau_w, o_w = (u, s, tau, o) if w is None else (u[w], s[w], tau[w], o[w])
            t_ = find_swichting_time(u_w, s_w, tau_w, o_w, alpha, beta, gamma) if t_ is None else t_
            t, tau, o = assign_timepoints(u, s, alpha, beta, gamma, t_)
        else:
            t_ = self.t_ if t_ is None else t_
        return t, t_, alpha, beta, gamma, scaling

    def get_ut(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return unspliced(tau, u0, alpha, beta) * scaling

    def get_st(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return spliced(tau, s0, u0, alpha, beta, gamma)

    def get_optimal_switch(self, alpha=None, beta=None, gamma=None):
        u, s, tau, o, w = self.u / self.scaling, self.s, self.tau, self.o, self.weights
        if w is not None:
            u, s, tau, o = (u[w], s[w], tau[w], o[w])
        return find_swichting_time(u, s, tau, o,
                                   self.alpha if alpha is None else alpha,
                                   self.beta if beta is None else beta,
                                   self.gamma if gamma is None else gamma)

    def get_optimal_alpha(self):
        u, s, tau, o, w = self.u / self.scaling, self.s, self.tau, self.o, self.weights
        if w is not None:
            u, s, tau, o = (u[w], s[w], tau[w], o[w])
        return fit_alpha(u, s, tau, o, self.beta, self.gamma)

    def get_optimal_scaling(self):
        u, t, w = self.u / self.scaling, self.t, self.weights
        if w is not None: u, t = (u[w], t[w])
        return fit_scaling(u, t, self.t_, self.alpha, self.beta)

    def get_time_assignment(self, t_=None, alpha=None, beta=None, gamma=None):
        t_ = self.get_optimal_switch(alpha, beta, gamma) if t_ is None else t_
        t, tau, o = assign_timepoints(self.u / self.scaling, self.s, self.alpha if alpha is None else alpha,
                                      self.beta if beta is None else beta, self.gamma if gamma is None else gamma, t_)
        return t, tau, o

    def get_loss(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False):
        alpha = self.alpha if alpha is None else alpha
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma
        scaling = self.scaling if scaling is None else scaling
        t = self.t if t is None else t
        u, s, tau, o, w = self.u, self.s, self.tau, self.o, self.weights
        if w is not None: u, s, t, tau, o = u[w], s[w], t[w], tau[w], o[w]

        if reassign_time:
            t_ = find_swichting_time(u / scaling, s, tau, o, alpha, beta, gamma) if t_ is None else t_
            t, tau, o = assign_timepoints(u / scaling, s, alpha, beta, gamma, t_)
        else:
            t_ = self.t_ if t_ is None else t_

        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)

        udiff = np.array(unspliced(tau, u0, alpha, beta) * scaling - u)
        sdiff = np.array(spliced(tau, s0, u0, alpha, beta, gamma) - s)
        loss = np.sqrt(np.sum(udiff ** 2 + sdiff ** 2) / len(udiff))
        return loss

    def get_likelihood(self):
        u, s, t, t_ = self.u, self.s, self.t, self.t_
        alpha, beta, gamma, scaling = self.alpha, self.beta, self.gamma, self.scaling
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)

        std_u, std_s, = np.std(u), np.std(s)
        corr = np.corrcoef(u, s)[0, 1]

        udiff = np.array(unspliced(tau, u0, alpha, beta) * scaling - u) / std_u
        sdiff = np.array(spliced(tau, s0, u0, alpha, beta, gamma) - s) / std_s

        denom = 2 * np.pi * std_u * std_s * np.sqrt(1 - corr**2)
        nom = -.5 / (1 - corr**2) * (np.sum(udiff ** 2) + np.sum(sdiff ** 2) - 2 * corr * np.sum(udiff * sdiff))
        likelihood = 1 / denom * np.exp(nom)

        return likelihood

    def update_loss(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False,
                    clip_loss=True, report=False):
        vals = [t_, alpha, beta, gamma, scaling]
        vals_prev = [self.t_, self.alpha, self.beta, self.gamma, self.scaling]
        vals_name = ['t_', 'alpha', 'beta', 'gamma', 'scaling']
        new_vals, new_vals_prev, new_vals_name = [], [], []
        loss_prev = self.loss[-1] if len(self.loss) > 0 else 1e6

        for val, val_prev, val_name in zip(vals, vals_prev, vals_name):
            if val is not None:
                new_vals.append(val)
                new_vals_prev.append(val_prev)
                new_vals_name.append(val_name)

        if reassign_time:
            t_ = self.get_optimal_switch(alpha, beta, gamma) if t_ is None else t_
            t, tau, o = self.get_time_assignment(t_, alpha, beta, gamma)

        loss = self.get_loss(t, t_, alpha, beta, gamma, scaling)
        perform_update = not clip_loss or loss < loss_prev

        if perform_update:
            if len(self.loss) > 0 and loss_prev - loss > loss_prev * .01 and report:  # improvement by at least 1%
                print('Update:',
                      ' '.join(map(str, new_vals_name)),
                      ' '.join(map(str, np.round(new_vals_prev, 2))), '-->',
                      ' '.join(map(str, np.round(new_vals, 2))))

                print('    loss:', np.round(loss_prev, 2), '-->', np.round(loss, 2))

            if 't_' in new_vals_name or reassign_time:
                self.t = t
                self.t_ = t_
                self.o = o = np.array(t < t_, dtype=bool)
                self.tau = t * o + (t - t_) * (1 - o)

            if 'alpha' in new_vals_name: self.alpha = alpha
            if 'beta' in new_vals_name: self.beta = beta
            if 'gamma' in new_vals_name: self.gamma = gamma
            if 'scaling' in new_vals_name: self.scaling = scaling

        self.pars = np.c_[self.pars, np.array([self.alpha, self.beta, self.gamma, self.t_, self.scaling])[:, None]]
        self.loss.append(loss if perform_update else loss_prev)

        return perform_update

    def scale_u(self, scaling=1):
        self.scaling *= scaling

    def uniform_weighting(self, n_regions=5, perc=95):
        from numpy import union1d as union
        from numpy import intersect1d as intersect
        u, s = self.u, self.s
        u_b = np.linspace(0, np.percentile(u, perc), n_regions)
        s_b = np.linspace(0, np.percentile(s, perc), n_regions)

        regions, weights = {}, np.ones(len(u))
        for i in range(n_regions):
            if i == 0:
                region = intersect(np.where(u < u_b[i + 1]), np.where(s < s_b[i + 1]))
            elif i < n_regions - 1:
                lower_cut = union(np.where(u > u_b[i]), np.where(s > s_b[i]))
                upper_cut = intersect(np.where(u < u_b[i + 1]), np.where(s < s_b[i + 1]))
                region = intersect(lower_cut, upper_cut)
            else:
                region = union(np.where(u > u_b[i]), np.where(s > s_b[i]))  # lower_cut for last region
            regions[i] = region
            if len(region) > 0:
                weights[region] = n_regions / len(region)
        # set weights accordingly such that each region has an equal overall contribution.
        self.weights = weights * len(u) / np.sum(weights)
        self.u_b, self.s_b = u_b, s_b

    def plot_phase(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False,
                   color=None, colorbar=False, optimal=False):
        multi = 0
        for param in [alpha, beta, gamma]:
            if param is not None and not np.isscalar(param):
                multi += 1

        if multi == 0:
            u, s = self.u, self.s
            t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
            ut = self.get_ut(t, t_, alpha, beta, gamma, scaling)
            st = self.get_st(t, t_, alpha, beta, gamma, scaling)

            idx_sorted = np.argsort(t)
            ut, st, t = ut[idx_sorted], st[idx_sorted], t[idx_sorted]

            o = np.array(t < t_, dtype=int)
            pl.scatter(s[o == 1], u[o == 1], color='lightblue', )
            pl.scatter(s[o == 0], u[o == 0], color='lightgrey')
            color = 'purple' if color is None else color
            linestyle = '-' if not optimal else '--'
            pl.plot(st, ut, color=color, linestyle=linestyle)
            pl.xlabel('s')
            pl.ylabel('u')
            if colorbar:
                ax = pl.gca()
                cax = inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0)
                import matplotlib as mpl
                cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "purple", "red"])
                cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical')
                cb1.set_ticks([0.1, 0.9])
                cb1.ax.set_yticklabels(['-', '+'])
                pl.sca(ax)
        elif multi == 1:
            a = alpha
            b = beta
            g = gamma
            if alpha is not None and not np.isscalar(alpha):
                n = len(alpha)
                loss = []
                for i in range(n):
                    colorbar = True if i == n - 1 else False
                    self.plot_phase(t, t_, alpha[i], b, g, scaling=scaling, reassign_time=reassign_time,
                                    color=[(i+1)/n, 0, 1-(i+1)/n], colorbar=colorbar)
                    loss.append(self.get_loss(t, t_, alpha[i], b, g, scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(t, t_, alpha[opt], b, g, scaling=scaling, reassign_time=reassign_time,
                                color=[0, 1, 0], colorbar=False, optimal=True)
            elif beta is not None and not np.isscalar(beta):
                n = len(beta)
                loss = []
                for i in range(n):
                    colorbar = True if i == n - 1 else False
                    self.plot_phase(t, t_, a, beta[i], g, scaling=scaling, reassign_time=reassign_time,
                                    color=[(i+1)/n, 0, 1-(i+1)/n], colorbar=colorbar)
                    loss.append(self.get_loss(t, t_, a, beta[i], g, scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(t, t_, a, beta[opt], g, scaling=scaling, reassign_time=reassign_time,
                                color=[0, 1, 0], colorbar=False, optimal=True)
            elif gamma is not None and not np.isscalar(gamma):
                n = len(gamma)
                loss = []
                for i in range(n):
                    colorbar = True if i == n - 1 else False
                    self.plot_phase(t, t_, a, b, gamma[i], scaling=scaling, reassign_time=reassign_time,
                                    color=[(i+1)/n, 0, 1-(i+1)/n], colorbar=colorbar)
                    loss.append(self.get_loss(t, t_, a, b, gamma[i], scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(t, t_, a, b, gamma[opt], scaling=scaling, reassign_time=reassign_time,
                                color=[0, 1, 0], colorbar=False, optimal=True)
        elif multi == 2:
            print('Too many varying Values. Only one varying parameter allowed.')


    def plot_regions(self):
        u, s, ut, st = self.u, self.s, self.ut, self.st
        u_b, s_b = self.u_b, self.s_b

        pl.figure(dpi=100)
        pl.scatter(s, u, color='grey')
        pl.xlim(0); pl.ylim(0); pl.xlabel('spliced'); pl.ylabel('unspliced')

        for i in range(len(s_b)):
            pl.plot([s_b[i], s_b[i], 0], [0, u_b[i], u_b[i]])

    def plot_derivatives(self):
        from .dynamical_model_utils import du, ds
        u, s = self.u, self.s
        alpha, beta, gamma = self.alpha, self.beta, self.gamma
        t, tau, o, t_ = self.t, self.tau, self.o, self.t_

        du0 = np.array(du(t_, alpha, beta))[:, None] * (1 - o)[None, :]
        ds0 = np.array(ds(t_, alpha, beta, gamma))[:, None] * (1 - o)[None, :]

        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        dt = np.array(dtau(u, s, alpha, beta, gamma, u0, s0))

        idx = np.argsort(t)
        t = np.sort(t)
        pl.plot(t, du(beta, tau, )[idx][0], label=r'$\partial u / \partial\alpha$')
        pl.plot(t, .2 * du_dbeta_o(alpha, beta, tau, o)[idx], label=r'$\partial u / \partial \beta$')
        pl.plot(t, ds_dalpha_o(beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \alpha$')
        pl.plot(t, ds_dbeta_o(alpha, beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \beta$')
        pl.plot(t, .2 * ds_dgamma_o(alpha, beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \gamma$')

        pl.legend()
        pl.xlabel('t')


def read_pars(adata, pars_names=['alpha', 'beta', 'gamma', 't_', 'scaling'], key='fit'):
    pars = []
    for name in pars_names:
        pkey = key + '_' + name
        par = adata.var[pkey].values if pkey in adata.var.keys() else np.zeros(adata.n_vars) * np.nan
        pars.append(par)
    return pars


def write_pars(adata, pars, pars_names=['alpha', 'beta', 'gamma', 't_', 'scaling'], add_key='fit'):
    for i, name in enumerate(pars_names):
        adata.var[add_key + '_' + name] = pars[i]


def recover_dynamics(data, var_names='all', n_iter=100, learning_rate=None, add_key='fit', t_max=None, use_raw=False,
                     reinitialize=True, return_model=False, plot_results=False, copy=False, **kwargs):
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

    idx = np.ones(adata.n_vars, dtype=bool) if var_names is 'all' else adata.var_names.isin(var_names)
    if 'velocity_genes' in var_names and 'velocity_genes' in adata.var.keys():
        idx = idx & np.array(adata.var.velocity_genes.values)
    idx = np.where(idx)[0]
    var_names = adata.var_names[idx]

    alpha, beta, gamma, t_, scaling = read_pars(adata)
    L, P, T = [], [], adata.layers['fit_t'] if 'fit_t' in adata.layers.keys() else np.zeros(adata.shape) * np.nan

    for i, gene in enumerate(var_names):
        dm = DynamicsRecovery(adata, gene, use_raw=use_raw, reinitialize=reinitialize)
        if n_iter > 1:
            dm.fit(n_iter, learning_rate, **kwargs)

        i = idx[i]
        alpha[i], beta[i], gamma[i], t_[i], scaling[i] = dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling
        T[:, i] = dm.t
        L.append(dm.loss)
        if plot_results and i < 4:
            P.append(dm.pars)

    m = t_max / T.max(0) if t_max is not None else np.ones(adata.n_vars)
    alpha, beta, gamma, T, t_ = alpha / m, beta / m, gamma / m, T * m, t_ * m

    write_pars(adata, [alpha, beta, gamma, t_, scaling])
    adata.layers['fit_t'] = T

    cur_len = adata.varm['loss'].shape[1] if 'loss' in adata.varm.keys() else 2
    max_len = max(np.max([len(l) for l in L]), cur_len)
    loss = np.ones((adata.n_vars, max_len)) * np.nan

    if 'loss' in adata.varm.keys():
        loss[:, :cur_len] = adata.varm['loss']

    loss[idx] = np.vstack([np.concatenate([l, np.ones(max_len-len(l)) * np.nan]) for l in L])
    adata.varm['loss'] = loss

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added \n' 
              '    \'' + add_key + '_pars' + '\', fitted parameters for splicing dynamics (adata.var)')

    if plot_results:  # Plot Parameter Stats
        figsize = [12, 5]  # rcParams['figure.figsize']
        fontsize = rcParams['font.size']
        fig, axes = pl.subplots(nrows=len(var_names[:4]), ncols=6, figsize=figsize)
        pl.subplots_adjust(wspace=0.7, hspace=0.5)
        for i, gene in enumerate(var_names[:4]):
            P[i] *= np.array([1 / m[i], 1 / m[i], 1 / m[i], m[i], 1])[:, None]
            for j, pij in enumerate(P[i]):
                # pij = pij / m[i] if j < 3 else pij * m[i] if j == 3 else pij
                axes[i][j].plot(pij)
            axes[i][len(P[i])].plot(L[i])
            if i == 0:
                for j, name in enumerate(['alpha', 'beta', 'gamma', 't_', 'scaling', 'loss']):
                    axes[i][j].set_title(name, fontsize=fontsize)

    return dm if return_model else adata if copy else None
