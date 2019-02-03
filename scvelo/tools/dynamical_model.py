from .. import settings
from .. import logging as logg
from .utils import make_dense
from .dynamical_model_utils import unspliced, spliced, derivatives, derivatives_o

import numpy as np
import matplotlib.pyplot as pl
from scipy.sparse import issparse


class DynamicsRecovery:
    def __init__(self, adata=None, gene=None, u=None, s=None, use_raw=False):
        # extract actual data
        if u is None or s is None:
            _layers = adata[:, gene].layers
            u = _layers['unspliced'] if use_raw or 'Mu' not in _layers.keys() else _layers['Mu']
            s = _layers['spliced'] if use_raw or 'Ms' not in _layers.keys() else _layers['Ms']
        self.s, self.u = make_dense(s).copy(), make_dense(u).copy()
        self.use_raw = use_raw

        # initialize
        self.u0, self.u0_ = None, None
        self.s0, self.s0_ = None, None
        self.ut, self.st = [], []

        self.alpha, self.beta, self.gamma, self.alpha_, self.t_,  = None, None, None, None, None
        self.t, self.tau, self.o = 0, 0, 0

        self.weights, self.u_b, self.s_b = None, None, None
        self.loss, self.dalpha, self.dbeta, self.dgamma = [], [], [], []
        self.m_dalpha, self.m_dbeta, self.m_dgamma = [0], [0], [0]
        self.v_dalpha, self.v_dbeta, self.v_dgamma = [0], [0], [0]

        self.initialize()

    def initialize(self, perc=95):
        u, s = self.u, self.s
        su_norm = s / np.clip(s.max(0), 1e-3, None) + u / np.clip(u.max(0), 1e-3, None)
        weights = su_norm >= np.percentile(su_norm, perc, axis=0)
        x = weights.multiply(s).tocsr() if issparse(weights) else weights * s
        y = weights.multiply(u).tocsr() if issparse(weights) else weights * u
        xx_ = x.multiply(x).sum(0) if issparse(x) else (x ** 2).sum(0)
        xy_ = x.multiply(y).sum(0) if issparse(x) else (x * y).sum(0)
        gamma = xy_ / xx_

        alpha_, beta = 0, 1
        alpha = u[u >= np.percentile(u, 90)].mean() * beta

        self.alpha, self.beta, self.gamma, self.alpha_ = alpha, beta, gamma, alpha_
        self.u0, self.u0_ = 0, alpha / beta
        self.s0, self.s0_ = 0, alpha / gamma
        self.update_state_dependent()

    def update_state_dependent(self):
        def log(x):  # to avoid invalid values for log.
            return np.log(np.clip(x, 1e-6, 1 - 1e-6))

        u, s = self.u, self.s
        alpha, beta, gamma, alpha_ = self.alpha, self.beta, self.gamma, self.alpha_
        u0, s0, u0_, s0_ = self.u0, self.s0, self.u0_, self.s0_

        # assign states (on if above gamma fit, else off)
        o = np.array(beta * u - gamma * s >= 0, dtype=int)

        # vectorize state-dependent vars
        alpha = alpha * o + alpha_ * (1 - o)
        u0 = u0 * o + u0_ * (1 - o)
        s0 = s0 * o + s0_ * (1 - o)

        # tau and t (explicitly given)
        u_inf = alpha / beta
        tau = - 1 / beta * log((u - u_inf) / (u0 - u_inf))
        t_ = np.max(tau * o)
        t = tau * o + (t_ + tau) * (1 - o)

        # estimated data via ODE
        ut = unspliced(tau, u0, alpha, beta)
        st = spliced(tau, s0, u0, alpha, beta, gamma)

        self.ut, self.st = ut, st
        self.t, self.tau, self.o, self.t_ = t.round(2), tau.round(2), o, t_

    def update_vars(self, r=1e-6, inits='ss', state_gd=True, add_noise_scale=0, b1=0.9, b2=0.999, eps=1e-8, n_regions=5):
        state_gd = False if 'learn' in inits else state_gd
        if self.weights is None:
            self.uniform_weighting(n_regions=n_regions, perc=95)
        dalpha, dbeta, dgamma, dalpha_, du0, ds0, du0_, ds0_, loss = self.get_derivatives(state_gd)
        self.loss = np.append(self.loss, loss)
        self.dalpha = np.append(self.dalpha, dalpha)
        self.dbeta = np.append(self.dbeta, dbeta)
        self.dgamma = np.append(self.dgamma, dgamma)

        noise_alpha = np.random.normal(scale=add_noise_scale) if add_noise_scale > 0 else 0
        noise_other = np.random.normal(scale=add_noise_scale * .1) if add_noise_scale > 0 else 0

        # Get first and second gradient moment estimations from previous step:
        m_dalpha, m_dbeta, m_dgamma = self.m_dalpha[-1], self.m_dbeta[-1], self.m_dgamma[-1]
        v_dalpha, v_dbeta, v_dgamma = self.v_dalpha[-1], self.v_dbeta[-1], self.v_dgamma[-1]

        t = len(self.v_dgamma)

        # update 1st gradient moments (and then bias corrected)
        m_dar = (b1 * m_dalpha + (1 - b1) * dalpha)
        m_dbr = (b1 * m_dbeta + (1 - b1) * dbeta)
        m_dgr = (b1 * m_dgamma + (1 - b1) * dgamma)

        m_da = m_dar / (1 - b1 ** t)
        m_db = m_dbr / (1 - b1 ** t)
        m_dg = m_dgr / (1 - b1 ** t)

        # update 2nd gradient moments (and then bias corrected)
        v_dar = (b2 * v_dalpha + (1 - b2) * dalpha ** 2)
        v_dbr = (b2 * v_dbeta + (1 - b2) * dbeta ** 2)
        v_dgr = (b2 * v_dgamma + (1 - b2) * dgamma ** 2)

        v_da = v_dar / (1 - b2 ** t)
        v_db = v_dbr / (1 - b2 ** t)
        v_dg = v_dgr / (1 - b2 ** t)

        # Adam Update of Parameters
        self.alpha -= r * m_da / (np.sqrt(v_da) + eps) + noise_alpha
        self.beta -= r * m_db / (np.sqrt(v_db) + eps) + noise_other
        self.gamma -= r * m_dg / (np.sqrt(v_dg) + eps) + noise_other

        # Save first gradient moments
        self.m_dalpha = np.append(self.m_dalpha, m_dar)
        self.m_dbeta = np.append(self.m_dbeta, m_dar)
        self.m_dgamma = np.append(self.m_dgamma, m_dar)
        # Save second gradient moments
        self.v_dalpha = np.append(self.v_dalpha, v_dar)
        self.v_dbeta = np.append(self.v_dbeta, v_dar)
        self.v_dgamma = np.append(self.v_dgamma, v_dar)

        if inits is 'ss':
            self.u0_ = np.array(self.alpha / self.beta)
            self.s0_ = np.array(self.alpha / self.gamma)
        elif inits is 'cont':
            self.u0_ = np.max(self.ut * self.o)
            self.s0_ = np.max(self.st * self.o)
        elif inits is 'learn':
            self.u0_ -= du0_ * r
            self.s0_ -= ds0_ * r
        elif inits is 'learn_all':
            self.alpha_ -= dalpha_ * r
            self.u0 -= du0 * r
            self.s0 -= ds0 * r
            self.u0_ -= du0_ * r
            self.s0_ -= ds0_ * r

    def get_derivatives(self, state_gd=True):
        alpha, beta, gamma, alpha_ = self.alpha, self.beta, self.gamma, self.alpha_
        t, tau, o = self.t, self.tau, self.o
        u0, s0, u0_, s0_ = self.u0, self.s0, self.u0_, self.s0_
        diffU = np.multiply(np.array(self.ut - self.u), self.weights)  # TODO check if this actually are np arrays
        diffS = np.multiply(np.array(self.st - self.s), self.weights)

        if state_gd:
            dL_dalpha, dL_dbeta, dL_dgamma, dL_dalpha_, dL_du0, dL_ds0, dL_du0_, dL_ds0_ = \
                derivatives_o(alpha, beta, gamma, tau, o, diffU, diffS, alpha_, u0, s0, u0_, s0_)
        else:
            dL_dalpha, dL_dbeta, dL_dgamma, dL_dalpha_, dL_du0, dL_ds0, dL_du0_, dL_ds0_ = \
                derivatives(alpha, beta, gamma, tau, o, diffU, diffS, alpha_, u0, s0, u0_, s0_)

        loss = np.sqrt(np.sum(diffU ** 2 + diffS ** 2) / len(diffU))
        return dL_dalpha, dL_dbeta, dL_dgamma, dL_dalpha_, dL_du0, dL_ds0, dL_du0_, dL_ds0_, loss

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
            weights[region] = n_regions / len(region)
        # set weights accordingly such that each region has an equal overall contribution.
        self.weights = weights * len(u) / np.sum(weights)
        self.u_b, self.s_b = u_b, s_b

    def plot_regions(self):
        u, s, ut, st = self.u, self.s, self.ut, self.st
        u_b, s_b = self.u_b, self.s_b

        pl.figure(dpi=100)
        pl.scatter(s, u, color='grey')
        pl.xlim(0); pl.ylim(0); pl.xlabel('spliced'); pl.ylabel('unspliced')

        for i in range(len(s_b)):
            pl.plot([s_b[i], s_b[i], 0], [0, u_b[i], u_b[i]])

    def plot_derivatives(self):
        from .dynamical_model_utils import du_dalpha_o, du_dbeta_o, ds_dalpha_o, ds_dbeta_o, ds_dgamma_o
        alpha, beta, gamma, = self.alpha, self.beta, self.gamma
        t, tau, o = self.t, self.tau, self.o
        idx = np.argsort(t)
        t = t[idx]
        pl.plot(t, du_dalpha_o(beta, tau, o)[idx], label=r'$\partial u / \partial\alpha$')
        pl.plot(t, .2 * du_dbeta_o(alpha, beta, tau, o)[idx], label=r'$\partial u / \partial \beta$')
        pl.plot(t, ds_dalpha_o(beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \alpha$')
        pl.plot(t, ds_dbeta_o(alpha, beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \beta$')
        pl.plot(t, .2 * ds_dgamma_o(alpha, beta, gamma, tau, o)[idx], label=r'$\partial s / \partial \gamma$')

        pl.legend()
        pl.xlabel('t')


def read_pars(adata, pars_names=['alpha', 'beta', 'gamma', 't_'], key='fit'):
    pars = []
    for name in pars_names:
        pkey = key + '_' + name
        par = adata.var[pkey] if pkey in adata.var.keys() else np.zeros(adata.n_vars) * np.nan
        pars.append(par)
    return pars


def write_pars(adata, pars, pars_names=['alpha', 'beta', 'gamma', 't_'], add_key='fit'):
    for i, name in enumerate(pars_names):
        adata.var[add_key + '_' + name] = pars[i]


def recover_dynamics(data, var_names='all', n_iter=100, learning_rate=1e-3, use_raw=False, add_key='fit', copy=False):
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
    if 'velocity_genes' in adata.var.keys(): idx = idx & adata.var.velocity_genes.values
    idx = np.where(idx)[0]
    var_names = adata.var_names[idx]

    alpha, beta, gamma, t_ = read_pars(adata)
    T = adata.layers['fit_t'] if 'fit_t' in adata.layers.keys() else np.zeros(adata.shape) * np.nan
    L = np.zeros(shape=(n_iter, adata.n_vars)) * np.nan

    for i, gene in enumerate(var_names):
        drec = DynamicsRecovery(adata, gene, use_raw=use_raw)

        for j in range(n_iter):
            drec.update_vars(r=learning_rate)
            drec.update_state_dependent()
        i = idx[i]
        alpha[i], beta[i], gamma[i], t_[i] = drec.alpha, drec.beta, drec.gamma, drec.t_
        T[:, i] = drec.t
        L[:, i] = drec.loss

    write_pars(adata, [alpha, beta, gamma, t_])
    adata.layers['fit_t'] = T
    adata.varm['loss'] = L.T

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added \n' 
              '    \'' + add_key + '_pars' + '\', fitted parameters for splicing dynamics (adata.var)')

    return adata if copy else None
