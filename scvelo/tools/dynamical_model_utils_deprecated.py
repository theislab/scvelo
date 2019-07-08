# DEPRECATED

from .. import settings
from .. import logging as logg
from ..preprocessing.moments import get_connectivities
from .utils import make_dense, make_unique_list, test_bimodality

import warnings
import matplotlib.pyplot as pl
from matplotlib import rcParams

import numpy as np
exp = np.exp


def log(x, eps=1e-6):  # to avoid invalid values for log.
    return np.log(np.clip(x, eps, 1 - eps))


def inv(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x_inv = 1 / x * (x != 0)
    return x_inv


def unspliced(tau, u0, alpha, beta):
    expu = exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


def spliced(tau, s0, u0, alpha, beta, gamma):
    c = (alpha - u0 * beta) * inv(gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


def mRNA(tau, u0, s0, alpha, beta, gamma):
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    u = u0 * expu + alpha / beta * (1 - expu)
    s = s0 * exps + alpha / gamma * (1 - exps) + (alpha - u0 * beta) * inv(gamma - beta) * (exps - expu)
    return u, s


def vectorize(t, t_, alpha, beta, gamma=None, alpha_=0, u0=0, s0=0, sorted=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        o = np.array(t < t_, dtype=int)
    tau = t * o + (t - t_) * (1 - o)

    u0_ = unspliced(t_, u0, alpha, beta)
    s0_ = spliced(t_, s0, u0, alpha, beta, gamma if gamma is not None else beta / 2)

    # vectorize u0, s0 and alpha
    u0 = u0 * o + u0_ * (1 - o)
    s0 = s0 * o + s0_ * (1 - o)
    alpha = alpha * o + alpha_ * (1 - o)

    if sorted:
        idx = np.argsort(t)
        tau, alpha, u0, s0 = tau[idx], alpha[idx], u0[idx], s0[idx]
    return tau, alpha, u0, s0


def tau_inv(u, s=None, u0=None, s0=None, alpha=None, beta=None, gamma=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inv_u = (gamma >= beta) if gamma is not None else True
        inv_us = np.invert(inv_u)
    any_invu = np.any(inv_u) or s is None
    any_invus = np.any(inv_us) and s is not None

    if any_invus:  # tau_inv(u, s)
        beta_ = beta * inv(gamma - beta)
        xinf = alpha / gamma - beta_ * (alpha / beta)
        tau = - 1 / gamma * log((s - beta_ * u - xinf) / (s0 - beta_ * u0 - xinf))

    if any_invu:  # tau_inv(u)
        uinf = alpha / beta
        tau_u = - 1 / beta * log((u - uinf) / (u0 - uinf))
        tau = tau_u * inv_u + tau * inv_us if any_invus else tau_u
    return tau


def find_swichting_time(u, s, tau, o, alpha, beta, gamma, plot=False):
    off, on = o == 0, o == 1
    t0_ = np.max(tau[on]) if on.sum() > 0 and np.max(tau[on]) > 0 else np.max(tau)

    if off.sum() > 0:
        u_, s_, tau_ = u[off], s[off], tau[off]

        beta_ = beta * inv(gamma - beta)
        ceta_ = alpha / gamma - beta_ * alpha / beta

        x = - ceta_ * exp(-gamma * tau_)
        y = s_ - beta_ * u_

        exp_t0_ = (y * x).sum() / (x ** 2).sum()
        if -1 < exp_t0_ < 0:
            t0_ = -1 / gamma * log(exp_t0_ + 1)
        if plot:
            pl.scatter(x, y)
    return t0_


def fit_alpha(u, s, tau, o, beta, gamma, fit_scaling=False):
    off, on = o == 0, o == 1
    if on.sum() > 0 or off.sum() > 0 or tau[on].min() == 0 or tau[off].min() == 0:
        alpha = None
    else:
        tau_on, tau_off = tau[on], tau[off]

        # 'on' state
        expu, exps = exp(-beta * tau_on), exp(-gamma * tau_on)

        # 'off' state
        t0_ = np.max(tau_on)
        expu_, exps_ = exp(-beta * tau_off), exp(-gamma * tau_off)
        expu0_, exps0_ = exp(-beta * t0_), exp(-gamma * t0_)

        # from unspliced dynamics
        c_beta = 1 / beta * (1 - expu)
        c_beta_ = 1 / beta * (1 - expu0_) * expu_

        # from spliced dynamics
        c_gamma = (1 - exps) / gamma + (exps - expu) * inv(gamma - beta)
        c_gamma_ = ((1 - exps0_) / gamma + (exps0_ - expu0_) * inv(gamma - beta)) * exps_ - (1 - expu0_) * (exps_ - expu_) * inv(gamma - beta)

        # concatenating together
        c = np.concatenate([c_beta, c_gamma, c_beta_, c_gamma_]).T
        x = np.concatenate([u[on], s[on], u[off], s[off]]).T

        alpha = (c * x).sum() / (c ** 2).sum()

        if fit_scaling:  # alternatively compute alpha and scaling simultaneously
            c = np.concatenate([c_gamma, c_gamma_]).T
            x = np.concatenate([s[on], s[off]]).T
            alpha = (c * x).sum() / (c ** 2).sum()

            c = np.concatenate([c_beta, c_beta_]).T
            x = np.concatenate([u[on], u[off]]).T
            scaling = (c * x).sum() / (c ** 2).sum() / alpha  # ~ alpha * z / alpha
            return alpha, scaling

    return alpha


def fit_scaling(u, t, t_, alpha, beta):
    tau, alpha, u0, _ = vectorize(t, t_, alpha, beta)
    ut = unspliced(tau, u0, alpha, beta)
    return (u * ut).sum() / (ut ** 2).sum()


def tau_s(s, s0, u0, alpha, beta, gamma, u=None, tau=None, eps=1e-2):
    if tau is None:
        tau = tau_inv(u, u0=u0, alpha=alpha, beta=beta) if u is not None else 1
    tau_prev, loss, n_iter, max_iter, mixed_states = 1e6, 1e6, 0, 10, np.any(alpha == 0)
    b0 = (alpha - beta * u0) * inv(gamma - beta)
    g0 = s0 - alpha / gamma + b0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        while np.abs(tau - tau_prev).max() > eps and loss > eps and n_iter < max_iter:
            tau_prev, n_iter = tau, n_iter + 1

            expu, exps = b0 * exp(-beta * tau), g0 * exp(-gamma * tau)
            f = exps - expu + alpha / gamma  # >0
            ft = - gamma * exps + beta * expu  # >0 if on else <0
            ftt = gamma ** 2 * exps - beta ** 2 * expu

            # ax^2 + bx + c = 0  <->  1/2 ftt x^2 + ft x + f = s, where x = (tau - tau_prev)
            a, b, c = ftt / 2, ft, f - s
            term = b ** 2 - 4 * a * c
            update = (-b + np.sqrt(term)) / (2 * a)
            if mixed_states:  # linear approx for off-state due of non-injectivity: tau = tau_prev - c / b
                update = np.nan_to_num(update) * (alpha > 0) + (- c / b) * (alpha <= 0)
            tau = np.nan_to_num(tau_prev + update) * (s != 0) if np.any(term > 0) else tau_prev / 10
            loss = np.abs(alpha / gamma + g0 * exp(-gamma * tau) - b0 * exp(-beta * tau) - s).max()

    return np.clip(tau, 0, None)


def assign_timepoints_projection(u, s, alpha, beta, gamma, t0_=None, u0_=None, s0_=None, n_timepoints=300):
    if t0_ is None:
        t0_ = tau_inv(u=u0_, u0=0, alpha=alpha, beta=beta)
    if u0_ is None or s0_ is None:
        u0_, s0_ = (unspliced(t0_, 0, alpha, beta), spliced(t0_, 0, 0, alpha, beta, gamma))

    tpoints = np.linspace(0, t0_, num=n_timepoints)
    tpoints_ = np.linspace(0, tau_inv(np.min(u[s > 0]), u0=u0_, alpha=0, beta=beta), num=n_timepoints)[1:]

    xt = np.vstack([unspliced(tpoints, 0, alpha, beta), spliced(tpoints, 0, 0, alpha, beta, gamma)]).T
    xt_ = np.vstack([unspliced(tpoints_, u0_, 0, beta), spliced(tpoints_, s0_, u0_, 0, beta, gamma)]).T
    x_obs = np.vstack([u, s]).T

    # assign time points (oth. projection onto 'on' and 'off' curve)
    tau, o, diff = np.zeros(len(u)), np.zeros(len(u), dtype=int), np.zeros(len(u))
    tau_alt, diff_alt = np.zeros(len(u)), np.zeros(len(u))
    for i, xi in enumerate(x_obs):
        diffs, diffs_ = np.linalg.norm((xt - xi), axis=1), np.linalg.norm((xt_ - xi), axis=1)
        idx, idx_ = np.argmin(diffs), np.argmin(diffs_)

        o[i] = np.argmin([diffs_[idx_], diffs[idx]])
        tau[i] = [tpoints_[idx_], tpoints[idx]][o[i]]
        diff[i] = [diffs_[idx_], diffs[idx]][o[i]]

        tau_alt[i] = [tpoints_[idx_], tpoints[idx]][1-o[i]]
        diff_alt[i] = [diffs_[idx_], diffs[idx]][1-o[i]]

    t = tau * o + (t0_ + tau) * (1 - o)

    # # remove meaningless jumps (reassign timepoints/states)
    # idx_ord = np.argsort(t)
    # t_ord = t[idx_ord]
    # dt_ord = t_ord - np.insert(t_ord[:-1], 0, 0)
    # dt = dt_ord[np.argsort(idx_ord)]
    # # Poisson with std = sqrt(mean) -> ~99.9% confidence
    # idx = np.where(dt > dt.mean() + 3 * np.sqrt(dt.mean()))[0]
    #
    # if len(idx) > 0:
    #     tvals = t[idx]
    #     idx_jumps = np.where(t * (1 - o) >= np.min(tvals[tvals > t0_]))[0] if np.any(tvals > t0_) else []
    #     idx_jumps_ = np.where(t * o >= np.min(tvals[tvals <= t0_]))[0] if np.any(tvals <= t0_) else []
    #     idx = np.array(np.concatenate([idx_jumps, idx_jumps_]), dtype=int)
    #
    #     # change if alternative is reasonable
    #     change = diff_alt[idx] < np.clip(2 * diff[idx], diff.mean() + 2 * diff.std(), None)
    #     tau[idx] = tau_alt[idx] * change + tau[idx] * (1 - change)
    #     o[idx] = (1 - o[idx]) * change + o[idx] * (1 - change)
    #
    #     t = tau * o + (t0_ + tau) * (1 - o)

    return t, tau, o


"""State-independent derivatives"""


# def du_du0(beta, tau):
#     return exp(-beta * tau)

# def ds_ds0(gamma, tau):
#     return exp(-gamma * tau)

# def ds_du0(beta, gamma, tau):
#     return - beta / (gamma - beta) * (exp(-gamma * tau) - exp(-beta * tau))

# def dus_u0s0(tau, beta, gamma):
#     du_u0 = exp(-beta * tau)
#     ds_s0 = exp(-gamma * tau)
#     ds_u0 = - beta / (gamma - beta) * (ds_s0 - du_u0)
#     return du_u0, ds_s0, ds_u0

# def dus_tau(tau, alpha, beta, gamma, u0=0, s0=0, du0_t0=0, ds0_t0=0):
#     expu, exps, cb, cc = exp(-beta * tau), exp(-gamma * tau), alpha - beta * u0, alpha - gamma * s0
#     du_tau = (cb - du0_t0) * expu
#     ds_tau = (cc - ds0_t0) * exps - cb / (gamma - beta) * (gamma * exps - beta * expu)
#     + du0_t0 * beta / (gamma - beta) * (exps - expu)
#     return du_tau, ds_tau


def dtau(u, s, alpha, beta, gamma, u0, s0, du0=[0, 0, 0], ds0=[0, 0, 0, 0]):
    a, b, g, gb, b0 = alpha, beta, gamma, gamma - beta, beta * inv(gamma - beta)

    cu = s - a/g - b0 * (u - a/b)
    c0 = s0 - a/g - b0 * (u0 - a/b)
    cu += cu == 0
    c0 += c0 == 0
    cu_, c0_ = 1 / cu, 1 / c0

    dtau_a = b0/g * (c0_ - cu_) + 1/g * c0_ * (ds0[0] - b0 * du0[0])
    dtau_b = 1/gb**2 * ((u - a/g) * cu_ - (u0 - a/g) * c0_)

    dtau_c = - a/g * (1/g**2 - 1/gb**2) * (cu_ - c0_) - b0/g/gb * (u*cu_ - u0*c0_)  # + 1/g**2 * np.log(cu/c0)

    return dtau_a, dtau_b, dtau_c


def du(tau, alpha, beta, u0=0, du0=[0, 0, 0], dtau=[0, 0, 0]):
    # du0 is the derivative du0 / d(alpha, beta, tau)
    expu, cb = exp(-beta * tau), alpha / beta
    du_a = du0[0] * expu + 1. / beta * (1 - expu) + (alpha - beta * u0) * dtau[0] * expu
    du_b = du0[1] * expu - cb / beta * (1 - expu) + (cb - u0) * tau * expu + (alpha - beta * u0) * dtau[1] * expu
    return du_a, du_b


def ds(tau, alpha, beta, gamma, u0=0, s0=0, du0=[0, 0, 0], ds0=[0, 0, 0, 0], dtau=[0, 0, 0]):
    # ds0 is the derivative ds0 / d(alpha, beta, gamma, tau)
    expu, exps, = exp(-beta * tau), exp(-gamma * tau)
    expus = exps - expu

    cbu = (alpha - beta * u0) * inv(gamma - beta)
    ccu = (alpha - gamma * u0) * inv(gamma - beta)
    ccs = alpha / gamma - s0 - cbu

    ds_a = ds0[0] * exps + 1. / gamma * (1 - exps) + 1 * inv(gamma - beta) * (1 - beta * du0[0]) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[0]
    ds_b = ds0[1] * exps + cbu * tau * expu + 1 * inv(gamma - beta) * (ccu - beta * du0[1]) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[1]
    ds_c = ds0[2] * exps + ccs * tau * exps - alpha / gamma**2 * (1 - exps) - cbu * inv(gamma - beta) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[2]

    return ds_a, ds_b, ds_c


def derivatives(u, s, t, t0_, alpha, beta, gamma, scaling=1, alpha_=0, u0=0, s0=0, weights=None):
    o = np.array(t < t0_, dtype=int)

    du0 = np.array(du(t0_, alpha, beta, u0))[:, None] * (1 - o)[None, :]
    ds0 = np.array(ds(t0_, alpha, beta, gamma, u0, s0))[:, None] * (1 - o)[None, :]

    tau, alpha, u0, s0 = vectorize(t, t0_, alpha, beta, gamma, alpha_, u0, s0)
    dt = np.array(dtau(u, s, alpha, beta, gamma, u0, s0, du0, ds0))

    # state-dependent derivatives:
    du_a, du_b = du(tau, alpha, beta, u0, du0, dt)
    du_a, du_b = du_a * scaling, du_b * scaling

    ds_a, ds_b, ds_c = ds(tau, alpha, beta, gamma, u0, s0, du0, ds0, dt)

    # evaluate derivative of likelihood:
    ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)

    # udiff = np.array(ut * scaling - u)
    udiff = np.array(ut - u / scaling)
    sdiff = np.array(st - s)

    if weights is not None:
        udiff = np.multiply(udiff, weights)
        sdiff = np.multiply(sdiff, weights)

    dl_a = (du_a * (1 - o)).dot(udiff) + (ds_a * (1 - o)).dot(sdiff)
    dl_a_ = (du_a * o).dot(udiff) + (ds_a * o).dot(sdiff)

    dl_b = du_b.dot(udiff) + ds_b.dot(sdiff)
    dl_c = ds_c.dot(sdiff)

    dl_tau, dl_t0_ = None, None
    return dl_a, dl_b, dl_c, dl_a_, dl_tau, dl_t0_


class BaseDynamics:
    def __init__(self, adata=None, u=None, s=None):
        self.s, self.u = s, u

        zeros, zeros3 = np.zeros(adata.n_obs), np.zeros((3, 1))
        self.u0, self.s0, self.u0_, self.s0_, self.t_, self.scaling = None, None, None, None, None, None
        self.t, self.tau, self.o, self.weights = zeros, zeros, zeros, zeros

        self.alpha, self.beta, self.gamma, self.alpha_, self.pars = None, None, None, None, None
        self.dpars, self.m_dpars, self.v_dpars, self.loss = zeros3, zeros3, zeros3, []

    def uniform_weighting(self, n_regions=5, perc=95):  # deprecated
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

    def plot_regions(self):
        u, s, ut, st = self.u, self.s, self.ut, self.st
        u_b, s_b = self.u_b, self.s_b

        pl.figure(dpi=100)
        pl.scatter(s, u, color='grey')
        pl.xlim(0);
        pl.ylim(0);
        pl.xlabel('spliced');
        pl.ylabel('unspliced')

        for i in range(len(s_b)):
            pl.plot([s_b[i], s_b[i], 0], [0, u_b[i], u_b[i]])

    def plot_derivatives(self):
        u, s = self.u, self.s
        alpha, beta, gamma = self.alpha, self.beta, self.gamma
        t, tau, o, t_ = self.t, self.tau, self.o, self.t_

        du0 = np.array(du(t_, alpha, beta))[:, None] * (1 - o)[None, :]
        ds0 = np.array(ds(t_, alpha, beta, gamma))[:, None] * (1 - o)[None, :]

        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        dt = np.array(dtau(u, s, alpha, beta, gamma, u0, s0))

        du_a, du_b = du(tau, alpha, beta, u0=u0, du0=du0, dtau=dt)
        ds_a, ds_b, ds_c = ds(tau, alpha, beta, gamma, u0=u0, s0=s0, du0=du0, ds0=ds0, dtau=dt)

        idx = np.argsort(t)
        t = np.sort(t)

        pl.plot(t, du_a[idx], label=r'$\partial u / \partial\alpha$')
        pl.plot(t, .2 * du_b[idx], label=r'$\partial u / \partial \beta$')
        pl.plot(t, ds_a[idx], label=r'$\partial s / \partial \alpha$')
        pl.plot(t, ds_b[idx], label=r'$\partial s / \partial \beta$')
        pl.plot(t, .2 * ds_c[idx], label=r'$\partial s / \partial \gamma$')

        pl.legend()
        pl.xlabel('t')


class DynamicsRecovery(BaseDynamics):
    def __init__(self, adata=None, gene=None, u=None, s=None, use_raw=False, load_pars=None, fit_scaling=False,
                 fit_time=True, fit_switching=True, fit_steady_states=True, fit_alpha=True, fit_connected_states=True):
        super(DynamicsRecovery, self).__init__(adata.n_obs)

        _layers = adata[:, gene].layers
        self.gene = gene
        self.use_raw = use_raw = use_raw or 'Ms' not in _layers.keys()

        # extract actual data
        if u is None or s is None:
            u = make_dense(_layers['unspliced']) if use_raw else make_dense(_layers['Mu'])
            s = make_dense(_layers['spliced']) if use_raw else make_dense(_layers['Ms'])
        self.s, self.u = s, u

        # set weights for fitting (exclude dropouts and extreme outliers)
        nonzero = np.ravel(s > 0) & np.ravel(u > 0)
        s_filter = np.ravel(s < np.percentile(s[nonzero], 98))
        u_filter = np.ravel(u < np.percentile(u[nonzero], 98))

        self.weights = s_filter & u_filter & nonzero
        self.fit_scaling = fit_scaling
        self.fit_time = fit_time
        self.fit_alpha = fit_alpha
        self.fit_switching = fit_switching
        self.fit_steady_states = fit_steady_states
        self.connectivities = get_connectivities(adata) if fit_connected_states is True else fit_connected_states

        if load_pars and 'fit_alpha' in adata.var.keys():
            self.load_pars(adata, gene)
        else:
            self.initialize()

    def initialize(self):
        # set weights
        u, s, w = self.u * 1., self.s * 1., self.weights
        u_w, s_w, perc = u[w], s[w], 98

        # initialize scaling
        self.std_u, self.std_s = np.std(u_w), np.std(s_w)
        scaling = self.std_u / self.std_s if isinstance(self.fit_scaling, bool) else self.fit_scaling
        u, u_w = u / scaling, u_w / scaling

        # initialize beta and gamma from extreme quantiles of s
        if True:
            weights_s = s_w >= np.percentile(s_w, perc, axis=0)
        else:
            us_norm = s_w / np.clip(np.max(s_w, axis=0), 1e-3, None) + u_w / np.clip(np.max(u_w, axis=0), 1e-3, None)
            weights_s = us_norm >= np.percentile(us_norm, perc, axis=0)

        beta, gamma = 1, linreg(convolve(u_w, weights_s), convolve(s_w, weights_s))

        u_inf, s_inf = u_w[weights_s].mean(), s_w[weights_s].mean()
        u0_, s0_ = u_inf, s_inf
        alpha = np.mean([s_inf * gamma, u_inf * beta])  # np.mean([s0_ * gamma, u0_ * beta])

        # initialize switching from u quantiles and alpha from s quantiles
        tstat_u, pval_u, means_u = test_bimodality(u_w, kde=True)
        tstat_s, pval_s, means_s = test_bimodality(s_w, kde=True)
        self.pval_steady = max(pval_u, pval_s)
        self.u_steady = means_u[1]
        self.s_steady = means_s[1]

        if self.pval_steady < .1:
            u_inf = np.mean([u_inf, self.u_steady])
            s_inf = np.mean([s_inf, self.s_steady])
            alpha = s_inf * gamma
            beta = alpha / u_inf

            weights_u = u_w >= np.percentile(u_w, perc, axis=0)
            u0_, s0_ = u_w[weights_u].mean(), s_w[weights_u].mean()

        # alpha, beta, gamma = np.array([alpha, beta, gamma]) * scaling
        t_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)

        # update object with initialized vars
        alpha_, u0, s0, = 0, 0, 0
        self.alpha, self.beta, self.gamma, self.alpha_, self.scaling = alpha, beta, gamma, alpha_, scaling
        self.u0, self.s0, self.u0_, self.s0_, self.t_ = u0, s0, u0_, s0_, t_
        self.pars = np.array([alpha, beta, gamma, self.t_, self.scaling])[:, None]

        # initialize time point assignment
        self.t, self.tau, self.o = self.get_time_assignment()
        self.loss = [self.get_loss()]

        self.update_scaling(sight=.5)
        self.update_scaling(sight=.1)

    def load_pars(self, adata, gene):
        idx = np.where(adata.var_names == gene)[0][0] if isinstance(gene, str) else gene
        self.alpha = adata.var['fit_alpha'][idx]
        self.beta = adata.var['fit_beta'][idx]
        self.gamma = adata.var['fit_gamma'][idx]
        self.scaling = adata.var['fit_scaling'][idx]
        self.t_ = adata.var['fit_t_'][idx]
        self.t = adata.layers['fit_t'][:, idx]
        self.o = self.t < self.t_
        self.tau = self.t * self.o + (self.t - self.t_) * (1 - self.o)
        self.pars = np.array([self.alpha, self.beta, self.gamma, self.t_, self.scaling])[:, None]

        self.u0, self.s0, self.alpha_ = 0, 0, 0
        self.u0_ = unspliced(self.t_, self.u0, self.alpha, self.beta)
        self.s0_ = spliced(self.t_, self.u0, self.s0, self.alpha, self.beta, self.gamma)

        self.update_state_dependent()

    def fit(self, max_iter=100, r=None, method=None, clip_loss=None, assignment_mode=None, min_loss=True):
        updated, idx_update = True, np.clip(int(max_iter / 10), 1, None)

        for i in range(max_iter):
            self.update_vars(r=r, method=method, clip_loss=clip_loss)
            if updated or (i % idx_update == 1) or i == max_iter - 1:
                updated = self.update_state_dependent()
            if i > 10 and (i % idx_update == 1):
                loss_prev, loss = np.max(self.loss[-10:]), self.loss[-1]
                if loss_prev - loss < loss_prev * 1e-3:
                    updated = self.shuffle_pars()
                    if not updated: break

        if self.fit_switching:
            self.update_switching()
        if min_loss:
            alpha, beta, gamma, t_, scaling = self.pars[:, np.argmin(self.loss)]
            up = self.update_loss(None, t_, alpha, beta, gamma, scaling, reassign_time=True)
        self.t, self.tau, self.o = self.get_time_assignment(assignment_mode=assignment_mode)

    def update_state_dependent(self):
        updated = False
        if self.fit_alpha:
            updated = self.update_alpha() | updated
        if self.fit_switching:
            updated = self.update_switching() | updated
        return updated

    def update_scaling(self, sight=.5):  # fit scaling and update if improved
        z_vals = self.scaling + np.linspace(-1, 1, num=5) * self.scaling * sight
        for z in z_vals:
            u0_ = self.u0_ * self.scaling / z
            beta = self.beta / self.scaling * z
            self.update_loss(scaling=z, beta=beta, u0_=u0_, s0_=self.s0_, reassign_time=True)

    def update_alpha(self):  # fit alpha (generalized lin.reg), update if improved
        updated = False
        alpha = self.get_optimal_alpha()
        gamma = self.gamma

        alpha_vals = alpha + np.linspace(-1, 1, num=5) * alpha / 30
        gamma_vals = gamma + np.linspace(-1, 1, num=4) * gamma / 30

        for alpha in alpha_vals:
            for gamma in gamma_vals:
                updated = self.update_loss(alpha=alpha, gamma=gamma, reassign_time=True) | updated
        return updated

    def update_switching(self):  # find optimal switching (generalized lin.reg) & assign timepoints/states (explicit)
        updated = False
        #t_ = self.t_
        t_ = self.get_optimal_switch()
        t_vals = t_ + np.linspace(-1, 1, num=3) * t_ / 5
        for t_ in t_vals:
            updated = self.update_loss(t_=t_, reassign_time=True) | updated

        if True:  # self.pval_steady > .1:
            z_vals = 1 + np.linspace(-1, 1, num=4) / 5
            for z in z_vals:
                beta, gamma = self.beta * z, self.gamma * z
                t, tau, o = self.get_time_assignment(beta=beta, gamma=gamma)
                t_ = np.max(t * o)
                if t_ > 0:
                    update = self.update_loss(t_=np.max(t * o), beta=beta, gamma=gamma, reassign_time=True)
                    updated |= update
                    if update:
                        self.update_loss(t_=self.get_optimal_switch(), reassign_time=True)
        return updated

    def update_vars(self, r=None, method=None, clip_loss=None):
        if r is None:
            r = 1e-2 if method is 'adam' else 1e-5
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
            # Parameters are restricted to be positive
            n_alpha = alpha - r * m_dpars[0] / (np.sqrt(v_dpars[0]) + eps)
            alpha = n_alpha if n_alpha > 0 else alpha
            n_beta = beta - r * m_dpars[1] / (np.sqrt(v_dpars[1]) + eps)
            beta = n_beta if n_beta > 0 else beta
            n_gamma = gamma - r * m_dpars[2] / (np.sqrt(v_dpars[2]) + eps)
            gamma = n_gamma if n_gamma > 0 else gamma

        else:
            # Parameters are restricted to be positive
            n_alpha = alpha - r * dalpha
            alpha = n_alpha if n_alpha > 0 else alpha
            n_beta = beta - r * dbeta
            beta = n_beta if n_beta > 0 else beta
            n_gamma = gamma - r * dgamma
            gamma = n_gamma if n_gamma > 0 else gamma

            # tau -= r * dtau
            # t_ -= r * dt_
            # t_ = np.max(self.tau * self.o)
            # t = tau * self.o + (tau + t_) * (1 - self.o)

        updated_vars = self.update_loss(alpha=alpha, beta=beta, gamma=gamma, clip_loss=clip_loss, reassign_time=False)

    def update_loss(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None,
                    reassign_time=False, clip_loss=True, report=False):
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
        if t_ is None:
            t_ = tau_inv(self.u0_ if u0_ is None else u0_,
                         self.s0_ if s0_ is None else s0_, 0, 0,
                         self.alpha if alpha is None else alpha,
                         self.beta if beta is None else beta,
                         self.gamma if gamma is None else gamma) if u0_ is not None else self.t_

        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling)

        if reassign_time:
            # t_ = self.get_optimal_switch(alpha, beta, gamma) if t_ is None else t_
            t, tau, o = self.get_time_assignment(t_, alpha, beta, gamma, scaling)

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
                if reassign_time: self.t = t
                self.t_ = t_
                self.o = o = np.array(self.t < t_, dtype=bool)
                self.tau = self.t * o + (self.t - t_) * (1 - o)

            if u0_ is not None:
                self.u0_ = u0_
                self.s0_ = s0_

            if 'alpha' in new_vals_name:
                self.alpha = alpha
            if 'beta' in new_vals_name:
                self.beta = beta
            if 'gamma' in new_vals_name:
                self.gamma = gamma
            if 'scaling' in new_vals_name:
                self.scaling = scaling
                # self.rescale_invariant()

            self.pars = np.c_[self.pars, np.array([self.alpha, self.beta, self.gamma, self.t_, self.scaling])[:, None]]
            self.loss.append(loss if perform_update else loss_prev)

        return perform_update

    def rescale_invariant(self, z=None):
        z = self.scaling / self.std_u * self.std_s if z is None else z
        self.alpha, self.beta, self.gamma = np.array([self.alpha, self.beta, self.gamma]) * z
        self.t, self.tau, self.t_ = self.t / z, self.tau / z, self.t_ / z

    def shuffle_pars(self, alpha_sight=[-.5, .5], gamma_sight=[-.5, .5], num=5):
        alpha_vals = np.linspace(alpha_sight[0], alpha_sight[1], num=num) * self.alpha + self.alpha
        gamma_vals = np.linspace(gamma_sight[0], gamma_sight[1], num=num) * self.gamma + self.gamma

        x, y = alpha_vals, gamma_vals
        f = lambda x, y: self.get_loss(alpha=x, gamma=y, reassign_time=self.fit_time)
        z = np.zeros((len(x), len(x)))

        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                z[i, j] = f(xi, yi)
        ix, iy = np.unravel_index(z.argmin(), z.shape)
        return self.update_loss(alpha=x[ix], gamma=y[ix], reassign_time=self.fit_time)


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


def recover_dynamics_deprecated(data, var_names='velocity_genes', max_iter=10, learning_rate=None, t_max=None, use_raw=False,
                     fit_scaling=True, fit_time=True, fit_switching=True, fit_steady_states=True, fit_alpha=True,
                     fit_connected_states=True, min_loss=True, assignment_mode=None, load_pars=None, add_key='fit',
                     method='adam', return_model=True, plot_results=False, copy=False, **kwargs):
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

    alpha, beta, gamma, t_, scaling = read_pars(adata)
    idx = []
    L, P, T = [], [], adata.layers['fit_t'] if 'fit_t' in adata.layers.keys() else np.zeros(adata.shape) * np.nan

    progress = logg.ProgressReporter(len(var_names))
    for i, gene in enumerate(var_names):
        dm = DynamicsRecovery(adata, gene, use_raw=use_raw, load_pars=load_pars, fit_time=fit_time, fit_alpha=fit_alpha,
                              fit_switching=fit_switching, fit_scaling=fit_scaling, fit_steady_states=fit_steady_states,
                              fit_connected_states=fit_connected_states)
        if max_iter > 1:
            dm.fit(max_iter, learning_rate, assignment_mode=assignment_mode, min_loss=min_loss, method=method, **kwargs)

        ix = np.where(adata.var_names == gene)[0][0]
        idx.append(ix)

        alpha[ix], beta[ix], gamma[ix], t_[ix], scaling[ix] = dm.pars[:, -1]
        T[:, ix] = dm.t
        L.append(dm.loss)
        if plot_results and i < 4:
            P.append(dm.pars)

        progress.update()
    progress.finish()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T_max = np.nanpercentile(T, 95, axis=0) - np.nanpercentile(T, 5, axis=0)
    m = t_max / T_max if t_max is not None else np.ones(adata.n_vars)
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
        n_rows, n_cols = len(var_names[:4]), 6
        figsize = [2 * n_cols, 1.5 * n_rows]  # rcParams['figure.figsize']
        fontsize = rcParams['font.size']
        fig, axes = pl.subplots(nrows=n_rows, ncols=6, figsize=figsize)
        pl.subplots_adjust(wspace=0.7, hspace=0.5)
        for i, gene in enumerate(var_names[:4]):
            P[i] *= np.array([1 / m[idx[i]], 1 / m[idx[i]], 1 / m[idx[i]], m[idx[i]], 1])[:, None]
            ax = axes[i] if n_rows > 1 else axes
            for j, pij in enumerate(P[i]):
                ax[j].plot(pij)
            ax[len(P[i])].plot(L[i])
            if i == 0:
                for j, name in enumerate(['alpha', 'beta', 'gamma', 't_', 'scaling', 'loss']):
                    ax[j].set_title(name, fontsize=fontsize)

    return dm if return_model else adata if copy else None
