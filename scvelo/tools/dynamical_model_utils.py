from ..preprocessing.moments import get_connectivities
from .utils import make_dense

import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

from scipy.sparse import issparse
import warnings
import numpy as np
exp = np.exp


"""Helper functions"""


def log(x, eps=1e-6):  # to avoid invalid values for log.
    return np.log(np.clip(x, eps, 1 - eps))


def inv(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x_inv = 1 / x * (x != 0)
    return x_inv


def normalize(X, axis=0, min_confidence=None):
    X_sum = np.sum(X, axis=axis)
    if min_confidence: X_sum += min_confidence
    X_sum += X_sum == 0
    return X / X_sum


def convolve(x, weights=None):
    return (weights.multiply(x).tocsr() if issparse(weights) else weights * x) if weights is not None else x


def linreg(u, s):  # linear regression fit
    ss_ = s.multiply(s).sum(0) if issparse(s) else (s ** 2).sum(0)
    us_ = s.multiply(u).sum(0) if issparse(s) else (s * u).sum(0)
    return us_ / ss_


def compute_dt(t, clipped=True, axis=0):
    prepend = np.min(t, axis=axis)[None, :]
    dt = np.diff(np.sort(t, axis=axis), prepend=prepend, axis=axis)
    m_dt = np.max([np.mean(dt, axis=axis), np.max(t, axis=axis) / len(t)], axis=axis)
    m_dt = np.clip(m_dt, 0, None)
    if clipped:  # Poisson upper bound
        ub = m_dt + 3 * np.sqrt(m_dt)
        dt = np.clip(dt, 0, ub)
    return dt


def root_time(t, root=None):
    nans = np.isnan(np.sum(t, axis=0))
    if np.any(nans): t = t[:, ~nans]

    t_root = 0 if root is None else t[root]
    o = np.array(t >= t_root, dtype=int)
    t_after = (t - t_root) * o
    t_origin = np.max(t_after, axis=0)
    t_before = (t + t_origin) * (1 - o)

    t_switch = np.min(t_before, axis=0)
    t_rooted = t_after + t_before
    return t_rooted, t_switch


def compute_shared_time(t, perc=None, norm=True):
    nans = np.isnan(np.sum(t, axis=0))
    if np.any(nans): t = np.array(t[:, ~nans])
    t -= np.min(t)

    tx_list = np.percentile(t, [15, 20, 25, 30, 35] if perc is None else perc, axis=1)
    tx_max = np.max(tx_list, axis=1)
    tx_max += tx_max == 0
    tx_list /= tx_max[:, None]

    mse = []
    for tx in tx_list:
        tx_ = np.sort(tx)
        linx = np.linspace(0, 1, num=len(tx_))
        mse.append(np.sum((tx_ - linx) ** 2))
    idx_best = np.argsort(mse)[:2]

    t_shared = tx_list[idx_best].sum(0)
    if norm: t_shared /= t_shared.max()

    return t_shared


"""Dynamics delineation"""


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


def adjust_increments(tau, tau_=None):
    tau_new = np.array(tau)
    tau_ord = np.sort(tau_new)
    dtau = np.diff(tau_ord, prepend=0)

    if tau_ is None:
        #m_dtau = np.max([np.mean(dtau), np.max(tau) / len(tau), 0])
        #ub = m_dtau + 3 * np.sqrt(m_dtau)  # Poisson with std = sqrt(mean) -> ~99.9% confidence
        ub = 3 * np.percentile(dtau, 99.5, axis=0)
    else:
        tau_new_ = np.array(tau_)
        tau_ord_ = np.sort(tau_new_)
        dtau_ = np.diff(tau_ord_, prepend=0)

        #m_dtau = np.min([m_dtau, np.max([np.mean(dtau_), np.max(tau_) / len(tau_), 0])])
        #ub = m_dtau + 3 * np.sqrt(m_dtau)
        ub = 3 * np.percentile(np.hstack([dtau, dtau_]), 99.5, axis=0)

        idx = np.where(dtau_ > ub)[0]
        for i in idx:
            ti, dti = tau_ord_[i], dtau_[i]  # - ub
            tau_new_[tau_ >= ti] -= dti

    idx = np.where(dtau > ub)[0]
    for i in idx:
        ti, dti = tau_ord[i], dtau[i]  # - ub
        tau_new[tau >= ti] -= dti

    return tau_new if tau_ is None else (tau_new, tau_new_)


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


def assign_tau(u, s, alpha, beta, gamma, t_=None, u0_=None, s0_=None, assignment_mode=None):
    if assignment_mode is 'projection' and beta < gamma:
        x_obs = np.vstack([u, s]).T
        t0 = tau_inv(np.min(u[s > 0]), u0=u0_, alpha=0, beta=beta)

        num = np.clip(int(len(u) / 5), 200, 500)
        tpoints = np.linspace(0, t_, num=num)
        tpoints_ = np.linspace(0, t0, num=num)[1:]

        xt = np.vstack(mRNA(tpoints, 0, 0, alpha, beta, gamma)).T
        xt_ = np.vstack(mRNA(tpoints_, u0_, s0_, 0, beta, gamma)).T

        # assign time points (oth. projection onto 'on' and 'off' curve)
        tau, tau_ = np.zeros(len(u)), np.zeros(len(u))
        for i, xi in enumerate(x_obs):
            diffx, diffx_ = np.sum((xt - xi)**2, 1), np.sum((xt_ - xi)**2, 1)
            tau[i] = tpoints[np.argmin(diffx)]
            tau_[i] = tpoints_[np.argmin(diffx_)]
    else:
        tau = tau_inv(u, s, 0, 0, alpha, beta, gamma)
        tau = np.clip(tau, 0, t_)

        tau_ = tau_inv(u, s, u0_, s0_, 0, beta, gamma)
        tau_ = np.clip(tau_, 0, np.max(tau_[s > 0]))

    return tau, tau_, t_


def compute_divergence(u, s, alpha, beta, gamma, scaling=1, t_=None, u0_=None, s0_=None, tau=None, tau_=None,
                       std_u=1, std_s=1, normalized=False, mode='distance', assignment_mode=None, var_scale=False,
                       kernel_width=None, fit_steady_states=True, connectivities=None, constraint_time_increments=True,
                       reg_time=None, reg_par=None, min_confidence=None, pval_steady=None, steady_u=None, steady_s=None,
                       noise_model='chi', time_connectivities=None, **kwargs):
    """Estimates the divergence (avaiable metrics: distance, mse, likelihood, loglikelihood) of ODE to observations

    Arguments
    ---------
    mode: `'distance'`, `'mse'`, `'likelihood'`, `'loglikelihood'`, `'soft_eval'` (default: `'distance'`)

    """
    # set tau, tau_
    if u0_ is None or s0_ is None:
        u0_, s0_ = mRNA(t_, 0, 0, alpha, beta, gamma)
    if tau is None or tau_ is None or t_ is None:
        tau, tau_, t_ = assign_tau(u, s, alpha, beta, gamma, t_, u0_, s0_, assignment_mode)

    std_u /= scaling

    # adjust increments of tau, tau_ to avoid meaningless jumps
    if constraint_time_increments:
        ut, st = mRNA(tau, 0, 0, alpha, beta, gamma)
        ut_, st_ = mRNA(tau_, u0_, s0_, 0, beta, gamma)

        distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
        dists, dists_ = (s - st) / std_s, (s - st_) / std_s

        res = np.array([distu_ ** 2 + dists_ ** 2, distu ** 2 + dists ** 2])
        if connectivities is not None and connectivities is not False:
            res = np.array([connectivities.dot(r) for r in res]) if res.ndim > 2 else connectivities.dot(res.T).T

        o = np.argmin(res, axis=0)

        off, on = o == 0, o == 1
        if np.any(on) and np.any(off): tau[on], tau_[off] = adjust_increments(tau[on], tau_[off])
        elif np.any(on): tau[on] = adjust_increments(tau[on])
        elif np.any(off): tau_[off] = adjust_increments(tau_[off])

    # compute induction/repression state distances
    ut, st = mRNA(tau, 0, 0, alpha, beta, gamma)
    ut_, st_ = mRNA(tau_, u0_, s0_, 0, beta, gamma)

    distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
    dists, dists_ = (s - st) / std_s, (s - st_) / std_s

    distx = distu ** 2 + dists ** 2
    distx_ = distu_ ** 2 + dists_ ** 2

    res, varx = np.array([distx_, distx]), 1  # default vals;

    if noise_model is 'normal':
        if var_scale:
            o = np.argmin([distx_, distx], axis=0)
            varu = np.nanvar(distu * o + distu_ + (1 - o), axis=0)
            vars = np.nanvar(dists * o + dists_ + (1 - o), axis=0)

            distx = distu ** 2 / varu + dists ** 2 / vars
            distx_ = distu_ ** 2 / varu + dists_ ** 2 / vars

            varx = varu * vars

            std_u *= np.sqrt(varu)
            std_s *= np.sqrt(vars)

    # compute steady state distances
    if fit_steady_states:
        distx_steady = ((u - alpha / beta) / std_u) ** 2 + ((s - alpha / gamma) / std_s) ** 2
        distx_steady_ = (u / std_u) ** 2 + (s / std_s) ** 2

        res = np.array([distx_, distx, distx_steady_, distx_steady])

    if connectivities is not None and connectivities is not False:
        res = np.array([connectivities.dot(r) for r in res]) if res.ndim > 2 else connectivities.dot(res.T).T

    # compute variances
    if noise_model is 'chi':
        if var_scale:
            o = np.argmin([distx_, distx], axis=0)
            dist = distx * o + distx_ * (1 - o)
            sign = np.sign(dists * o + dists_ * (1 - o))
            varx = np.mean(dist, axis=0) - np.mean(sign * np.sqrt(dist), axis=0) ** 2
            if kernel_width is not None: varx *= kernel_width ** 2
            res /= varx
        elif kernel_width is not None:
            res /= kernel_width ** 2

    if reg_time is not None and len(reg_time) == len(distu_):
        o = np.argmin(res, axis=0)
        t_max = (t_ + tau_) * (o == 0)
        t_max /= np.max(t_max, axis=0)
        reg_time /= np.max(reg_time)

        dist_tau = (tau - reg_time[:, None]) ** 2
        dist_tau_ = (tau_ + t_ - reg_time[:, None]) ** 2
        mu_res = np.mean(res, axis=1)
        if reg_par is not None: mu_res *= reg_par

        res[0] += dist_tau_ * mu_res[0]
        res[1] += dist_tau * mu_res[1]
        if fit_steady_states:
            res[2] += dist_tau * mu_res[1]
            res[3] += dist_tau_ * mu_res[0]

    if mode is 'tau':
        res = [tau, tau_]

    elif mode is 'likelihood':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)

    elif mode is 'nll':
        res = np.log(2 * np.pi * np.sqrt(varx)) + .5 * res
        if normalized: res = normalize(res, min_confidence=min_confidence)

    elif mode is 'confidence':
        res = np.array([res[0], res[1]])
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        res = np.median(np.max(res, axis=0) - (np.sum(res, axis=0) - np.max(res, axis=0)), axis=1)

    elif mode is 'soft_eval' or mode is 'soft':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)

        o_, o = res[0], res[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode is 'hardsoft_eval' or mode is 'hardsoft':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        o = np.argmax(res, axis=0)
        o_, o = (o == 0) * res[0], (o == 1) * res[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode is 'hard_eval' or mode is 'hard':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        o = np.argmax(res, axis=0)
        o_, o = o == 0, o == 1
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode is 'soft_state':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        res = res[1] - res[0]

    elif mode is 'hard_state':
        res = np.argmin(res, axis=0)

    elif mode is 'steady_state':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        res = res[2] + res[3]

    elif mode is 'assign_timepoints' or mode is 'time':
        o = np.argmin(res, axis=0)

        tau_ *= (o == 0)
        tau  *= (o == 1)

        if 2 in o: o[o == 2] = 1
        if 3 in o: o[o == 3] = 0

        t = tau * (o == 1) + (tau_ + t_) * (o == 0)
        res = [t, tau, o] if mode is 'assign_timepoints' else t

    elif mode is 'gene_likelihood':
        o = np.argmin(res, axis=0)

        tau_ *= (o == 0)
        tau *= (o == 1)

        if 2 in o: o[o == 2] = 1
        if 3 in o: o[o == 3] = 0

        distu = distu * (o == 1) + distu_ * (o == 0)
        dists = dists * (o == 1) + dists_ * (o == 0)

        if True:
            idx = (u > np.max(u) / 3) & (s > np.max(s) / 3)
            if np.sum(idx) > 0: distu, dists = distu[idx], dists[idx]

        distx = distu ** 2 + dists ** 2
        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        varx = np.mean(distx) - np.mean(np.sign(dists) * np.sqrt(distx)) ** 2
        n = np.clip(len(distu) - len(distu) * .01, 2, None)
        ll = - 1 / 2 / n * np.sum(distx) / varx - 1 / 2 * np.log(2 * np.pi * varx)
        res = np.exp(ll)

    elif mode is 'velocity':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        res = np.array([res[0], res[1], res[2], res[3], 1e-6 * np.ones(res[0].shape)]) if res.ndim > 2 \
            else np.array([res[0], res[1], min_confidence * np.ones(res[0].shape)])
        if normalized: res = normalize(res, min_confidence=min_confidence)
        res = np.argmax(res, axis=0)
        o_, o = res == 0, res == 1
        t = tau * o + (tau_ + t_) * o_
        if time_connectivities:
            if time_connectivities is True: time_connectivities = connectivities
            t = time_connectivities.dot(t)
            o = (res < 2) * (t < t_)
            o_ = (res < 2) * (t >= t_)
            tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
            ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)
            ut_, st_ = ut, st

        ut = ut * o + ut_ * o_
        st = st * o + st_ * o_
        alpha = alpha * o

        vt = (ut * beta - st * gamma)  # ds/dt
        wt = (alpha - beta * ut) * scaling  # du/dt

        vt, wt = np.clip(vt, -s, None), np.clip(wt, -u * scaling, None)

        res = [vt, wt]

    elif mode is 'soft_velocity':
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-.5 * res)
        if normalized: res = normalize(res, min_confidence=min_confidence)
        o_, o = res[0], res[1]
        ut = ut * o + ut_ * o_
        st = st * o + st_ * o_
        alpha = alpha * o
        vt = (ut * beta - st * gamma)  # ds/dt
        wt = (alpha - beta * ut) * scaling  # du/dt
        res = [vt, wt]

    return res


def assign_timepoints(**kwargs):
    return compute_divergence(**kwargs, mode='assign_timepoints')


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


def curve_dists(u, s, alpha, beta, gamma, t_=None, u0_=None, s0_=None, std_u=1, std_s=1, scaling=1, num=None):
    if u0_ is None or s0_ is None:
        u0_, s0_ = mRNA(t_, 0, 0, alpha, beta, gamma)

    x_obs = np.vstack([u, s]).T
    std_x = np.vstack([std_u / scaling, std_s]).T
    t0 = tau_inv(np.min(u[s > 0]), u0=u0_, alpha=0, beta=beta)

    num = np.clip(int(len(u) / 10), 50, 200) if num is None else num
    tpoints = np.linspace(0, t_, num=num)
    tpoints_ = np.linspace(0, t0, num=num)[1:]

    curve_t = np.vstack(mRNA(tpoints, 0, 0, alpha, beta, gamma)).T
    curve_t_ = np.vstack(mRNA(tpoints_, u0_, s0_, 0, beta, gamma)).T

    # match each curve point to nearest observation
    dist, dist_ = np.zeros(len(curve_t)), np.zeros(len(curve_t_))
    for i, ci in enumerate(curve_t):
        dist[i] = np.min(np.sum((x_obs - ci)**2 / std_x**2, 1))
    for i, ci in enumerate(curve_t_):
        dist_[i] = np.min(np.sum((x_obs - ci)**2 / std_x**2, 1))

    return dist, dist_


"""Base Class for Dynamics Recovery"""


class BaseDynamics:
    def __init__(self, adata=None, gene=None, u=None, s=None, use_raw=False, perc=99, max_iter=10, fit_time=True,
                 fit_scaling=True, fit_steady_states=True, fit_connected_states=True, fit_basal_transcription=None,
                 high_pars_resolution=False, steady_state_prior=None):
        self.s, self.u, self.use_raw = None, None, None

        _layers = adata[:, gene].layers
        self.gene = gene
        self.use_raw = use_raw or 'Ms' not in _layers.keys()

        # extract actual data
        if u is None or s is None:
            u = _layers['unspliced'] if self.use_raw else _layers['Mu']
            s = _layers['spliced'] if self.use_raw else _layers['Ms']
        self.s, self.u = make_dense(s), make_dense(u)

        # Basal transcription
        if fit_basal_transcription:
            self.u0, self.s0 = np.min(u), np.min(s)
            self.u -= self.u0
            self.s -= self.s0
        else:
            self.u0, self.s0 = 0, 0

        self.alpha, self.beta, self.gamma, self.scaling, self.t_, self.alpha_ = None, None, None, None, None, None
        self.u0_, self.s0_, self.weights, self.pars = None, None, None, None
        self.t, self.tau, self.o, self.tau_, self.likelihood, self.loss = None, None, None, None, None, None

        self.max_iter = max_iter
        # partition to total of 5 fitting procedures (t_ and alpha, scaling, rates, t_, all together)
        self.simplex_kwargs = {'method': 'Nelder-Mead', 'options': {'maxiter': int(self.max_iter / 5)}}

        self.perc = perc
        self.recoverable = True
        self.initialize_weights()

        self.refit_time = fit_time

        self.assignment_mode = None
        self.steady_state_ratio = None
        self.steady_state_prior = steady_state_prior

        self.fit_scaling = fit_scaling
        self.fit_steady_states = fit_steady_states
        self.fit_connected_states = fit_connected_states
        self.connectivities = get_connectivities(adata) if self.fit_connected_states is True else self.fit_connected_states
        self.high_pars_resolution = high_pars_resolution

    def initialize_weights(self):
        nonzero_s = np.ravel(self.s > 0)
        nonzero_u = np.ravel(self.u > 0)
        filter_by_s = np.sum(nonzero_s) > .1 * len(self.s)
        filter_by_u = np.sum(nonzero_u) > .1 * len(self.u)

        nonzero = np.ones(len(self.s), dtype=bool)
        if filter_by_s: nonzero &= nonzero_s
        if filter_by_u: nonzero &= nonzero_u

        weights = np.array(nonzero, dtype=bool)
        if filter_by_s: weights &= np.ravel(self.s <= np.percentile(self.s[nonzero], self.perc))
        if filter_by_u: weights &= np.ravel(self.u <= np.percentile(self.u[nonzero], self.perc))

        self.weights = weights
        self.std_u = np.std(self.u[weights])
        self.std_s = np.std(self.s[weights])
        self.recoverable = (np.sum(nonzero_s) > 0) & (np.sum(nonzero_u) > 0)

    def load_pars(self, adata, gene):
        idx = np.where(adata.var_names == gene)[0][0] if isinstance(gene, str) else gene
        self.alpha = adata.var['fit_alpha'][idx]
        self.beta = adata.var['fit_beta'][idx] * adata.var['fit_scaling'][idx]
        self.gamma = adata.var['fit_gamma'][idx]
        self.scaling = adata.var['fit_scaling'][idx]
        self.t_ = adata.var['fit_t_'][idx]
        self.steady_state_ratio = self.gamma / self.beta

        self.alpha_ = 0
        self.u0_, self.s0_ = mRNA(self.t_, 0, 0, self.alpha, self.beta, self.gamma)
        self.pars = np.array([self.alpha, self.beta, self.gamma, self.t_, self.scaling])[:, None]

        t = adata.obs['shared_time'] if 'shared_time' in adata.obs.keys() else adata.layers['fit_t'][:, idx]

        if isinstance(self.refit_time, bool):
            self.t, self.tau, self.o = self.get_time_assignment(t=t)
        else:
            self.t = adata.obs[self.refit_time].values if isinstance(self.refit_time, str) else self.refit_time
            self.refit_time = False
            steady_states = t == self.t_
            if np.any(steady_states):
                self.t_ = np.mean(self.t[steady_states])
            self.t, self.tau, self.o = self.get_time_assignment(t=self.t)

        self.loss = [self.get_loss()]

    def get_reads(self, scaling=None, weighted=False):
        scaling = self.scaling if scaling is None else scaling
        u, s = self.u / scaling, self.s
        if weighted:
            u, s = u[self.weights], s[self.weights]
        return u, s

    def get_vars(self, alpha=None, beta=None, gamma=None, scaling=None, t_=None, u0_=None, s0_=None):
        alpha = self.alpha if alpha is None else alpha
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma
        scaling = self.scaling if scaling is None else scaling
        if t_ is None or t_ == 0:
            t_ = self.t_ if u0_ is None else tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)
        return alpha, beta, gamma, scaling, t_

    def get_divergence(self, alpha=None, beta=None, gamma=None, scaling=None, t_=None, u0_=None, s0_=None, mode=None, **kwargs):
        alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
        res = compute_divergence(self.u / scaling, self.s, alpha, beta, gamma, scaling, t_, u0_, s0_, mode=mode,
                                 std_u=self.std_u, std_s=self.std_s, assignment_mode=self.assignment_mode,
                                 connectivities=self.connectivities, fit_steady_states=self.fit_steady_states, **kwargs)
        return res

    def get_time_assignment(self, alpha=None, beta=None, gamma=None, scaling=None, t_=None, u0_=None, s0_=None,
                             t=None, refit_time=None, rescale_factor=None, weighted=False):
        if refit_time is None:
            refit_time = self.refit_time

        if t is not None:
            t_ = self.t_ if t_ is None else t_
            o = np.array(t < t_, dtype=int)
            tau = t * o + (t - t_) * (1 - o)
        elif refit_time:
            if rescale_factor is not None:
                alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)

                u0_ = self.u0_ if u0_ is None else u0_
                s0_ = self.s0_ if s0_ is None else s0_
                rescale_factor *= gamma / beta

                scaling *= rescale_factor
                beta *= rescale_factor

                u0_ /= rescale_factor
                t_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)

            t, tau, o = self.get_divergence(alpha, beta, gamma, scaling, t_, u0_, s0_, mode='assign_timepoints')
            if rescale_factor is not None:
                t *= self.t_ / t_
                tau *= self.t_ / t_
        else:
            t, tau, o = self.t, self.tau, self.o

        if weighted and self.weights is not None:
            t, tau, o = t[self.weights], tau[self.weights], o[self.weights]

        return t, tau, o

    def get_vals(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None,
                 refit_time=None):
        alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
        t, tau, o = self.get_time_assignment(alpha, beta, gamma, scaling, t_, u0_, s0_, t, refit_time)
        return t, t_, alpha, beta, gamma, scaling

    def get_dists(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None,
                  refit_time=None, weighted=True):
        u, s = self.get_reads(scaling, weighted=weighted)

        alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
        t, tau, o = self.get_time_assignment(alpha, beta, gamma, scaling, t_, u0_, s0_, t, refit_time, weighted=weighted)

        if weighted is 'dynamical':
            idx = (u > np.max(u) / 3) & (s > np.max(s) / 3)
            if np.sum(idx) > 0: u, s, t = u[idx], s[idx], t[idx]

        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)

        udiff = np.array(ut - u) / self.std_u * scaling
        sdiff = np.array(st - s) / self.std_s
        reg = (gamma / beta - self.steady_state_ratio) * s / self.std_s if self.steady_state_ratio is not None else 0
        return udiff, sdiff, reg

    def get_residuals(self, **kwargs):
        udiff, sdiff, reg = self.get_dists(**kwargs)
        return np.sign(sdiff) * np.sqrt(udiff ** 2 + sdiff ** 2)

    def get_se(self, **kwargs):
        udiff, sdiff, reg = self.get_dists(**kwargs)
        return np.sum(udiff ** 2 + sdiff ** 2 + reg ** 2)

    def get_mse(self, **kwargs):
        return self.get_se(**kwargs) / np.sum(self.weights)

    def get_loss(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None, refit_time=None):
        return self.get_se(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma, scaling=scaling, u0_=u0_, s0_=s0_, refit_time=refit_time)

    def get_loglikelihood(self, varx=None, **kwargs):
        if 'weighted' not in kwargs: kwargs.update({'weighted': 'dynamical'})
        udiff, sdiff, reg = self.get_dists(**kwargs)
        distx = udiff ** 2 + sdiff ** 2 + reg ** 2
        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        varx = np.mean(distx) - np.mean(np.sign(sdiff) * np.sqrt(distx))**2 if varx is None else varx
        varx += varx == 0  # edge case of mRNAs levels to be the same across all cells
        n = np.clip(len(distx) - len(self.u) * .01, 2, None)
        return - 1 / 2 / n * np.sum(distx) / varx - 1 / 2 * np.log(2 * np.pi * varx)

    def get_likelihood(self, **kwargs):
        if 'weighted' not in kwargs: kwargs.update({'weighted': 'dynamical'})
        likelihood = np.exp(self.get_loglikelihood(**kwargs))
        return likelihood

    def get_curve_likelihood(self):
        alpha, beta, gamma, scaling, t_ = self.get_vars()
        u, s = self.get_reads(scaling, weighted=False)
        varx = self.get_variance()

        dist, dist_ = curve_dists(u, s, alpha, beta, gamma, t_, std_u=self.std_u, std_s=self.std_s, scaling=scaling)
        l = - 1 / 2 / len(dist) * np.sum(dist) / varx - 1 / 2 * np.log(2 * np.pi * varx)
        l_ = - 1 / 2 / len(dist_) * np.sum(dist_) / varx - 1 / 2 * np.log(2 * np.pi * varx)
        likelihood = np.exp(np.max([l, l_]))
        return likelihood

    def get_variance(self, **kwargs):
        if 'weighted' not in kwargs: kwargs.update({'weighted': 'dynamical'})
        udiff, sdiff, reg = self.get_dists(**kwargs)
        distx = udiff ** 2 + sdiff ** 2
        return np.mean(distx) - np.mean(np.sign(sdiff) * np.sqrt(distx)) ** 2

    def get_ut(self, **kwargs):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(**kwargs)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return unspliced(tau, u0, alpha, beta)

    def get_st(self, **kwargs):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(**kwargs)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return spliced(tau, s0, u0, alpha, beta, gamma)

    def get_vt(self, mode='soft_eval'):
        alpha, beta, gamma, scaling, _ = self.get_vars()
        o_, o, ut, st = self.get_divergence(mode=mode)
        return ut * beta - st * gamma

    def plot_phase(self, adata=None, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, refit_time=None,
                   show=True, show_assignments=None, **kwargs):
        from ..plotting.scatter import scatter
        if np.all([x is None for x in [alpha, beta, gamma, scaling, t_]]): refit_time = False
        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, refit_time=refit_time)
        ut = self.get_ut(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma) * scaling
        st = self.get_st(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma)

        idx_sorted = np.argsort(t)
        ut, st, t = ut[idx_sorted], st[idx_sorted], t[idx_sorted]

        args = {"color_map": 'RdYlGn', "vmin": -.1, "vmax": 1.1}
        args.update(kwargs)

        ax = scatter(adata, x=self.s, y=self.u, colorbar=False, show=False, **kwargs)

        pl.plot(st, ut, color='purple', linestyle='-')
        pl.xlabel('spliced'); pl.ylabel('unspliced');

        if show_assignments is not None and show_assignments is not False:
            ax.plot(np.array([self.s[idx_sorted], st]),
                    np.array([self.u[idx_sorted], ut]), color='grey', linewidth=.1 * 1)
        return ax if show is False else None

    def plot_profile_contour(self, xkey='gamma', ykey='alpha', x_sight=.5, y_sight=.5, num=20, contour_levels=4,
                             fontsize=12, refit_time=None, ax=None, color_map='RdGy', figsize=None, dpi=None,
                             vmin=None, vmax=None, horizontal_ylabels=True, show_path=False, show=True,
                             return_color_scale=False, **kwargs):
        from ..plotting.utils import update_axes
        x_var = eval('self.' + xkey)
        y_var = eval('self.' + ykey)

        x = np.linspace(-x_sight, x_sight, num=num) * x_var + x_var
        y = np.linspace(-y_sight, y_sight, num=num) * y_var + y_var

        assignment_mode = self.assignment_mode
        self.assignment_mode = None

        fp = lambda x, y: self.get_likelihood(**{xkey: x, ykey: y}, refit_time=refit_time)

        zp = np.zeros((len(x), len(x)))
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                zp[i, j] = fp(xi, yi)
        log_zp = np.log1p(zp.T)

        if vmin is None: vmin = np.min(log_zp)
        if vmax is None: vmax = np.max(log_zp)

        x_label = r'$'+'\\'+xkey+'$' if xkey in ['gamma', 'alpha', 'beta'] else xkey
        y_label = r'$' + '\\' + ykey + '$' if ykey in ['gamma', 'alpha', 'beta'] else ykey

        if ax is None:
            figsize = rcParams['figure.figsize'] if figsize is None else figsize
            ax = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi).gca()

        ax.contourf(x, y, log_zp, levels=num, cmap=color_map, vmin=vmin, vmax=vmax)
        if contour_levels is not 0:
            contours = ax.contour(x, y, log_zp, levels=contour_levels, colors='k', linewidths=.5)
            fmt = '%1.1f' if np.isscalar(contour_levels) else '%1.0f'
            ax.clabel(contours, fmt=fmt, inline=True, fontsize=fontsize * .75)

        ax.scatter(x=x_var, y=y_var, s=50, c='purple', zorder=3, **kwargs)
        ax.set_xlabel(x_label, fontsize=fontsize)
        ax.set_ylabel(y_label, fontsize=fontsize, rotation=0 if horizontal_ylabels else 90)
        update_axes(ax, fontsize=fontsize, frameon=True)

        if show_path:
            axis = ax.axis()
            x_hist = self.pars[['alpha', 'beta', 'gamma', 't_', 'scaling'].index(xkey)]
            y_hist = self.pars[['alpha', 'beta', 'gamma', 't_', 'scaling'].index(ykey)]
            ax.plot(x_hist, y_hist)
            ax.axis(axis)

        self.assignment_mode = assignment_mode
        if return_color_scale:
            return np.min(log_zp), np.max(log_zp)
        elif not show: return ax

    def plot_profile_hist(self, xkey='gamma', sight=.5, num=20, dpi=None, fontsize=12, ax=None, figsize=None,
                          color_map='RdGy', vmin=None, vmax=None, show=True):
        from ..plotting.utils import update_axes
        x_var = eval('self.' + xkey)

        x = np.linspace(-sight, sight, num=num) * x_var + x_var

        assignment_mode = self.assignment_mode
        self.assignment_mode = None

        fp = lambda x: self.get_likelihood(**{xkey: x}, refit_time=True)

        zp = np.zeros((len(x)))
        for i, xi in enumerate(x):
            zp[i] = fp(xi)

        log_zp = np.log1p(zp.T)
        if vmin is None: vmin = np.min(log_zp)
        if vmax is None: vmax = np.max(log_zp)

        x_label = r'$'+'\\'+xkey+'$' if xkey in ['gamma', 'alpha', 'beta'] else xkey
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        if ax is None:
            fig = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi)
            ax = fig.gca()

        xp = np.linspace(x.min(), x.max(), 1000)
        yp = np.interp(xp, x, log_zp)
        ax.scatter(xp, yp, c=yp, cmap=color_map, edgecolor='none', vmin=vmin, vmax=vmax)
        ax.set_xlabel(x_label, fontsize=fontsize)
        ax.set_ylabel('likelihood', fontsize=fontsize)
        update_axes(ax, fontsize=fontsize, frameon=True)

        self.assignment_mode = assignment_mode
        if not show: return ax

    def plot_profiles(self, params=['alpha', 'beta', 'gamma'], contour_levels=0, sight=.5, num=20, fontsize=12,
                      color_map='RdGy', vmin=None, vmax=None, figsize=None, dpi=None, **kwargs):
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        fig = pl.figure(constrained_layout=True, dpi=dpi, figsize=figsize)
        n = len(params)
        gs = gridspec.GridSpec(n, n, figure=fig)

        for i in range(len(params)):
            for j in range(n-1, i-1, -1):
                xkey = params[j]
                ykey = params[i]
                ax = fig.add_subplot(gs[n - 1 - i, n - 1 - j])
                if xkey == ykey:
                    ax = self.plot_profile_hist(xkey, figsize=figsize, ax=ax, color_map=color_map, num=num,
                                                sight=sight if np.isscalar(sight) else sight[j],
                                                fontsize=fontsize, vmin=vmin, vmax=vmax, show=False)
                    if i == 0 & j == 0:
                        cax = inset_axes(ax, width="7%", height="100%", loc='lower left',
                                         bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0)
                        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=color_map, norm=norm)
                else:
                    vmin_, vmax_ = self.plot_profile_contour(xkey, ykey, figsize=figsize, ax=ax, color_map=color_map,
                                                             contour_levels=contour_levels, vmin=vmin, vmax=vmax,
                                                             x_sight=sight if np.isscalar(sight) else sight[j],
                                                             y_sight=sight if np.isscalar(sight) else sight[i],
                                                             num=num, fontsize=fontsize, return_color_scale=True, **kwargs)
                    if vmin is None or vmax is None: vmin, vmax = vmin_, vmax_  # scaled to first contour plot

                if i != 0:
                    ax.set_xlabel('')
                    ax.set_xticks([])
                if j - n + 1 != 0:
                    ax.set_ylabel('')
                    ax.set_yticks([])

    def plot_state_likelihoods(self, num=300, dpi=None, figsize=None, color_map=None, color_map_steady=None,
                               continuous=True, common_color_scale=True, var_scale=True, kernel_width=None,
                               normalized=None, transitions=None, colorbar=False, alpha_=0.5, linewidths=3,
                               padding_u=.1, padding_s=.1, fontsize=12, title=None, ax=None, **kwargs):
        from ..plotting.utils import update_axes
        from ..plotting.utils import rgb_custom_colormap
        if color_map is None:
            color_map = rgb_custom_colormap(['royalblue', 'white', 'seagreen'], alpha=[1, .5, 1])
        if color_map_steady is None:
            color_map_steady = rgb_custom_colormap(colors=3*['sienna'], alpha=[0, .5, 1])

        alpha, beta, gamma, scaling, t_ = self.get_vars()
        u, s = self.u / scaling, self.s
        padding_u *= np.max(u) - np.min(u)
        padding_s *= np.max(s) - np.min(s)
        uu = np.linspace(np.min(u) - padding_u, np.max(u) + padding_u, num=num)
        ss = np.linspace(np.min(s) - padding_s, np.max(s) + padding_s, num=num)

        grid_u, grid_s = np.meshgrid(uu, ss)
        grid_u = grid_u.flatten()
        grid_s = grid_s.flatten()

        if var_scale: var_scale = self.get_variance()

        dkwargs = {'alpha': alpha, 'beta': beta, 'gamma': gamma, 'scaling': scaling, 't_': t_, 'kernel_width': kernel_width,
                   'std_u': self.std_u, 'std_s': self.std_s, 'var_scale': var_scale, 'normalized': normalized,
                   'fit_steady_states': True, 'constraint_time_increments': False, 'assignment_mode': 'projection'}

        likelihoods = compute_divergence(u, s, mode='soft_state', **dkwargs)
        likelihoods_steady = compute_divergence(u, s, mode='steady_state', **dkwargs)

        likelihoods_grid = compute_divergence(grid_u, grid_s, mode='soft_state', **dkwargs)
        likelihoods_grid_steady = compute_divergence(grid_u, grid_s, mode='steady_state', **dkwargs)

        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        if ax is None:
            fig = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi)
            ax = fig.gca()

        ax.scatter(x=s, y=u, s=50, c=likelihoods_steady, zorder=3, cmap=color_map_steady, edgecolors='black', **kwargs)
        ax.scatter(x=s, y=u, s=50, c=likelihoods, zorder=3, cmap=color_map, edgecolors='black', **kwargs)

        # Grid scatter for test
        # ax.scatter(x=grid_s, y=grid_u, s=50, c=likelihoods, zorder=3, cmap=color_map,
        #            edgecolors='black', vmax=vmax, **kwargs)
        # ax.scatter(x=grid_s, y=grid_u, s=50, c=likelihoods_steady, zorder=3, cmap=color_map_steady,
        #            edgecolors='black', vmax=vmax, **kwargs)

        l_grid, l_grid_steady = likelihoods_grid.reshape(num, num).T, likelihoods_grid_steady.reshape(num, num).T

        if common_color_scale:
            vmax = vmax_steady = np.max([np.abs(likelihoods_grid), np.abs(likelihoods_grid_steady)])
        else:
            vmax, vmax_steady = np.max(np.abs(likelihoods_grid)), None

        if continuous:
            contf_steady = ax.imshow(l_grid_steady, cmap=color_map_steady, alpha=alpha_, vmin=0, vmax=vmax_steady,
                                     aspect='auto', origin='lower', extent=(min(ss), max(ss), min(uu), max(uu)))
            contf = ax.imshow(l_grid, cmap=color_map, alpha=alpha_, vmin=-vmax, vmax=vmax,
                              aspect='auto', origin='lower', extent=(min(ss), max(ss), min(uu), max(uu)))
        else:
            contf_steady = ax.contourf(ss, uu, l_grid_steady, vmin=0, vmax=vmax_steady, levels=30, cmap=color_map_steady)
            contf = ax.contourf(ss, uu, l_grid, vmin=-vmax, vmax=vmax, levels=30, cmap=color_map)

        # Contour lines
        if transitions is not None:
            trans_width = np.max(likelihoods_grid) - np.min(likelihoods_grid)
            transitions = np.multiply(np.array(transitions), [np.min(likelihoods_grid), np.max(likelihoods_grid)])# trans_width
            ax.contour(ss, uu, likelihoods_grid.reshape(num, num).T, transitions, linestyles='solid', colors='k', linewidths=linewidths)

        if colorbar:
            pl.colorbar(contf, ax=ax)
            pl.colorbar(contf_steady, ax=ax)
        ax.set_xlabel('spliced', fontsize=fontsize)
        ax.set_ylabel('unspliced', fontsize=fontsize)
        title = '' if title is None else title
        ax.set_title(title, fontsize=fontsize)
        update_axes(ax, fontsize=fontsize, frameon=True)

        return ax


def get_reads(adata, key='fit', scaled=True, use_raw=False):
    if 'Ms' not in adata.layers.keys(): use_raw = True
    s = make_dense(adata.layers['spliced' if use_raw else 'Ms'])
    u = make_dense(adata.layers['unspliced'if use_raw else 'Mu'])
    if scaled: u /= adata.var[key + '_scaling'].values
    return u, s


def get_vars(adata, scaled=True, key='fit'):
    alpha = adata.var[key + '_alpha'].values if key + '_alpha' in adata.var.keys() else 1
    beta = adata.var[key + '_beta'].values if key + '_beta' in adata.var.keys() else 1
    gamma = adata.var[key + '_gamma'].values
    scaling = adata.var[key + '_scaling'].values if key + '_scaling' in adata.var.keys() else 1
    t_ = adata.var[key + '_t_'].values
    return alpha, beta * scaling if scaled else beta, gamma, scaling, t_


def get_latent_vars(adata, scaled=True, key='fit'):
    scaling = adata.var[key + '_scaling'].values
    std_u = adata.var[key + '_std_u'].values
    std_s = adata.var[key + '_std_s'].values
    u0 = adata.var[key + '_u0'].values
    s0 = adata.var[key + '_s0'].values
    pval_steady = adata.var[key + '_pval_steady'].values if key + '_pval_steady' in adata.var.keys() else None
    steady_u = adata.var[key + '_steady_u'].values if key + '_steady_u' in adata.var.keys() else None
    steady_s = adata.var[key + '_steady_s'].values if key + '_steady_s' in adata.var.keys() else None
    return std_u, std_s, u0 / scaling if scaled else u0, s0, pval_steady, steady_u, steady_s


def get_divergence(adata, mode='soft', use_latent_time=None, use_connectivities=None, **kwargs):
    vdata = adata[:, ~np.isnan(adata.var['fit_alpha'].values)].copy()
    alpha, beta, gamma, scaling, t_ = get_vars(vdata)
    std_u, std_s, u0, s0, pval_steady, steady_u, steady_s = get_latent_vars(vdata)

    kwargs_ = {'kernel_width': None, 'normalized': True, 'var_scale': True, 'reg_par': None, 'min_confidence': 1e-2,
               'constraint_time_increments': False, 'fit_steady_states': True, 'fit_basal_transcription': None,
               'std_u': std_u, 'std_s': std_s, 'pval_steady': pval_steady, 'steady_u': steady_u, 'steady_s': steady_s}
    kwargs_.update(adata.uns['recover_dynamics'])
    kwargs_.update(**kwargs)

    reg_time = None
    if use_latent_time is True: use_latent_time = 'latent_time'
    if isinstance(use_latent_time, str) and use_latent_time in adata.obs.keys():
        reg_time = adata.obs[use_latent_time].values
    u, s = get_reads(vdata, use_raw=kwargs_['use_raw'])
    if kwargs_['fit_basal_transcription']: u, s = u - u0, s - s0
    tau = np.array(vdata.layers['fit_tau']) if 'fit_tau' in vdata.layers.keys() else None
    tau_ = np.array(vdata.layers['fit_tau_']) if 'fit_tau_' in vdata.layers.keys() else None

    res = compute_divergence(u, s, alpha, beta, gamma, scaling, t_, tau=tau, tau_=tau_, reg_time=reg_time, mode=mode,
                             connectivities=get_connectivities(adata) if use_connectivities else None, **kwargs_)
    return res
