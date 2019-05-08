from .utils import make_dense

import matplotlib.pyplot as pl
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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


def convolve(x, weights=None):
    return (weights.multiply(x).tocsr() if issparse(weights) else weights * x) if weights is not None else x


def apply_weight(arrays, w=None):
    return arrays if w is None else [array[w] for array in arrays]


def linreg(u, s):  # linear regression fit
    ss_ = s.multiply(s).sum(0) if issparse(s) else (s ** 2).sum(0)
    us_ = s.multiply(u).sum(0) if issparse(s) else (s * u).sum(0)
    return us_ / ss_


"""Dynamics delineation"""


def unspliced(tau, u0, alpha, beta):
    expu = exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


def spliced(tau, s0, u0, alpha, beta, gamma):
    c = (alpha - u0 * beta) * inv(gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


def unspliced_spliced(tau, s0, u0, alpha, beta, gamma):
    c = (alpha - u0 * beta) * inv(gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    u = u0 * expu + alpha / beta * (1 - expu)
    s = s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)
    return u, s


def tau_u(u, u0, alpha, beta):
    u_ratio = (u - alpha / beta) / (u0 - alpha / beta)
    return - 1 / beta * log(u_ratio)


def tau_inv(u, s=None, u0=None, s0=None, alpha=None, beta=None, gamma=None):
    if s is None:
        c0 = u0 - alpha / beta
        cu = u - alpha / beta

        tau = - 1 / beta * log(cu / c0)
    else:
        beta_ = beta * inv(gamma - beta)
        ceta_ = alpha / gamma - beta_ * (alpha / beta)

        c0 = s0 - beta_ * u0 - ceta_
        cs = s - beta_ * u - ceta_

        tau = - 1 / gamma * log(cs / c0)
    return tau


def find_swichting_time(u, s, tau, o, alpha, beta, gamma):
    off, on = o == 0, o == 1
    if off.sum() > 0:
        u_, s_, tau_ = u[off], s[off], tau[off]

        beta_ = beta * inv(gamma - beta)
        ceta_ = alpha / gamma - beta_ * alpha / beta

        x = - ceta_ * exp(-gamma * tau_)
        y = s_ - beta_ * u_

        exp_t0_ = (y * x).sum() / (x ** 2).sum()
        t0_ = -1 / gamma * log(exp_t0_ + 1) if -1 < exp_t0_ < 0 else np.max(tau[on]) if on.sum() > 0 else np.max(tau)
    else:
        t0_ = np.max(tau)
    return t0_


def compute_state_likelihoods(u, s, alpha, beta, gamma, t0_=None, u0_=None, s0_=None, adjust_variance=False,
                              normalized=True, mode='hard'):
    def nanvar(x, axis=0):
        return np.nanvar(x, axis) if np.isnan(x).any() else np.var(x, axis)

    if t0_ is None:
        t0_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)
    if u0_ is None or s0_ is None:
        u0_, s0_ = (unspliced(t0_, 0, alpha, beta), spliced(t0_, 0, 0, alpha, beta, gamma))

    # compute inverse timepoints
    if mode is 'projection':
        t0 = tau_u(np.min(u[s > 0]), u0_, 0, beta)
        tpoints = np.linspace(0, t0_, num=200)
        tpoints_ = np.linspace(0, t0, num=200)[1:]

        x_obs = np.vstack([u, s]).T
        xt = np.vstack([unspliced(tpoints, 0, alpha, beta), spliced(tpoints, 0, 0, alpha, beta, gamma)]).T
        xt_ = np.vstack([unspliced(tpoints_, u0_, 0, beta), spliced(tpoints_, s0_, u0_, 0, beta, gamma)]).T

        # assign time points (oth. projection onto 'on' and 'off' curve)
        tau, tau_ = np.zeros(len(u)), np.zeros(len(u))
        for i, xi in enumerate(x_obs):
            diffx, diffx_ = ((xt - xi)**2).sum(1), ((xt_ - xi)**2).sum(1)
            tau[i] = tpoints[np.argmin(diffx)]
            tau_[i] = tpoints_[np.argmin(diffx_)]
    else:
        tau = tau_inv(u, s, 0, 0, alpha, beta, gamma)
        tau = np.clip(tau, 0, t0_)

        tau_ = tau_inv(u, s, u0_, s0_, 0, beta, gamma)
        tau_ = np.clip(tau_, 0, np.max(tau_[s > 0]))

    # compute distances from states (induction state, repression state, steady state)
    distu, distu_,  = u - unspliced(tau, 0, alpha, beta), u - unspliced(tau_, u0_, 0, beta)
    dists, dists_,  = s - spliced(tau, 0, 0, alpha, beta, gamma), s - spliced(tau_, u0_, s0_, 0, beta, gamma)

    distu_steady, distu_steady_ = u - alpha / beta, u
    dists_steady, dists_steady_ = s - alpha / gamma, s

    # compute variances of distances
    varu, varu_ = nanvar(distu), nanvar(distu_)
    vars, vars_ = nanvar(dists), nanvar(dists_)

    varu_steady, varu_steady_ = nanvar(distu_steady), nanvar(distu_steady_)
    vars_steady, vars_steady_ = nanvar(dists_steady), nanvar(dists_steady_)

    # compute variance weighted distances
    distx = distu ** 2 / varu + dists ** 2 / vars
    distx_ = distu_ ** 2 / varu_ + dists_ ** 2 / vars_
    distx_steady = distu_steady ** 2 / varu_steady + dists_steady ** 2 / vars_steady
    distx_steady_ = distu_steady_ ** 2 / varu_steady_ + dists_steady_ ** 2 / vars_steady_

    if adjust_variance:  # recompute variance weighted distances
        id_state = np.argmin([distx_, distx, distx_steady_, distx_steady], axis=0)

        on, off, steady, steady_  = (id_state == 1), (id_state == 0), (id_state == 3), (id_state == 2)
        on, off, steady, steady_ = on / on, off / off, steady / steady, steady_ / steady_

        varu, varu_ = np.nanvar(distu * on, 0), np.nanvar(distu_ * off, 0)
        vars, vars_ = np.nanvar(dists * on, 0), np.nanvar(dists_ * off, 0)

        varu_steady, varu_steady_ = np.nanvar(distu_steady * steady, 0), np.nanvar(distu_steady_ * steady_, 0)
        vars_steady, vars_steady_ = np.nanvar(dists_steady * steady, 0), np.nanvar(dists_steady_ * steady_, 0)

        distx = distu ** 2 / varu + dists ** 2 / vars
        distx_ = distu_ ** 2 / varu_ + dists_ ** 2 / vars_
        distx_steady = distu_steady ** 2 / varu_steady + dists_steady ** 2 / vars_steady
        distx_steady_ = distu_steady_ ** 2 / varu_steady_ + dists_steady_ ** 2 / vars_steady_

    if mode is 'soft':
        div = 1 / (2 * np.pi)
        varx = np.sqrt(varu * vars)
        varx_ = np.sqrt(varu_ * vars_)
        varx_steady = np.sqrt(varu_steady * vars_steady)
        varx_steady_ = np.sqrt(varu_steady_ * vars_steady_)

        l = div / varx * np.exp(-.5 * distx)
        l_ = div / varx_ * np.exp(-.5 * distx_)
        l_steady = div / varx_steady * np.exp(-.5 * distx_steady)
        l_steady_ = div / varx_steady_ * np.exp(-.5 * distx_steady_)

        if normalized:
            l, l_, l_steady, l_steady_ = np.array([l, l_, l_steady, l_steady_]) / (l + l_ + l_steady + l_steady_)

        return l, l_, l_steady, l_steady_

    else:
        o = 1 if mode is 'on' else 0 if mode is 'off' else np.argmin([distx_, distx, distx_steady_, distx_steady], axis=0)
        tau = tau * (o == 1) + tau_ * (1 - o)
        t = tau * o + (tau_ + t0_) * (1 - o)

        return t, tau, o


def assign_timepoints(u, s, alpha, beta, gamma, t0_=None, u0_=None, s0_=None, mode='hard'):
    if t0_ is None:
        t0_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)
    if u0_ is None or s0_ is None:
        u0_, s0_ = (unspliced(t0_, 0, alpha, beta), spliced(t0_, 0, 0, alpha, beta, gamma))

    x_obs = np.vstack([u, s]).T

    if mode is 'projection':
        t0 = tau_u(np.min(u[s > 0]), u0_, 0, beta)
        tpoints = np.linspace(0, t0_, num=200)
        tpoints_ = np.linspace(0, t0, num=200)[1:]

        xt = np.vstack([unspliced(tpoints, 0, alpha, beta), spliced(tpoints, 0, 0, alpha, beta, gamma)]).T
        xt_ = np.vstack([unspliced(tpoints_, u0_, 0, beta), spliced(tpoints_, s0_, u0_, 0, beta, gamma)]).T

        # assign time points (oth. projection onto 'on' and 'off' curve)
        tau, tau_ = np.zeros(len(u)), np.zeros(len(u))
        for i, xi in enumerate(x_obs):
            diffx, diffx_ = ((xt - xi)**2).sum(1), ((xt_ - xi)**2).sum(1)
            tau[i] = tpoints[np.argmin(diffx)]
            tau_[i] = tpoints_[np.argmin(diffx_)]

    else:
        tau = tau_inv(u, s, 0, 0, alpha, beta, gamma)
        tau = np.clip(tau, 0, t0_)

        tau_ = tau_inv(u, s, u0_, s0_, 0, beta, gamma)
        tau_ = np.clip(tau_, 0, np.max(tau_[s > 0]))

    # l, l_, l_steady, l_steady_ = compute_state_likelihoods(u, s, alpha, beta, gamma, t0_, u0_, s0_)

    xt = np.vstack([unspliced(tau, 0, alpha, beta), spliced(tau, 0, 0, alpha, beta, gamma)]).T
    xt_ = np.vstack([unspliced(tau_, u0_, 0, beta), spliced(tau_, s0_, u0_, 0, beta, gamma)]).T

    diffx = ((xt - x_obs)**2).sum(1)
    diffx_ = ((xt_ - x_obs)**2).sum(1)

    o = 1 if mode is 'on' else 0 if mode is 'off' else np.argmax([diffx_, diffx], axis=0)
    tau = tau * o + tau_ * (1 - o)
    t = tau * o + (tau_ + t0_) * (1 - o)

    if mode is 'soft':
        var = np.var(s)
        l = np.exp(- .5 * diffx / var) + .1
        l_ = np.exp(- .5 * diffx_ / var) + .1
        o = l / (l + l_)

    return t, tau, o


def fit_alpha(u, s, tau, o, beta, gamma, fit_scaling=False):
    off, on = o == 0, o == 1

    # 'on' state
    expu, exps = exp(-beta * tau[on]), exp(-gamma * tau[on])

    # 'off' state
    t0_ = np.max(tau * o)
    expu_, exps_ = exp(-beta * tau[off]), exp(-gamma * tau[off])
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

    alpha = (c * x).sum() / (c ** 2).sum() if (on.sum() > 0 and off.sum() > 0) else None

    if fit_scaling:  # alternatively compute alpha and scaling simultaneously
        c = np.concatenate([c_gamma, c_gamma_]).T
        x = np.concatenate([s[on], s[off]]).T
        alpha = (c * x).sum() / (c ** 2).sum()

        c = np.concatenate([c_beta, c_beta_]).T
        x = np.concatenate([u[on], u[off]]).T
        scaling = (c * x).sum() / (c ** 2).sum() / alpha  # ~ alpha * z / alpha
        return alpha, scaling

    else:
        return alpha


def fit_scaling(u, t, t_, alpha, beta):
    tau, alpha, u0, _ = vectorize(t, t_, alpha, beta)
    ut = unspliced(tau, u0, alpha, beta)
    return (u * ut).sum() / (ut ** 2).sum()


def vectorize(t, t_, alpha, beta, gamma=None, alpha_=0, u0=0, s0=0):
    o = np.array(t <= t_, dtype=int)
    tau = t * o + (t - t_) * (1 - o)

    u0_ = unspliced(t_, u0, alpha, beta)
    s0_ = spliced(t_, s0, u0, alpha, beta, gamma if gamma is not None else beta / 2)

    # vectorize u0, s0 and alpha
    u0 = u0 * o + u0_ * (1 - o)
    s0 = s0 * o + s0_ * (1 - o)
    alpha = alpha * o + alpha_ * (1 - o)

    return tau, alpha, u0, s0


"""State-independent derivatives"""


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
    o = np.array(t <= t0_, dtype=int)

    du0 = np.array(du(t0_, alpha, beta, u0))[:, None] * (1 - o)[None, :]
    ds0 = np.array(ds(t0_, alpha, beta, gamma, u0, s0))[:, None] * (1 - o)[None, :]

    tau, alpha, u0, s0 = vectorize(t, t0_, alpha, beta, gamma, alpha_, u0, s0)
    dt = np.array(dtau(u, s, alpha, beta, gamma, u0, s0, du0, ds0))

    # state-dependent derivatives:
    du_a, du_b = du(tau, alpha, beta, u0, du0, dt)
    du_a, du_b = du_a * scaling, du_b * scaling

    ds_a, ds_b, ds_c = ds(tau, alpha, beta, gamma, u0, s0, du0, ds0, dt)

    # evaluate derivative of likelihood:
    udiff = np.array(unspliced(tau, u0, alpha, beta) * scaling - u)
    sdiff = np.array(spliced(tau, s0, u0, alpha, beta, gamma) - s)

    if weights is not None:
        udiff = np.multiply(udiff, weights)
        sdiff = np.multiply(sdiff, weights)

    dl_a = (du_a * (1 - o)).dot(udiff) + (ds_a * (1 - o)).dot(sdiff)
    dl_a_ = (du_a * o).dot(udiff) + (ds_a * o).dot(sdiff)

    dl_b = du_b.dot(udiff) + ds_b.dot(sdiff)
    dl_c = ds_c.dot(sdiff)

    dl_tau, dl_t0_ = None, None
    return dl_a, dl_b, dl_c, dl_a_, dl_tau, dl_t0_


"""Base Class for Dynamics Recovery"""


class BaseDynamics:
    def __init__(self, n_obs):
        self.s, self.u, self.use_raw = None, None, None

        zeros, zeros3 = np.zeros(n_obs), np.zeros((3, 1))
        self.u0, self.s0, self.u0_, self.s0_, self.t_, self.scaling = None, None, None, None, None, None
        self.t, self.tau, self.o, self.weights = zeros, zeros, zeros, zeros

        self.alpha, self.beta, self.gamma, self.alpha_, self.pars = None, None, None, None, None
        self.dpars, self.m_dpars, self.v_dpars, self.loss = zeros3, zeros3, zeros3, []

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

    def get_ut(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False, mode=None):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma) if mode is None \
            else (self.get_time_assignment(mode='on')[1], self.alpha, self.u0, self.s0) if mode is 'on' \
            else (self.get_time_assignment(mode='off')[1], self.alpha_, self.u0_, self.s0_)
        return unspliced(tau, u0, alpha, beta)

    def get_st(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, reassign_time=False, mode=None):
        t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma) if mode is None \
            else (self.get_time_assignment(mode='on')[1], self.alpha, self.u0, self.s0) if mode is 'on' \
            else (self.get_time_assignment(mode='off')[1], self.alpha_, self.u0_, self.s0_)
        return spliced(tau, s0, u0, alpha, beta, gamma)

    def get_vt(self, mode='soft'):
        t, t_, alpha, beta, gamma, scaling = self.get_vals()
        _, _, o = self.get_time_assignment(mode='soft')
        u, s, u_, s_ = self.get_ut(mode='on'), self.get_st(mode='on'), self.get_ut(mode='off'), self.get_st(mode='off')
        return (beta * u - gamma * s) * o + (beta * u_ - gamma * s_) * (1 - o)

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

    def get_time_assignment(self, t_=None, alpha=None, beta=None, gamma=None, mode='hard'):
        t, tau, o = assign_timepoints(self.u / self.scaling, self.s,
                                      self.alpha if alpha is None else alpha,
                                      self.beta if beta is None else beta,
                                      self.gamma if gamma is None else gamma,
                                      self.get_optimal_switch(alpha, beta, gamma) if t_ is None else t_, mode=mode)
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

        loss = np.sum(udiff ** 2 + sdiff ** 2) / len(udiff)
        return loss

    def get_likelihood(self):
        u, s, t, t_ = self.u, self.s, self.t, self.t_
        alpha, beta, gamma, scaling = self.alpha, self.beta, self.gamma, self.scaling
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)

        std_u, std_s, = np.std(u), np.std(s)
        corr = np.corrcoef(u, s)[0, 1]

        udiff = np.array(unspliced(tau, u0, alpha, beta) * scaling - u) / std_u
        sdiff = np.array(spliced(tau, s0, u0, alpha, beta, gamma) - s) / std_s

        denom = 2 * np.pi * std_u * std_s * np.sqrt(1 - corr ** 2)
        nom = -.5 / (1 - corr ** 2) * (np.sum(udiff ** 2) + np.sum(sdiff ** 2) - 2 * corr * np.sum(udiff * sdiff)) / (
                    2 * len(udiff))
        likelihood = 1 / np.sqrt(denom) * np.exp(nom)

        return likelihood

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

    def plot_phase(self, adata=None, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None,
                   reassign_time=False, optimal=False, mode='soft', **kwargs):
        from ..plotting.scatter import scatter
        from ..plotting.simulation import show_full_dynamics

        multi = 0
        for param in [alpha, beta, gamma, t_]:
            if param is not None and not np.isscalar(param):
                multi += 1

        if multi == 0:
            t, t_, alpha, beta, gamma, scaling = self.get_vals(t, t_, alpha, beta, gamma, scaling, reassign_time)
            u, s = self.u, self.s
            ut = self.get_ut(t, t_, alpha, beta, gamma) * scaling
            st = self.get_st(t, t_, alpha, beta, gamma)

            idx_sorted = np.argsort(t)
            ut, st, t = ut[idx_sorted], st[idx_sorted], t[idx_sorted]

            args = {"color": self.get_time_assignment(mode=mode)[2], "color_map": 'RdYlGn', "vmin": -.1, "vmax": 1.1}
            args.update(kwargs)

            ax = scatter(adata, x=s, y=u, colorbar=False, show=False, **kwargs)

            linestyle = '-' if not optimal else '--'
            pl.plot(st, ut, color='purple', linestyle=linestyle)
            pl.xlabel('s')
            pl.ylabel('u')
            if False:  # colorbar
                ax = pl.gca()
                cax = inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0)
                import matplotlib as mpl
                cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "purple", "red"])
                cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical')
                cb1.set_ticks([0.1, 0.9])
                cb1.ax.set_yticklabels(['-', '+'])
                pl.sca(ax)
        elif multi == 1:
            a, b, g = alpha, beta, gamma
            if alpha is not None and not np.isscalar(alpha):
                n = len(alpha)
                loss = []
                for i in range(n):
                    self.plot_phase(adata, t, t_, alpha[i], b, g, scaling=scaling, reassign_time=reassign_time, **kwargs)
                    loss.append(self.get_loss(t, t_, alpha[i], b, g, scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(adata, t, t_, alpha[opt], b, g, scaling=scaling, reassign_time=reassign_time,
                                optimal=True, **kwargs)
            elif beta is not None and not np.isscalar(beta):
                n = len(beta)
                loss = []
                for i in range(n):
                    self.plot_phase(t, t_, a, beta[i], g, scaling=scaling, reassign_time=reassign_time,
                                    color=[(i+1)/n, 0, 1-(i+1)/n])
                    loss.append(self.get_loss(t, t_, a, beta[i], g, scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(adata, t, t_, a, beta[opt], g, scaling=scaling, reassign_time=reassign_time,
                                color=[0, 1, 0], optimal=True)
            elif gamma is not None and not np.isscalar(gamma):
                n = len(gamma)
                loss = []
                for i in range(n):
                    self.plot_phase(adata, t, t_, a, b, gamma[i], scaling=scaling, reassign_time=reassign_time,
                                    color=[(i+1)/n, 0, 1-(i+1)/n])
                    loss.append(self.get_loss(t, t_, a, b, gamma[i], scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(t, t_, a, b, gamma[opt], scaling=scaling, reassign_time=reassign_time,
                                color=[0, 1, 0], optimal=True)
            elif t_ is not None and not np.isscalar(t_):
                n = len(t_)
                loss = []
                for i in range(n):
                    self.plot_phase(adata, t, t_[i], a, b, g, scaling=scaling, reassign_time=reassign_time, **kwargs)
                    loss.append(self.get_loss(t, t_[i], a, b, g, scaling=scaling, reassign_time=reassign_time))
                opt = np.argmin(loss)
                self.plot_phase(adata, t, t_[opt], a, b, g, scaling=scaling, reassign_time=reassign_time, optimal=True, **kwargs)
        elif multi == 2:
            print('Too many varying Values. Only one varying parameter allowed.')

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

    def plot_contours(self, xkey='gamma', ykey='alpha', x_sight=[-.9, .9], y_sight=[-.9, .9], num=20, dpi=None, fontsize=8, **kwargs):
        from ..plotting.utils import update_axes
        x_var = eval('self.' + xkey)
        y_var = eval('self.' + ykey)

        x = np.linspace(x_sight[0], x_sight[1], num=num) * x_var + x_var
        y = np.linspace(y_sight[0], y_sight[1], num=num) * y_var + y_var

        f0 = lambda x, y: self.get_loss(**{xkey: x, ykey: y}, reassign_time=False)
        fp = lambda x, y: self.get_loss(**{xkey: x, ykey: y}, reassign_time=True)
        z0, zp = np.zeros((len(x), len(x))), np.zeros((len(x), len(x)))

        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                z0[i, j] = f0(xi, yi)
                zp[i, j] = fp(xi, yi)

        # ix, iy = np.unravel_index(zp.argmin(), zp.shape)
        # gamma_opt, alpha_opt = x[ix],  y[iy]

        x_label = r'$'+'\\'+xkey+'$' if xkey in ['gamma', 'alpha', 'beta'] else xkey
        y_label = r'$' + '\\' + ykey + '$' if ykey in ['gamma', 'alpha', 'beta'] else ykey
        figsize = rcParams['figure.figsize']
        fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(figsize[0], figsize[1] / 2), dpi=dpi)
        ax1.contourf(x, y, np.log1p(zp.T), levels=20, cmap='RdGy_r')
        contours = ax1.contour(x, y, np.log1p(zp.T), 4, colors='k', linewidths=.5)
        ax1.clabel(contours, inline=True, fontsize=fontsize * .75)
        ax1.scatter(x=x_var, y=y_var, s=50, c='purple', zorder=3, **kwargs)
        # ax1.quiver(self.gamma, self.alpha, 0, self.get_optimal_alpha() - self.alpha, color='k', zorder=3,
        #            headlength=4, headwidth=3, headaxislength=3, alpha=.5)
        ax1.set_xlabel(x_label, fontsize=fontsize)
        ax1.set_ylabel(y_label, fontsize=fontsize)
        ax1.set_title('MSE (profiled)', fontsize=fontsize)
        ax1 = update_axes(ax1, fontsize, frameon=True)

        ax2.contourf(x, y, np.log1p(z0.T), levels=20, cmap='RdGy_r')
        contours = ax2.contour(x, y, np.log1p(z0.T), 4, colors='k', linewidths=.5)
        ax2.clabel(contours, inline=True, fontsize=fontsize * .75)
        ax2.scatter(x=x_var, y=y_var, s=50, c='purple', zorder=3, **kwargs)
        ax2.set_xlabel(x_label, fontsize=fontsize)
        ax2.set_ylabel(y_label, fontsize=fontsize)
        ax2.set_title('MSE', fontsize=fontsize)
        ax2 = update_axes(ax2, fontsize, frameon=True)

        ix, iy = np.unravel_index(zp.argmin(), zp.shape)
        x_opt, y_opt = x[ix].round(2), y[ix].round(2)

        return x_opt, y_opt

