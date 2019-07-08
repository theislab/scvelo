from ..preprocessing.moments import get_connectivities
from .utils import make_dense

import matplotlib.pyplot as pl
from matplotlib import rcParams

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


def normalize(X, axis=0):
    X_sum = np.sum(X, axis=axis)
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


def root_time(t, root=0):
    t_root = t[root]
    o = t >= t_root  # True if after 'root'
    t_after = (t - t_root) * o
    t_origin = np.max(t_after, axis=0)
    t_before = (t + t_origin) * (1 - o)

    t_switch = np.min(t_before, axis=0)
    t_rooted = t_after + t_before
    return t_rooted, t_switch


def compute_shared_time(t, perc=None):
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
    t_shared /= t_shared.max()

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


def adjust_increments(tau):
    tau_new = tau * 1.
    tau_ord = np.sort(tau_new)
    dtau = np.diff(tau_ord, prepend=0)
    m_dtau = np.max([np.mean(dtau), np.max(tau) / len(tau), 0])

    # Poisson with std = sqrt(mean) -> ~99.9% confidence
    ub = m_dtau + 6 * np.sqrt(m_dtau)
    idx = np.where(dtau > ub)[0]

    for i in idx:
        ti, dti = tau_ord[i], dtau[i]  # - ub
        tau_new[tau >= ti] -= dti

    return tau_new


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
            diffx, diffx_ = ((xt - xi)**2).sum(1), ((xt_ - xi)**2).sum(1)
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
                       fit_steady_states=True, connectivities=None, constraint_time_increments=True):
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

    # adjust increments of tau, tau_ to avoid meaningless jumps
    if constraint_time_increments:
        ut, st = mRNA(tau, 0, 0, alpha, beta, gamma)
        ut_, st_ = mRNA(tau_, u0_, s0_, 0, beta, gamma)

        std_u /= scaling
        distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
        dists, dists_ = (s - st) / std_s, (s - st_) / std_s

        o = np.argmin(np.array([distu_ ** 2 + dists_ ** 2, distu ** 2 + dists ** 2]), axis=0)

        off, on = o == 0, o == 1
        if np.any(on): tau[on] = adjust_increments(tau[on])
        if np.any(off): tau_[off] = adjust_increments(tau_[off])

    # compute distances from states (induction state, repression state, steady state)
    ut, st = mRNA(tau, 0, 0, alpha, beta, gamma)
    ut_, st_ = mRNA(tau_, u0_, s0_, 0, beta, gamma)

    distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
    dists, dists_ = (s - st) / std_s, (s - st_) / std_s

    distx = distu ** 2 + dists ** 2
    distx_ = distu_ ** 2 + dists_ ** 2

    #distx_ *= 10

    if var_scale:
        o = np.argmin([distx_, distx], axis=0)
        varu = np.nanvar(distu * o + distu_ + (1 - o))
        vars = np.nanvar(dists * o + dists_ + (1 - o))

        distx = distu ** 2 / varu  + dists ** 2 / vars
        distx_ = distu_ ** 2 / varu + dists_ ** 2 / vars
    else:
        varu, vars = 1, 1

    if fit_steady_states:
        distx_steady = (u - alpha / beta) ** 2 / varu + (s - alpha / gamma) ** 2 / vars
        distx_steady_ = u ** 2 / varu + s ** 2 / vars
        res = np.array([distx_, distx, distx_steady_, distx_steady])
    else:
        res = np.array([distx_, distx])

    if connectivities is not None and connectivities is not False:
        if res.ndim > 2:
            res = np.array([connectivities.dot(r) for r in res])
        else:
            res = connectivities.dot(res.T).T

    if mode is 'tau':
        res = [tau, tau_]

    elif mode is 'likelihood':
        res = 1 / (2 * np.pi * np.sqrt(varu * vars)) * np.exp(-.5 * res)
        if normalized: res = normalize(res)

    elif mode is 'nll':
        res = np.log(2 * np.pi * np.sqrt(varu * vars)) + .5 * res
        if normalized: res = normalize(res)

    elif mode is 'soft_eval':
        res = normalize(1 / (2 * np.pi * np.sqrt(varu * vars)) * np.exp(-.5 * res))
        o_, o = res[0], res[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode is 'hard_eval':
        o = np.argmin(res, axis=0)
        o_, o = o[0], o[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode is 'assign_timepoints':
        o = np.argmin(res, axis=0)
        if fit_steady_states:
            idx_steady = (distx_steady < 1)
            tau_[idx_steady], o[idx_steady] = 0, 0

        tau_ *= (o == 0)
        tau  *= (o == 1)

        if 2 in o: o[o == 2] = 1
        if 3 in o: o[o == 3] = 0

        t = tau * (o == 1) + (tau_ + t_) * (o == 0)
        res = [t, tau, o]

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


"""Base Class for Dynamics Recovery"""


class BaseDynamics:
    def __init__(self, adata=None, gene=None, u=None, s=None, use_raw=False, perc=99, max_iter=10, fit_time=True,
                 fit_scaling=True, fit_steady_states=True, fit_connected_states=True):
        self.s, self.u, self.use_raw = None, None, None

        _layers = adata[:, gene].layers
        self.gene = gene
        self.use_raw = use_raw or 'Ms' not in _layers.keys()

        # extract actual data
        if u is None or s is None:
            u = _layers['unspliced'] if self.use_raw else _layers['Mu']
            s = _layers['spliced'] if self.use_raw else _layers['Ms']
        self.s, self.u = make_dense(s), make_dense(u)

        self.alpha, self.beta, self.gamma, self.scaling, self.t_, self.alpha_ = None, None, None, None, None, None
        self.u0, self.s0, self.u0_, self.s0_, self.weights, self.pars = None, None, None, None, None, None
        self.t, self.tau, self.o, self.tau_, self.likelihood, self.loss = None, None, None, None, None, None

        self.max_iter = max_iter
        self.simplex_kwargs = {'method': 'Nelder-Mead', 'options': {'maxiter': self.max_iter}}

        self.perc = perc
        self.initialize_weights()

        self.refit_time = fit_time

        self.assignment_mode = None
        self.steady_state_ratio = None

        self.fit_scaling = fit_scaling
        self.fit_steady_states = fit_steady_states
        self.fit_connected_states = fit_connected_states
        self.connectivities = get_connectivities(adata) if self.fit_connected_states is True else self.fit_connected_states

    def initialize_weights(self):
        nonzero = np.ravel(self.s > 0) & np.ravel(self.u > 0)
        s_filter = np.ravel(self.s < np.percentile(self.s[nonzero], self.perc))
        u_filter = np.ravel(self.u < np.percentile(self.u[nonzero], self.perc))

        self.weights = w = s_filter & u_filter & nonzero
        self.std_u = np.std(self.u[w])
        self.std_s = np.std(self.s[w])

    def load_pars(self, adata, gene):
        idx = np.where(adata.var_names == gene)[0][0] if isinstance(gene, str) else gene
        self.alpha = adata.var['fit_alpha'][idx]
        self.beta = adata.var['fit_beta'][idx]
        self.gamma = adata.var['fit_gamma'][idx]
        self.scaling = adata.var['fit_scaling'][idx]
        self.t_ = adata.var['fit_t_'][idx]
        self.steady_state_ratio = self.gamma / self.beta

        self.u0, self.s0, self.alpha_ = 0, 0, 0
        self.u0_, self.s0_ = mRNA(self.t_, self.u0, self.s0, self.alpha, self.beta, self.gamma)
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
            print(self.t_)
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

    def get_divergence(self, alpha=None, beta=None, gamma=None, scaling=None, t_=None, u0_=None, s0_=None, mode=None):
        alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
        res = compute_divergence(self.u / scaling, self.s, alpha, beta, gamma, scaling, t_, u0_, s0_, mode=mode,
                                 std_u=self.std_u, std_s=self.std_s, assignment_mode=self.assignment_mode,
                                 connectivities=self.connectivities, fit_steady_states=self.fit_steady_states)
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
            u, s, t = u[idx], s[idx], t[idx]
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
        return np.sum(udiff ** 2) + np.sum(sdiff ** 2) + np.sum(reg ** 2)

    def get_mse(self, **kwargs):
        return self.get_se(**kwargs) / np.sum(self.weights)

    def get_loss(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None, refit_time=None):
        return self.get_se(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma, scaling=scaling, u0_=u0_, s0_=s0_, refit_time=refit_time)

    def get_loglikelihood(self, **kwargs):
        if 'weighted' not in kwargs: kwargs.update({'weighted': 'dynamical'})
        udiff, sdiff, reg = self.get_dists(**kwargs)
        distx = udiff ** 2 + sdiff ** 2
        varx = np.mean(distx) - np.mean(np.sign(sdiff) * np.sqrt(distx))**2  # np.var(np.sign(sdiff) * np.sqrt(distx))
        return - 1 / 2 / len(distx) * np.sum(distx) / varx - 1 / 2 * np.log(2 * np.pi * varx)

    def get_likelihood(self, **kwargs):
        if 'weighted' not in kwargs: kwargs.update({'weighted': 'dynamical'})
        likelihood = np.exp(self.get_loglikelihood(**kwargs))
        return likelihood

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
                   show=True, **kwargs):
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
        return ax if show is False else None

    def plot_contours(self, xkey='gamma', ykey='alpha', x_sight=.5, y_sight=.5, num=20, dpi=None,
                      fontsize=8, show_unprofiled=False, refit_time=None, **kwargs):
        from ..plotting.utils import update_axes
        x_var = eval('self.' + xkey)
        y_var = eval('self.' + ykey)

        x = np.linspace(-x_sight, x_sight, num=num) * x_var + x_var
        y = np.linspace(-y_sight, y_sight, num=num) * y_var + y_var

        assignment_mode = self.assignment_mode
        self.assignment_mode = None

        fp = lambda x, y: self.get_loss(**{xkey: x, ykey: y}, refit_time=refit_time)

        zp = np.zeros((len(x), len(x)))
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
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
        ax1 = update_axes(ax1, fontsize=fontsize, frameon=True)

        if show_unprofiled:
            x, y, refit_time = x, y, False
        else:
            x = np.linspace(-x_sight / 5, x_sight / 5, num=num) * x_var + x_var
            y = np.linspace(-y_sight / 5, y_sight / 5, num=num) * y_var + y_var
        f0 = lambda x, y: self.get_loss(**{xkey: x, ykey: y}, refit_time=refit_time)

        z0 = np.zeros((len(x), len(x)))
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                z0[i, j] = f0(xi, yi)

        ax2.contourf(x, y, np.log1p(z0.T), levels=20, cmap='RdGy_r')
        contours = ax2.contour(x, y, np.log1p(z0.T), 4, colors='k', linewidths=.5)
        ax2.clabel(contours, inline=True, fontsize=fontsize * .75)
        ax2.scatter(x=x_var, y=y_var, s=50, c='purple', zorder=3, **kwargs)
        ax2.set_xlabel(x_label, fontsize=fontsize)
        ax2.set_ylabel(y_label, fontsize=fontsize)
        ax2.set_title('MSE', fontsize=fontsize)
        ax2 = update_axes(ax2, fontsize=fontsize, frameon=True)

        ix, iy = np.unravel_index(zp.argmin(), zp.shape)
        x_opt, y_opt = x[ix].round(2), y[ix].round(2)

        self.assignment_mode = assignment_mode
        return x_opt, y_opt
