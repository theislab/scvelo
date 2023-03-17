import warnings

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats.distributions import chi2, norm

import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib import gridspec, rcParams

from scvelo import logging as logg
from scvelo.core import clipped_log, invert, SplicingDynamics
from scvelo.preprocessing.moments import get_connectivities
from .utils import make_dense, round

exp = np.exp


"""Helper functions"""


# TODO: Add docstrings
def normalize(X, axis=0, min_confidence=None):
    """TODO."""
    X_sum = np.sum(X, axis=axis)
    if min_confidence:
        X_sum += min_confidence
    X_sum += X_sum == 0
    return X / X_sum


# TODO: Add docstrings
def convolve(x, weights=None):
    """TODO."""
    if weights is None:
        return x
    else:
        return weights.multiply(x).tocsr() if issparse(weights) else weights * x


# TODO: Add docstrings
def linreg(u, s):  # linear regression fit
    """TODO."""
    ss_ = s.multiply(s).sum(0) if issparse(s) else (s**2).sum(0)
    us_ = s.multiply(u).sum(0) if issparse(s) else (s * u).sum(0)
    return us_ / ss_


# TODO: Add docstrings
def compute_dt(t, clipped=True, axis=0):
    """TODO."""
    prepend = np.min(t, axis=axis)[None, :]
    dt = np.diff(np.sort(t, axis=axis), prepend=prepend, axis=axis)
    m_dt = np.max([np.mean(dt, axis=axis), np.max(t, axis=axis) / len(t)], axis=axis)
    m_dt = np.clip(m_dt, 0, None)
    if clipped:  # Poisson upper bound
        ub = m_dt + 3 * np.sqrt(m_dt)
        dt = np.clip(dt, 0, ub)
    return dt


# TODO: Add docstrings
def root_time(t, root=None):
    """TODO."""
    nans = np.isnan(np.sum(t, axis=0))
    if np.any(nans):
        t = t[:, ~nans]

    t_root = 0 if root is None else t[root]
    o = np.array(t >= t_root, dtype=int)
    t_after = (t - t_root) * o
    t_origin = np.max(t_after, axis=0)
    t_before = (t + t_origin) * (1 - o)

    t_switch = np.min(t_before, axis=0)
    t_rooted = t_after + t_before
    return t_rooted, t_switch


# TODO: Add docstrings
def compute_shared_time(t, perc=None, norm=True):
    """TODO."""
    nans = np.isnan(np.sum(t, axis=0))
    if np.any(nans):
        t = np.array(t[:, ~nans])
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
    if norm:
        t_shared /= t_shared.max()

    return t_shared


"""Dynamics delineation"""


# TODO: Add docstrings
def unspliced(tau, u0, alpha, beta):
    """TODO."""
    expu = exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


# TODO: Add docstrings
def spliced(tau, s0, u0, alpha, beta, gamma):
    """TODO."""
    c = (alpha - u0 * beta) * invert(gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


# TODO: Add docstrings
def adjust_increments(tau, tau_=None):
    """TODO."""
    tau_new = np.array(tau)
    tau_ord = np.sort(tau_new)
    dtau = np.diff(tau_ord, prepend=0)

    if tau_ is None:
        # m_dtau = np.max([np.mean(dtau), np.max(tau) / len(tau), 0])
        # ub = m_dtau + 3 * np.sqrt(m_dtau)  # Poisson with std = sqrt(mean)
        ub = 3 * np.percentile(dtau, 99.5, axis=0)
    else:
        tau_new_ = np.array(tau_)
        tau_ord_ = np.sort(tau_new_)
        dtau_ = np.diff(tau_ord_, prepend=0)

        # m_dtaus = [m_dtau, np.max([np.mean(dtau_), np.max(tau_) / len(tau_), 0])]
        # m_dtau = np.min(m_dtaus)
        # ub = m_dtau + 3 * np.sqrt(m_dtau)  # Poisson with std = sqrt(mean)
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


# TODO: Add docstrings
def tau_inv(u, s=None, u0=None, s0=None, alpha=None, beta=None, gamma=None):
    """TODO."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inv_u = (gamma >= beta) if gamma is not None else True
        inv_us = np.invert(inv_u)
    any_invu = np.any(inv_u) or s is None
    any_invus = np.any(inv_us) and s is not None

    if any_invus:  # tau_inv(u, s)
        beta_ = beta * invert(gamma - beta)
        xinf = alpha / gamma - beta_ * (alpha / beta)
        tau = (
            -1 / gamma * clipped_log((s - beta_ * u - xinf) / (s0 - beta_ * u0 - xinf))
        )

    if any_invu:  # tau_inv(u)
        uinf = alpha / beta
        tau_u = -1 / beta * clipped_log((u - uinf) / (u0 - uinf))
        tau = tau_u * inv_u + tau * inv_us if any_invus else tau_u
    return tau


# TODO: Add docstrings
def assign_tau(
    u, s, alpha, beta, gamma, t_=None, u0_=None, s0_=None, assignment_mode=None
):
    """TODO."""
    if assignment_mode in {"full_projection", "partial_projection"} or (
        assignment_mode == "projection" and beta < gamma
    ):
        x_obs = np.vstack([u, s]).T
        t0 = tau_inv(np.min(u[s > 0]), u0=u0_, alpha=0, beta=beta)

        num = np.clip(int(len(u) / 5), 200, 500)
        tpoints = np.linspace(0, t_, num=num)
        tpoints_ = np.linspace(0, t0, num=num)[1:]
        xt = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(tpoints)
        xt_ = SplicingDynamics(
            alpha=0, beta=beta, gamma=gamma, initial_state=[u0_, s0_]
        ).get_solution(tpoints_)

        # assign time points (oth. projection onto 'on' and 'off' curve)
        tau = tpoints[
            ((xt[None, :, :] - x_obs[:, None, :]) ** 2).sum(axis=2).argmin(axis=1)
        ]
        tau_ = tpoints_[
            ((xt_[None, :, :] - x_obs[:, None, :]) ** 2).sum(axis=2).argmin(axis=1)
        ]
    else:
        tau = tau_inv(u, s, 0, 0, alpha, beta, gamma)
        tau = np.clip(tau, 0, t_)

        tau_ = tau_inv(u, s, u0_, s0_, 0, beta, gamma)
        tau_ = np.clip(tau_, 0, np.max(tau_[s > 0]))

    return tau, tau_, t_


def compute_divergence(
    u,
    s,
    alpha,
    beta,
    gamma,
    scaling=1,
    t_=None,
    u0_=None,
    s0_=None,
    tau=None,
    tau_=None,
    std_u=1,
    std_s=1,
    normalized=False,
    mode="distance",
    assignment_mode=None,
    var_scale=False,
    kernel_width=None,
    fit_steady_states=True,
    connectivities=None,
    constraint_time_increments=True,
    reg_time=None,
    reg_par=None,
    min_confidence=None,
    pval_steady=None,
    steady_u=None,
    steady_s=None,
    noise_model="chi",
    time_connectivities=None,
    clusters=None,
    **kwargs,
):
    """Estimates the divergence of ODE to observations (avaiable metrics: distance, mse, likelihood, loglikelihood).

    Arguments:
    ---------
    mode: `'distance'`, `'mse'`, `'likelihood'` (default: `'distance'`)

    """
    # set tau, tau_
    if u0_ is None or s0_ is None:
        u0_, s0_ = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
            t_, stacked=False
        )
    if tau is None or tau_ is None or t_ is None:
        tau, tau_, t_ = assign_tau(
            u, s, alpha, beta, gamma, t_, u0_, s0_, assignment_mode
        )

    std_u /= scaling

    # adjust increments of tau, tau_ to avoid meaningless jumps
    if constraint_time_increments:
        ut, st = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
            tau, stacked=False
        )
        ut_, st_ = SplicingDynamics(
            alpha=0, beta=beta, gamma=gamma, initial_state=[u0_, s0_]
        ).get_solution(tau_, stacked=False)

        distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
        dists, dists_ = (s - st) / std_s, (s - st_) / std_s

        res = np.array([distu_**2 + dists_**2, distu**2 + dists**2])
        if connectivities is not None and connectivities is not False:
            res = (
                np.array([connectivities.dot(r) for r in res])
                if res.ndim > 2
                else connectivities.dot(res.T).T
            )

        o = np.argmin(res, axis=0)

        off, on = o == 0, o == 1
        if np.any(on) and np.any(off):
            tau[on], tau_[off] = adjust_increments(tau[on], tau_[off])
        elif np.any(on):
            tau[on] = adjust_increments(tau[on])
        elif np.any(off):
            tau_[off] = adjust_increments(tau_[off])

    # compute induction/repression state distances
    ut, st = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
        tau, stacked=False
    )
    ut_, st_ = SplicingDynamics(
        alpha=0, beta=beta, gamma=gamma, initial_state=[u0_, s0_]
    ).get_solution(tau_, stacked=False)

    if ut.ndim > 1 and ut.shape[1] == 1:
        ut = np.ravel(ut)
        st = np.ravel(st)
    if ut_.ndim > 1 and ut_.shape[1] == 1:
        ut_ = np.ravel(ut_)
        st_ = np.ravel(st_)

    distu, distu_ = (u - ut) / std_u, (u - ut_) / std_u
    dists, dists_ = (s - st) / std_s, (s - st_) / std_s

    if mode == "unspliced_dists":
        return distu, distu_

    elif mode == "outside_of_trajectory":
        return np.sign(distu) * np.sign(distu_) == 1

    distx = distu**2 + dists**2
    distx_ = distu_**2 + dists_**2

    res, varx = np.array([distx_, distx]), 1  # default vals;

    if noise_model == "normal":
        if var_scale:
            o = np.argmin([distx_, distx], axis=0)
            varu = np.nanvar(distu * o + distu_ + (1 - o), axis=0)
            vars = np.nanvar(dists * o + dists_ + (1 - o), axis=0)

            distx = distu**2 / varu + dists**2 / vars
            distx_ = distu_**2 / varu + dists_**2 / vars

            varx = varu * vars

            std_u *= np.sqrt(varu)
            std_s *= np.sqrt(vars)

    # compute steady state distances
    if fit_steady_states:
        u_inf, s_inf = alpha / beta, alpha / gamma
        distx_steady = ((u - u_inf) / std_u) ** 2 + ((s - s_inf) / std_s) ** 2
        distx_steady_ = (u / std_u) ** 2 + (s / std_s) ** 2
        res = np.array([distx_, distx, distx_steady_, distx_steady])

    if connectivities is not None and connectivities is not False:
        res = (
            np.array([connectivities.dot(r) for r in res])
            if res.ndim > 2
            else connectivities.dot(res.T).T
        )

    # compute variances
    if noise_model == "chi":
        if var_scale:
            o = np.argmin([distx_, distx], axis=0)
            dist = distx * o + distx_ * (1 - o)
            sign = np.sign(dists * o + dists_ * (1 - o))
            varx = np.mean(dist, axis=0) - np.mean(sign * np.sqrt(dist), axis=0) ** 2
            if kernel_width is not None:
                varx *= kernel_width**2
            res /= varx
        elif kernel_width is not None:
            res /= kernel_width**2

    if reg_time is not None and len(reg_time) == len(distu_):
        o = np.argmin(res, axis=0)
        t_max = (t_ + tau_) * (o == 0)
        t_max /= np.max(t_max, axis=0)
        reg_time /= np.max(reg_time)

        dist_tau = (tau - reg_time[:, None]) ** 2
        dist_tau_ = (tau_ + t_ - reg_time[:, None]) ** 2
        mu_res = np.mean(res, axis=1)
        if reg_par is not None:
            mu_res *= reg_par

        res[0] += dist_tau_ * mu_res[0]
        res[1] += dist_tau * mu_res[1]
        if fit_steady_states:
            res[2] += dist_tau * mu_res[1]
            res[3] += dist_tau_ * mu_res[0]

    if mode == "tau":
        res = [tau, tau_]

    elif mode == "likelihood":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)

    elif mode == "nll":
        res = np.log(2 * np.pi * np.sqrt(varx)) + 0.5 * res
        if normalized:
            res = normalize(res, min_confidence=min_confidence)

    elif mode == "confidence":
        res = np.array([res[0], res[1]])
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        res = np.median(
            np.max(res, axis=0) - (np.sum(res, axis=0) - np.max(res, axis=0)), axis=1
        )

    elif mode in {"soft_eval", "soft"}:
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)

        o_, o = res[0], res[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode in {"hardsoft_eval", "hardsoft"}:
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        o = np.argmax(res, axis=0)
        o_, o = (o == 0) * res[0], (o == 1) * res[1]
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode in {"hard_eval", "hard"}:
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        o = np.argmax(res, axis=0)
        o_, o = o == 0, o == 1
        res = np.array([o_, o, ut * o + ut_ * o_, st * o + st_ * o_])

    elif mode == "soft_state":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        res = res[1] - res[0]

    elif mode == "hard_state":
        res = np.argmin(res, axis=0)

    elif mode == "steady_state":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        res = res[2] + res[3]

    elif mode in {"assign_timepoints", "time"}:
        o = np.argmin(res, axis=0)

        tau_ *= o == 0
        tau *= o == 1

        if 2 in o:
            o[o == 2] = 1
        if 3 in o:
            o[o == 3] = 0

        t = tau * (o == 1) + (tau_ + t_) * (o == 0)
        res = [t, tau, o] if mode == "assign_timepoints" else t

    elif mode == "dists":
        o = np.argmin(res, axis=0)

        tau_ *= o == 0
        tau *= o == 1

        if 2 in o:
            o[o == 2] = 1
        if 3 in o:
            o[o == 3] = 0

        distu = distu * (o == 1) + distu_ * (o == 0)
        dists = dists * (o == 1) + dists_ * (o == 0)
        res = distu, dists

    elif mode == "distx":
        o = np.argmin(res, axis=0)

        tau_ *= o == 0
        tau *= o == 1

        if 2 in o:
            o[o == 2] = 1
        if 3 in o:
            o[o == 3] = 0

        distu = distu * (o == 1) + distu_ * (o == 0)
        dists = dists * (o == 1) + dists_ * (o == 0)
        res = distu**2 + dists**2

    elif mode == "gene_likelihood":
        o = np.argmin(res, axis=0)

        tau_ *= o == 0
        tau *= o == 1

        if 2 in o:
            o[o == 2] = 1
        if 3 in o:
            o[o == 3] = 0

        distu = distu * (o == 1) + distu_ * (o == 0)
        dists = dists * (o == 1) + dists_ * (o == 0)

        idx = np.array((u > np.max(u, 0) / 5) & (s > np.max(s, 0) / 5), dtype=int)
        idx = idx / idx
        distu *= idx
        dists *= idx

        distx = distu**2 + dists**2

        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        varx = (
            np.nanmean(distx, 0) - np.nanmean(np.sign(dists) * np.sqrt(distx), 0) ** 2
        )

        if clusters is not None:
            res = []
            for cat in clusters.cat.categories:
                idx_cat = np.array(clusters == cat)
                distx_cat = distu[idx_cat] ** 2 + dists[idx_cat] ** 2
                distx_sum = np.nansum(distx_cat, 0)

                # penalize if very low count number
                n = np.sum(np.invert(np.isnan(distx_cat)), 0) - len(distx_cat) * 0.01
                n = np.clip(n, 2, None)
                distx_sum[n < np.nanmax(n) / 5] = np.nanmax(distx_sum)
                ll = -1 / 2 / n * distx_sum / varx - 1 / 2 * np.log(2 * np.pi * varx)
                ll[distx_sum == 0] = np.nan
                res.append(ll)
            res = np.exp(res)
        else:
            n = np.clip(len(distu) - len(distu) * 0.01, 2, None)
            ll = -1 / 2 / n * np.nansum(distx, 0) / varx
            ll -= 1 / 2 * np.log(2 * np.pi * varx)
            res = np.exp(ll)

    elif mode == "velocity":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        res = (
            np.array([res[0], res[1], res[2], res[3], 1e-6 * np.ones(res[0].shape)])
            if res.ndim > 2
            else np.array([res[0], res[1], min_confidence * np.ones(res[0].shape)])
        )
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        res = np.argmax(res, axis=0)
        o_, o = res == 0, res == 1
        t = tau * o + (tau_ + t_) * o_
        if time_connectivities:
            if time_connectivities is True:
                time_connectivities = connectivities
            t = time_connectivities.dot(t)
            o = (res < 2) * (t < t_)
            o_ = (res < 2) * (t >= t_)
            tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
            ut, st = SplicingDynamics(
                alpha=alpha, beta=beta, gamma=gamma, initial_state=[u0, s0]
            ).get_solution(tau, stacked=False)
            ut_, st_ = ut, st

        ut = ut * o + ut_ * o_
        st = st * o + st_ * o_
        alpha = alpha * o

        vt = ut * beta - st * gamma  # ds/dt
        wt = (alpha - beta * ut) * scaling  # du/dt

        vt, wt = np.clip(vt, -s, None), np.clip(wt, -u * scaling, None)
        if vt.ndim == 1:
            vt = vt.reshape(len(vt), 1)
            wt = wt.reshape(len(wt), 1)

        res = [vt, wt]

    elif mode == "velocity_residuals":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        res = (
            np.array([res[0], res[1], res[2], res[3], 1e-6 * np.ones(res[0].shape)])
            if res.ndim > 2
            else np.array([res[0], res[1], min_confidence * np.ones(res[0].shape)])
        )
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        res = np.argmax(res, axis=0)
        o_, o = res == 0, res == 1
        t = tau * o + (tau_ + t_) * o_
        if time_connectivities:
            if time_connectivities is True:
                time_connectivities = connectivities
            t = time_connectivities.dot(t)
            o = (res < 2) * (t < t_)
            o_ = (res < 2) * (t >= t_)
            tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
            ut, st = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
                tau, stacked=False
            )
            ut_, st_ = ut, st

        alpha = alpha * o

        vt = u * beta - s * gamma  # ds/dt
        wt = (alpha - beta * u) * scaling  # du/dt

        vt, wt = np.clip(vt, -s, None), np.clip(wt, -u * scaling, None)

        res = [vt, wt]

    elif mode == "soft_velocity":
        res = 1 / (2 * np.pi * np.sqrt(varx)) * np.exp(-0.5 * res)
        if normalized:
            res = normalize(res, min_confidence=min_confidence)
        o_, o = res[0], res[1]
        ut = ut * o + ut_ * o_
        st = st * o + st_ * o_
        alpha = alpha * o
        vt = ut * beta - st * gamma  # ds/dt
        wt = (alpha - beta * ut) * scaling  # du/dt
        res = [vt, wt]

    return res


# TODO: Add docstrings
def assign_timepoints(**kwargs):
    """TODO."""
    return compute_divergence(**kwargs, mode="assign_timepoints")


# TODO: Add docstrings
def vectorize(t, t_, alpha, beta, gamma=None, alpha_=0, u0=0, s0=0, sorted=False):
    """TODO."""
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


# TODO: Add docstrings
def curve_dists(
    u,
    s,
    alpha,
    beta,
    gamma,
    t_=None,
    u0_=None,
    s0_=None,
    std_u=1,
    std_s=1,
    scaling=1,
    num=None,
):
    """TODO."""
    if u0_ is None or s0_ is None:
        u0_, s0_ = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
            t_, stacked=False
        )

    x_obs = np.vstack([u, s]).T
    std_x = np.vstack([std_u / scaling, std_s]).T
    t0 = tau_inv(np.min(u[s > 0]), u0=u0_, alpha=0, beta=beta)

    num = np.clip(int(len(u) / 10), 50, 200) if num is None else num
    tpoints = np.linspace(0, t_, num=num)
    tpoints_ = np.linspace(0, t0, num=num)[1:]

    curve_t = SplicingDynamics(alpha=alpha, beta=beta, gamma=gamma).get_solution(
        tpoints
    )
    curve_t_ = SplicingDynamics(
        alpha=0, beta=beta, gamma=gamma, initial_state=[u0_, s0_]
    ).get_solution(tpoints_)

    # match each curve point to nearest observation
    dist, dist_ = np.zeros(len(curve_t)), np.zeros(len(curve_t_))
    for i, ci in enumerate(curve_t):
        dist[i] = np.min(np.sum((x_obs - ci) ** 2 / std_x**2, 1))
    for i, ci in enumerate(curve_t_):
        dist_[i] = np.min(np.sum((x_obs - ci) ** 2 / std_x**2, 1))

    return dist, dist_


"""Base Class for Dynamics Recovery"""


# TODO: Add docstrings
class BaseDynamics:
    """TODO."""

    def __init__(
        self,
        adata,
        gene,
        u=None,
        s=None,
        use_raw=False,
        perc=99,
        max_iter=10,
        fit_time=True,
        fit_scaling=True,
        fit_steady_states=True,
        fit_connected_states=True,
        fit_basal_transcription=None,
        high_pars_resolution=False,
        steady_state_prior=None,
        init_vals=None,
    ):
        self.s, self.u, self.use_raw = None, None, None

        _layers = adata[:, gene].layers
        self.gene = gene
        self.use_raw = use_raw or "Ms" not in _layers.keys()

        # extract actual data
        if u is None or s is None:
            u = _layers["unspliced"] if self.use_raw else _layers["Mu"]
            s = _layers["spliced"] if self.use_raw else _layers["Ms"]
        self.s, self.u = make_dense(s), make_dense(u)

        # Basal transcription
        if fit_basal_transcription:
            self.u0, self.s0 = np.min(u), np.min(s)
            self.u -= self.u0
            self.s -= self.s0
        else:
            self.u0, self.s0 = 0, 0

        self.alpha, self.beta, self.gamma = None, None, None
        self.scaling, self.t_, self.alpha_ = None, None, None
        self.u0_, self.s0_, self.weights = None, None, None
        self.weights_outer, self.weights_upper = None, None
        self.t, self.tau, self.o, self.tau_ = None, None, None, None
        self.likelihood, self.loss, self.pars = None, None, None

        self.max_iter = max_iter
        # partition to total of 5 fitting procedures
        # (t_ and alpha, scaling, rates, t_, all together)
        self.simplex_kwargs = {
            "method": "Nelder-Mead",
            "options": {"maxiter": int(self.max_iter / 5)},
        }

        self.perc = perc
        self.recoverable = True
        # TODO: Refactor to handle failed instantiation properly
        try:
            self.initialize_weights()
        except Warning:
            self.recoverable = False
            logg.warn(f"Model for {self.gene} could not be instantiated.")

        self.refit_time = fit_time

        self.assignment_mode = None
        self.steady_state_ratio = None
        self.steady_state_prior = steady_state_prior

        self.fit_scaling = fit_scaling
        self.fit_steady_states = fit_steady_states
        self.fit_connected_states = fit_connected_states
        self.connectivities = (
            get_connectivities(adata)
            if self.fit_connected_states is True
            else self.fit_connected_states
        )
        self.high_pars_resolution = high_pars_resolution
        self.init_vals = init_vals

        # for differential kinetic test
        self.clusters, self.cats, self.varx, self.orth_beta = None, None, None, None
        self.diff_kinetics, self.pval_kinetics, self.pvals_kinetics = None, None, None

    # TODO: Add docstrings
    def initialize_weights(self, weighted=True):
        """TODO."""
        nonzero_s = np.ravel(self.s > 0)
        nonzero_u = np.ravel(self.u > 0)

        weights = np.array(nonzero_s & nonzero_u, dtype=bool)
        self.recoverable = np.sum(weights) > 2

        if self.recoverable:
            if weighted:
                ub_s = np.percentile(self.s[weights], self.perc)
                ub_u = np.percentile(self.u[weights], self.perc)
                if ub_s > 0:
                    weights &= np.ravel(self.s <= ub_s)
                if ub_u > 0:
                    weights &= np.ravel(self.u <= ub_u)

            self.weights = weights
            u, s = self.u[weights], self.s[weights]
            self.std_u = np.std(u)
            self.std_s = np.std(s)

            self.weights_upper = np.array(weights)
            if np.any(weights):
                w_upper = (self.u > np.max(u) / 3) & (self.s > np.max(s) / 3)
                self.weights_upper &= w_upper

    # TODO: Add docstrings
    def load_pars(self, adata, gene):
        """TODO."""
        idx = adata.var_names.get_loc(gene) if isinstance(gene, str) else gene
        self.alpha = adata.var["fit_alpha"][idx]
        self.beta = adata.var["fit_beta"][idx] * adata.var["fit_scaling"][idx]
        self.gamma = adata.var["fit_gamma"][idx]
        self.scaling = adata.var["fit_scaling"][idx]
        self.t_ = adata.var["fit_t_"][idx]
        self.steady_state_ratio = self.gamma / self.beta

        if "fit_steady_u" in adata.var.keys():
            self.steady_u = adata.var["fit_steady_u"][idx]
        if "fit_steady_s" in adata.var.keys():
            self.steady_s = adata.var["fit_steady_s"][idx]
        if "fit_pval_steady" in adata.var.keys():
            self.pval_steady = adata.var["fit_pval_steady"][idx]

        self.alpha_ = 0
        self.u0_, self.s0_ = SplicingDynamics(
            alpha=self.alpha,
            beta=self.beta,
            gamma=self.gamma,
        ).get_solution(self.t_, stacked=False)
        self.pars = [self.alpha, self.beta, self.gamma, self.t_, self.scaling]
        self.pars = np.array(self.pars)[:, None]

        lt = "latent_time"
        t = adata.obs[lt] if lt in adata.obs.keys() else adata.layers["fit_t"][:, idx]

        if isinstance(self.refit_time, bool):
            self.t, self.tau, self.o = self.get_time_assignment(t=t)
        else:
            tkey = self.refit_time
            self.t = adata.obs[tkey].values if isinstance(tkey, str) else tkey
            self.refit_time = False
            steady_states = t == self.t_
            if np.any(steady_states):
                self.t_ = np.mean(self.t[steady_states])
            self.t, self.tau, self.o = self.get_time_assignment(t=self.t)

        self.loss = [self.get_loss()]

    # TODO: Add docstrings
    def get_weights(self, weighted=None, weights_cluster=None):
        """TODO."""
        weights = (
            np.array(
                self.weights_outer
                if weighted == "outer"
                else self.weights_upper
                if weighted == "upper"
                else self.weights
            )
            if weighted
            else np.ones(len(self.weights), bool)
        )
        if weights_cluster is not None and len(weights) == len(weights_cluster):
            weights &= weights_cluster
        return weights

    # TODO: Add docstrings
    def get_reads(self, scaling=None, weighted=None, weights_cluster=None):
        """TODO."""
        scaling = self.scaling if scaling is None else scaling
        u, s = self.u / scaling, self.s
        if weighted or weights_cluster is not None:
            weights = self.get_weights(
                weighted=weighted, weights_cluster=weights_cluster
            )
            u, s = u[weights], s[weights]
        return u, s

    # TODO: Add docstrings
    def get_vars(
        self,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        t_=None,
        u0_=None,
        s0_=None,
    ):
        """TODO."""
        alpha = self.alpha if alpha is None else alpha
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma
        scaling = self.scaling if scaling is None else scaling
        if t_ is None or t_ == 0:
            t_ = self.t_ if u0_ is None else tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)
        return alpha, beta, gamma, scaling, t_

    # TODO: Add docstrings
    def get_divergence(
        self,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        t_=None,
        u0_=None,
        s0_=None,
        mode=None,
        **kwargs,
    ):
        """TODO."""
        alpha, beta, gamma, scaling, t_ = self.get_vars(
            alpha, beta, gamma, scaling, t_, u0_, s0_
        )
        u, s = self.u / scaling, self.s
        kwargs.update(
            {"t_": t_, "u0_": u0_, "s0_": s0_, "std_u": self.std_u, "std_s": self.std_s}
        )
        kwargs.update(
            {
                "mode": mode,
                "assignment_mode": self.assignment_mode,
                "connectivities": self.connectivities,
                "fit_steady_states": self.fit_steady_states,
            }
        )
        res = compute_divergence(u, s, alpha, beta, gamma, scaling, **kwargs)
        return res

    # TODO: Add docstrings
    def get_time_assignment(
        self,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        t_=None,
        u0_=None,
        s0_=None,
        t=None,
        refit_time=None,
        rescale_factor=None,
        weighted=None,
        weights_cluster=None,
    ):
        """TODO."""
        if refit_time is None:
            refit_time = self.refit_time

        if t is not None:
            t_ = self.t_ if t_ is None else t_
            o = np.array(t < t_, dtype=int)
            tau = t * o + (t - t_) * (1 - o)
        elif refit_time:
            if rescale_factor is not None:
                alpha, beta, gamma, scaling, t_ = self.get_vars(
                    alpha, beta, gamma, scaling, t_, u0_, s0_
                )

                u0_ = self.u0_ if u0_ is None else u0_
                s0_ = self.s0_ if s0_ is None else s0_
                rescale_factor *= gamma / beta

                scaling *= rescale_factor
                beta *= rescale_factor

                u0_ /= rescale_factor
                t_ = tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)

            t, tau, o = self.get_divergence(
                alpha, beta, gamma, scaling, t_, u0_, s0_, mode="assign_timepoints"
            )
            if rescale_factor is not None:
                t *= self.t_ / t_
                tau *= self.t_ / t_
        else:
            t, tau, o = self.t, self.tau, self.o

        if weighted or weights_cluster is not None:
            weights = self.get_weights(
                weighted=weighted, weights_cluster=weights_cluster
            )
            t, tau, o = t[weights], tau[weights], o[weights]
        return t, tau, o

    # TODO: Add docstrings
    def get_vals(
        self,
        t=None,
        t_=None,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        u0_=None,
        s0_=None,
        refit_time=None,
    ):
        """TODO."""
        alpha, beta, gamma, scaling, t_ = self.get_vars(
            alpha, beta, gamma, scaling, t_, u0_, s0_
        )
        t, tau, o = self.get_time_assignment(
            alpha, beta, gamma, scaling, t_, u0_, s0_, t, refit_time
        )
        return t, t_, alpha, beta, gamma, scaling

    # TODO: Add docstrings
    def get_dists(
        self,
        t=None,
        t_=None,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        u0_=None,
        s0_=None,
        refit_time=None,
        weighted=True,
        weights_cluster=None,
        reg=None,
    ):
        """TODO."""
        weight_args = {"weighted": weighted, "weights_cluster": weights_cluster}
        u, s = self.get_reads(scaling, **weight_args)

        alpha, beta, gamma, scaling, t_ = self.get_vars(
            alpha, beta, gamma, scaling, t_, u0_, s0_
        )
        t, tau, o = self.get_time_assignment(
            alpha, beta, gamma, scaling, t_, u0_, s0_, t, refit_time, **weight_args
        )

        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        ut, st = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=[u0, s0]
        ).get_solution(tau, stacked=False)

        udiff = np.array(ut - u) / self.std_u * scaling
        sdiff = np.array(st - s) / self.std_s
        if reg is None:
            reg = 0
            if self.steady_state_ratio is not None:
                reg = (gamma / beta - self.steady_state_ratio) * s / self.std_s
        return udiff, sdiff, reg

    # TODO: Add docstrings
    def get_residuals_linear(self, **kwargs):
        """TODO."""
        udiff, sdiff, reg = self.get_dists(**kwargs)
        return udiff, sdiff

    # TODO: Add docstrings
    def get_residuals(self, **kwargs):
        """TODO."""
        udiff, sdiff, reg = self.get_dists(**kwargs)
        return np.sign(sdiff) * np.sqrt(udiff**2 + sdiff**2)

    # TODO: Add docstrings
    def get_distx(self, noise_model="normal", regularize=True, **kwargs):
        """TODO."""
        udiff, sdiff, reg = self.get_dists(**kwargs)
        distx = udiff**2 + sdiff**2
        if regularize:
            distx += reg**2
        return np.sqrt(distx) if noise_model == "laplace" else distx

    # TODO: Add docstrings
    def get_se(self, **kwargs):
        """TODO."""
        return np.sum(self.get_distx(**kwargs))

    # TODO: Add docstrings
    def get_mse(self, **kwargs):
        """TODO."""
        return np.mean(self.get_distx(**kwargs))

    # TODO: Add docstrings
    def get_loss(
        self,
        t=None,
        t_=None,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        u0_=None,
        s0_=None,
        refit_time=None,
    ):
        """TODO."""
        kwargs = {
            "t": t,
            "t_": t_,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "scaling": scaling,
        }
        kwargs.update({"u0_": u0_, "s0_": s0_, "refit_time": refit_time})
        return self.get_se(**kwargs)

    # TODO: Add docstrings
    def get_loglikelihood(self, varx=None, noise_model="normal", **kwargs):
        """TODO."""
        if "weighted" not in kwargs:
            kwargs.update({"weighted": "upper"})
        udiff, sdiff, reg = self.get_dists(**kwargs)

        distx = udiff**2 + sdiff**2 + reg**2
        eucl_distx = np.sqrt(distx)
        n = np.clip(len(distx) - len(self.u) * 0.01, 2, None)

        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        if varx is None:
            varx = np.mean(distx) - np.mean(np.sign(sdiff) * eucl_distx) ** 2
        varx += varx == 0  # edge case of mRNAs levels to be the same across all cells

        if noise_model == "normal":
            loglik = -1 / 2 / n * np.sum(distx) / varx
            loglik -= 1 / 2 * np.log(2 * np.pi * varx)
        elif noise_model == "laplace":
            loglik = -1 / np.sqrt(2) / n * np.sum(eucl_distx) / np.sqrt(varx)
            loglik -= 1 / 2 * np.log(2 * varx)
        else:
            raise ValueError("That noise model is not supported.")
        return loglik

    # TODO: Add docstrings
    def get_likelihood(self, **kwargs):
        """TODO."""
        if "weighted" not in kwargs:
            kwargs.update({"weighted": "upper"})
        likelihood = np.exp(self.get_loglikelihood(**kwargs))
        return likelihood

    # TODO: Add docstrings
    def get_curve_likelihood(self):
        """TODO."""
        alpha, beta, gamma, scaling, t_ = self.get_vars()
        u, s = self.get_reads(scaling, weighted=False)
        varx = self.get_variance()

        kwargs = {"std_u": self.std_u, "std_s": self.std_s, "scaling": scaling}
        dist, dist_ = curve_dists(u, s, alpha, beta, gamma, t_, **kwargs)
        log_likelihood = -0.5 / len(dist) * np.sum(dist) / varx - 0.5 * np.log(
            2 * np.pi * varx
        )
        log_likelihood_ = -0.5 / len(dist_) * np.sum(dist_) / varx - 0.5 * np.log(
            2 * np.pi * varx
        )
        likelihood = np.exp(np.max([log_likelihood, log_likelihood_]))
        return likelihood

    # TODO: Add docstrings
    def get_variance(self, **kwargs):
        """TODO."""
        if "weighted" not in kwargs:
            kwargs.update({"weighted": "upper"})
        udiff, sdiff, reg = self.get_dists(**kwargs)
        distx = udiff**2 + sdiff**2
        return np.mean(distx) - np.mean(np.sign(sdiff) * np.sqrt(distx)) ** 2

    # TODO: Add docstrings
    def get_ut(self, **kwargs):
        """TODO."""
        t, t_, alpha, beta, gamma, scaling = self.get_vals(**kwargs)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return unspliced(tau, u0, alpha, beta)

    # TODO: Add docstrings
    def get_st(self, **kwargs):
        """TODO."""
        t, t_, alpha, beta, gamma, scaling = self.get_vals(**kwargs)
        tau, alpha, u0, s0 = vectorize(t, t_, alpha, beta, gamma)
        return spliced(tau, s0, u0, alpha, beta, gamma)

    # TODO: Add docstrings
    def get_vt(self, mode="soft_eval"):
        """TODO."""
        alpha, beta, gamma, scaling, _ = self.get_vars()
        o_, o, ut, st = self.get_divergence(mode=mode)
        return ut * beta - st * gamma

    # TODO: Add docstrings
    def plot_phase(
        self,
        adata=None,
        t=None,
        t_=None,
        alpha=None,
        beta=None,
        gamma=None,
        scaling=None,
        refit_time=None,
        show=True,
        show_assignments=None,
        **kwargs,
    ):
        """TODO."""
        from scvelo.plotting.scatter import scatter

        if np.all([x is None for x in [alpha, beta, gamma, scaling, t_]]):
            refit_time = False
        t, t_, alpha, beta, gamma, scaling = self.get_vals(
            t, t_, alpha, beta, gamma, scaling, refit_time=refit_time
        )
        ut = self.get_ut(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma) * scaling
        st = self.get_st(t=t, t_=t_, alpha=alpha, beta=beta, gamma=gamma)

        idx_sorted = np.argsort(t)
        ut, st, t = ut[idx_sorted], st[idx_sorted], t[idx_sorted]

        args = {"color_map": "RdYlGn", "vmin": -0.1, "vmax": 1.1}
        args.update(kwargs)

        ax = scatter(adata, x=self.s, y=self.u, colorbar=False, show=False, **kwargs)

        pl.plot(st, ut, color="purple", linestyle="-")
        pl.xlabel("spliced")
        pl.ylabel("unspliced")

        if show_assignments:
            xnew = (np.array([self.s[idx_sorted], st]),)
            ynew = np.array([self.u[idx_sorted], ut])
            ax.plot(xnew, ynew, color="grey", linewidth=0.1)

        return ax if not show else None

    # TODO: Add docstrings
    def plot_profile_contour(
        self,
        xkey="gamma",
        ykey="alpha",
        x_sight=0.5,
        y_sight=0.5,
        num=20,
        contour_levels=4,
        fontsize=12,
        refit_time=None,
        ax=None,
        color_map="RdGy",
        figsize=None,
        dpi=None,
        vmin=None,
        vmax=None,
        horizontal_ylabels=True,
        show_path=False,
        show=True,
        return_color_scale=False,
        **kwargs,
    ):
        """TODO."""
        from scvelo.plotting.utils import update_axes

        x_var = getattr(self, xkey)
        y_var = getattr(self, ykey)

        x = np.linspace(-x_sight, x_sight, num=num) * x_var + x_var
        y = np.linspace(-y_sight, y_sight, num=num) * y_var + y_var

        assignment_mode = self.assignment_mode
        self.assignment_mode = None

        # TODO: Check if list comprehension can be used
        zp = np.zeros((len(x), len(x)))
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                zp[i, j] = self.get_likelihood(
                    **{xkey: xi, ykey: yi}, refit_time=refit_time
                )
        log_zp = np.log1p(zp.T)

        if vmin is None:
            vmin = np.min(log_zp)
        if vmax is None:
            vmax = np.max(log_zp)

        x_label = r"$" + f"\\{xkey}$" if xkey in ["gamma", "alpha", "beta"] else xkey
        y_label = r"$" + f"\\{ykey}$" if ykey in ["gamma", "alpha", "beta"] else ykey

        if ax is None:
            figsize = rcParams["figure.figsize"] if figsize is None else figsize
            ax = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi).gca()

        ax.contourf(x, y, log_zp, levels=num, cmap=color_map, vmin=vmin, vmax=vmax)
        if contour_levels != 0:
            contours = ax.contour(
                x, y, log_zp, levels=contour_levels, colors="k", linewidths=0.5
            )
            fmt = "%1.1f" if np.isscalar(contour_levels) else "%1.0f"
            ax.clabel(contours, fmt=fmt, inline=True, fontsize=fontsize * 0.75)

        ax.scatter(x=x_var, y=y_var, s=50, c="purple", zorder=3, **kwargs)
        ax.set_xlabel(x_label, fontsize=fontsize)
        rotation = 0 if horizontal_ylabels else 90
        ax.set_ylabel(y_label, fontsize=fontsize, rotation=rotation)
        update_axes(ax, fontsize=fontsize, frameon=True)

        if show_path:
            axis = ax.axis()
            x_hist = self.pars[["alpha", "beta", "gamma", "t_", "scaling"].index(xkey)]
            y_hist = self.pars[["alpha", "beta", "gamma", "t_", "scaling"].index(ykey)]
            ax.plot(x_hist, y_hist)
            ax.axis(axis)

        self.assignment_mode = assignment_mode
        if return_color_scale:
            return np.min(log_zp), np.max(log_zp)
        elif not show:
            return ax

    # TODO: Add docstrings
    def plot_profile_hist(
        self,
        xkey="gamma",
        sight=0.5,
        num=20,
        dpi=None,
        fontsize=12,
        ax=None,
        figsize=None,
        color_map="RdGy",
        vmin=None,
        vmax=None,
        show=True,
    ):
        """TODO."""
        from scvelo.plotting.utils import update_axes

        x_var = getattr(self, xkey)
        x = np.linspace(-sight, sight, num=num) * x_var + x_var

        assignment_mode = self.assignment_mode
        self.assignment_mode = None

        # TODO: Check if list comprehension can be used
        zp = np.zeros(len(x))
        for i, xi in enumerate(x):
            zp[i] = self.get_likelihood(**{xkey: xi}, refit_time=True)

        log_zp = np.log1p(zp.T)
        if vmin is None:
            vmin = np.min(log_zp)
        if vmax is None:
            vmax = np.max(log_zp)

        x_label = r"$" + f"\\{xkey}$" if xkey in ["gamma", "alpha", "beta"] else xkey
        figsize = rcParams["figure.figsize"] if figsize is None else figsize
        if ax is None:
            fig = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi)
            ax = fig.gca()

        xp = np.linspace(x.min(), x.max(), 1000)
        yp = np.interp(xp, x, log_zp)
        ax.scatter(xp, yp, c=yp, cmap=color_map, edgecolor="none", vmin=vmin, vmax=vmax)
        ax.set_xlabel(x_label, fontsize=fontsize)
        ax.set_ylabel("likelihood", fontsize=fontsize)
        update_axes(ax, fontsize=fontsize, frameon=True)

        self.assignment_mode = assignment_mode
        if not show:
            return ax

    # TODO: Add docstrings
    def plot_profiles(
        self,
        params=None,
        contour_levels=0,
        sight=0.5,
        num=20,
        fontsize=12,
        color_map="RdGy",
        vmin=None,
        vmax=None,
        figsize=None,
        dpi=None,
        **kwargs,
    ):
        """TODO."""
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        if params is None:
            params = ["alpha", "beta", "gamma"]
        fig = pl.figure(constrained_layout=True, dpi=dpi, figsize=figsize)
        n = len(params)
        gs = gridspec.GridSpec(n, n, figure=fig)

        pkwargs = {
            "color_map": color_map,
            "vmin": vmin,
            "vmax": vmax,
            "figsize": figsize,
        }

        for i in range(len(params)):
            for j in range(n - 1, i - 1, -1):
                xkey = params[j]
                ykey = params[i]
                ax = fig.add_subplot(gs[n - 1 - i, n - 1 - j])
                if xkey == ykey:
                    ax = self.plot_profile_hist(
                        xkey,
                        ax=ax,
                        num=num,
                        sight=sight if np.isscalar(sight) else sight[j],
                        fontsize=fontsize,
                        show=False,
                        **pkwargs,
                    )
                    if i == 0 & j == 0:
                        cax = inset_axes(
                            ax,
                            width="7%",
                            height="100%",
                            loc="lower left",
                            bbox_to_anchor=(1.05, 0.0, 1, 1),
                            bbox_transform=ax.transAxes,
                            borderpad=0,
                        )
                        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                        cmap = mpl.cm.get_cmap(color_map)
                        _ = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
                else:
                    vmin_, vmax_ = self.plot_profile_contour(
                        xkey,
                        ykey,
                        ax=ax,
                        contour_levels=contour_levels,
                        x_sight=sight if np.isscalar(sight) else sight[j],
                        y_sight=sight if np.isscalar(sight) else sight[i],
                        num=num,
                        fontsize=fontsize,
                        return_color_scale=True,
                        **pkwargs,
                        **kwargs,
                    )
                    if vmin is None or vmax is None:
                        vmin, vmax = vmin_, vmax_  # scaled to first contour plot

                if i != 0:
                    ax.set_xlabel("")
                    ax.set_xticks([])
                if j - n + 1 != 0:
                    ax.set_ylabel("")
                    ax.set_yticks([])

    # TODO: Add docstrings
    def plot_state_likelihoods(
        self,
        num=300,
        dpi=None,
        figsize=None,
        color_map=None,
        color_map_steady=None,
        continuous=True,
        common_color_scale=True,
        var_scale=True,
        kernel_width=None,
        normalized=None,
        transitions=None,
        colorbar=False,
        alpha_=0.5,
        linewidths=3,
        padding_u=0.1,
        padding_s=0.1,
        fontsize=12,
        title=None,
        ax=None,
        **kwargs,
    ):
        """TODO."""
        from scvelo.plotting.utils import rgb_custom_colormap, update_axes

        if color_map is None:
            color_map = rgb_custom_colormap(
                ["royalblue", "white", "seagreen"], alpha=[1, 0.5, 1]
            )
        if color_map_steady is None:
            color_map_steady = rgb_custom_colormap(
                colors=3 * ["sienna"], alpha=[0, 0.5, 1]
            )

        alpha, beta, gamma, scaling, t_ = self.get_vars()
        u, s = self.u / scaling, self.s
        padding_u *= np.max(u) - np.min(u)
        padding_s *= np.max(s) - np.min(s)
        uu = np.linspace(np.min(u) - padding_u, np.max(u) + padding_u, num=num)
        ss = np.linspace(np.min(s) - padding_s, np.max(s) + padding_s, num=num)

        grid_u, grid_s = np.meshgrid(uu, ss)
        grid_u = grid_u.flatten()
        grid_s = grid_s.flatten()

        if var_scale:
            var_scale = self.get_variance()

        dkwargs = {
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "scaling": scaling,
            "t_": t_,
            "kernel_width": kernel_width,
            "std_u": self.std_u,
            "std_s": self.std_s,
            "var_scale": var_scale,
            "normalized": normalized,
            "fit_steady_states": True,
            "constraint_time_increments": False,
            "assignment_mode": "projection",
        }

        likelihoods = compute_divergence(u, s, mode="soft_state", **dkwargs)
        likelihoods_steady = compute_divergence(u, s, mode="steady_state", **dkwargs)

        likelihoods_grid = compute_divergence(
            grid_u, grid_s, mode="soft_state", **dkwargs
        )
        likelihoods_grid_steady = compute_divergence(
            grid_u, grid_s, mode="steady_state", **dkwargs
        )

        figsize = rcParams["figure.figsize"] if figsize is None else figsize
        if ax is None:
            fig = pl.figure(figsize=(figsize[0], figsize[1]), dpi=dpi)
            ax = fig.gca()

        ax.scatter(
            x=s,
            y=u,
            s=50,
            c=likelihoods_steady,
            zorder=3,
            cmap=color_map_steady,
            edgecolors="black",
            **kwargs,
        )
        ax.scatter(
            x=s,
            y=u,
            s=50,
            c=likelihoods,
            zorder=3,
            cmap=color_map,
            edgecolors="black",
            **kwargs,
        )

        l_grid, l_grid_steady = (
            likelihoods_grid.reshape(num, num).T,
            likelihoods_grid_steady.reshape(num, num).T,
        )

        if common_color_scale:
            vmax = vmax_steady = np.max(
                [np.abs(likelihoods_grid), np.abs(likelihoods_grid_steady)]
            )
        else:
            vmax, vmax_steady = np.max(np.abs(likelihoods_grid)), None

        if continuous:
            extent = (min(ss), max(ss), min(uu), max(uu))
            contf_steady = ax.imshow(
                l_grid_steady,
                cmap=color_map_steady,
                alpha=alpha_,
                vmin=0,
                vmax=vmax_steady,
                aspect="auto",
                origin="lower",
                extent=extent,
            )
            contf = ax.imshow(
                l_grid,
                cmap=color_map,
                alpha=alpha_,
                vmin=-vmax,
                vmax=vmax,
                aspect="auto",
                origin="lower",
                extent=extent,
            )
        else:
            cmap = color_map_steady
            contf_steady = ax.contourf(
                ss, uu, l_grid_steady, vmin=0, vmax=vmax_steady, levels=30, cmap=cmap
            )
            contf = ax.contourf(
                ss, uu, l_grid, vmin=-vmax, vmax=vmax, levels=30, cmap=color_map
            )

        # Contour lines
        if transitions is not None:
            transitions = np.multiply(
                np.array(transitions),
                [np.min(likelihoods_grid), np.max(likelihoods_grid)],
            )  # trans_width
            ax.contour(
                ss,
                uu,
                likelihoods_grid.reshape(num, num).T,
                transitions,
                linestyles="solid",
                colors="k",
                linewidths=linewidths,
            )

        if colorbar:
            pl.colorbar(contf, ax=ax)
            pl.colorbar(contf_steady, ax=ax)
        ax.set_xlabel("spliced", fontsize=fontsize)
        ax.set_ylabel("unspliced", fontsize=fontsize)
        title = "" if title is None else title
        ax.set_title(title, fontsize=fontsize)
        update_axes(ax, fontsize=fontsize, frameon=True)

        return ax

    # TODO: Add docstrings
    # for differential kinetic test
    def initialize_diff_kinetics(self, clusters):
        """TODO."""
        # after fitting dyn. model
        if self.varx is None:
            self.varx = self.get_variance()
        self.initialize_weights(weighted=False)
        self.steady_state_ratio = None
        self.clusters = clusters
        self.cats = pd.Categorical(clusters).categories
        self.weights_outer = np.array(self.weights) & self.get_divergence(
            mode="outside_of_trajectory"
        )

    # TODO: Add docstrings
    def get_orth_fit(self, **kwargs):
        """TODO."""
        kwargs["weighted"] = True  # include inner vals for orthogonal regression
        u, s = self.get_reads(**kwargs)
        a, b = np.sum(s * u), np.sum(u**2 - s**2)
        orth_beta = (b + ((b**2 + 4 * a**2) ** 0.5)) / (2 * a)
        return orth_beta

    # TODO: Add docstrings
    def get_orth_distx(self, orth_beta=None, **kwargs):
        """TODO."""
        if "weighted" not in kwargs:
            kwargs["weighted"] = "outer"
        u, s = self.get_reads(**kwargs)
        if orth_beta is None:
            orth_beta = self.get_orth_fit(**kwargs)
        s_real = np.array((s + (orth_beta * u)) / (1 + orth_beta**2))
        sdiff = np.array(s_real - s) / self.std_s
        udiff = np.array(orth_beta * s_real - u) / self.std_u * self.scaling
        return udiff**2 + sdiff**2

    # TODO: Add docstrings
    def get_pval(self, model="dynamical", **kwargs):
        """TODO."""
        # assuming var-scaled udiff, sdiff follow N(0,1),
        # the sum of errors for the cluster follows chi2(df=2n)
        if "weighted" not in kwargs:
            kwargs["weighted"] = "outer"
        distx = (
            self.get_orth_distx(**kwargs)
            if model == "orthogonal"
            else self.get_distx(**kwargs) / 2
        )
        return chi2.sf(df=2 * len(distx), x=np.sum(distx) / self.varx)

    def get_pval_diff_kinetics(self, orth_beta=None, min_cells=10, **kwargs):
        """Calculates the p-value for the likelihood ratio using the asymptotic property of the chi^2 distr.

        Derivation:
        - dists_dynamical and dists_orthogonal are squared N(0,1) distributed residuals
        - X1 = sum(dists_dynamical) / variance ~ chi2(df=2n)
        - X2 = sum(dists_orthogonal) / variance ~ chi2(df=2n)
        - Y1 = (X1 - df) / sqrt(2*df) ~ N(0,1) for large df
        - Y2 = (X2 - df) / sqrt(2*df) ~ N(0,1) for large df
        - since Y1~N(0,1) and Y2~N(0,1), Y1 - Y2 ~ N(0,2) or (Y1 -Y2) / sqrt(2) ~ N(0,1)
        - thus Z = (X1 - X2) / sqrt(4*df) ~ N(0,1) for large df

        Parameters
        ----------
        indices: "bool" array
            bool array for cluster of interest
        orth_beta: "float"
            orthogonal line fit beta

        Returns
        -------
        p-value
        """
        if (
            "weights_cluster" in kwargs
            and np.sum(kwargs["weights_cluster"]) < min_cells
        ):
            return 1
        if "weighted" not in kwargs:
            kwargs["weighted"] = "outer"
        distx = self.get_distx(**kwargs) / 2  # due to convolved assignments (tbd)
        orth_distx = self.get_orth_distx(orth_beta=orth_beta, **kwargs)
        denom = self.varx * np.sqrt(4 * 2 * len(distx))
        pval = norm.sf(
            (np.sum(distx) - np.sum(orth_distx)) / denom
        )  # see derivation above
        return pval

    # TODO: Add docstrings
    def get_cluster_mse(self, clusters=None, min_cells=10, weighted="outer"):
        """TODO."""
        if self.clusters is None or clusters is not None:
            self.initialize_diff_kinetics(clusters)
        mse = np.array(
            [
                self.get_mse(weights_cluster=self.clusters == c, weighted=weighted)
                for c in self.cats
            ]
        )
        if min_cells is not None:
            w = (
                self.weights_outer
                if weighted == "outer"
                else self.weights_upper
                if weighted == "upper"
                else self.weights
            )
            mse[
                np.array([np.sum(w & (self.clusters == c)) for c in self.cats])
                < min_cells
            ] = 0
        return mse

    # TODO: Add docstrings
    def get_cluster_pvals(self, clusters=None, model=None, orth_beta=None, **kwargs):
        """TODO."""
        if self.clusters is None or clusters is not None:
            self.initialize_diff_kinetics(clusters)
        pvals = np.array(
            [
                self.get_pval_diff_kinetics(
                    weights_cluster=self.clusters == c, orth_beta=orth_beta, **kwargs
                )
                if model is None
                else self.get_pval(
                    model=model, weights_cluster=self.clusters == c, **kwargs
                )
                for c in self.cats
            ]
        )
        return pvals

    # TODO: Add docstrings
    def differential_kinetic_test(
        self, clusters, as_df=None, min_cells=10, weighted="outer"
    ):
        """TODO."""
        # after fitting dyn. model
        self.initialize_diff_kinetics(clusters)
        mse = self.get_cluster_mse(weighted=weighted, min_cells=min_cells)

        weights_cluster = self.clusters == self.cats[np.argmax(mse)]
        self.orth_beta = self.get_orth_fit(
            weights_cluster=weights_cluster, weighted=False
        )  # include inner vals
        pval = self.get_pval_diff_kinetics(
            weights_cluster=weights_cluster, orth_beta=self.orth_beta, weighted=weighted
        )

        if pval > 1e-2:
            weighted = "upper"
            mse = self.get_cluster_mse(weighted=weighted, min_cells=min_cells)
            self.orth_beta = self.get_orth_fit(
                weights_cluster=self.clusters == self.cats[np.argmax(mse)]
            )

        self.pvals_kinetics = self.get_cluster_pvals(
            orth_beta=self.orth_beta, weighted=weighted
        )
        self.diff_kinetics = ",".join(
            [c for (c, p) in zip(self.cats, self.pvals_kinetics) if p < 1e-2]
        )
        if np.any(self.pvals_kinetics < 1e-2):
            self.pval_kinetics = np.max(self.pvals_kinetics[self.pvals_kinetics < 1e-2])

        if as_df:
            return pd.DataFrame(
                np.array(round(self.pvals_kinetics, 2))[None, :],
                columns=self.cats,
                index=["pval"],
            )


# TODO: Add docstrings
def get_reads(adata, key="fit", scaled=True, use_raw=False):
    """TODO."""
    if "Ms" not in adata.layers.keys():
        use_raw = True
    s = make_dense(adata.layers["spliced" if use_raw else "Ms"])
    u = make_dense(adata.layers["unspliced" if use_raw else "Mu"])
    if scaled:
        u /= adata.var[f"{key}_scaling"].values
    return u, s


# TODO: Add docstrings
def get_vars(adata, scaled=True, key="fit"):
    """TODO."""
    alpha = (
        adata.var[f"{key}_alpha"].values if f"{key}_alpha" in adata.var.keys() else 1
    )
    beta = adata.var[f"{key}_beta"].values if f"{key}_beta" in adata.var.keys() else 1
    gamma = adata.var[f"{key}_gamma"].values
    scaling = (
        adata.var[f"{key}_scaling"].values
        if f"{key}_scaling" in adata.var.keys()
        else 1
    )
    t_ = adata.var[f"{key}_t_"].values
    return alpha, beta * scaling if scaled else beta, gamma, scaling, t_


# TODO: Add docstrings
def get_latent_vars(adata, scaled=True, key="fit"):
    """TODO."""
    scaling = adata.var[f"{key}_scaling"].values
    std_u = adata.var[f"{key}_std_u"].values
    std_s = adata.var[f"{key}_std_s"].values
    u0 = adata.var[f"{key}_u0"].values
    s0 = adata.var[f"{key}_s0"].values
    pval_steady = (
        adata.var[f"{key}_pval_steady"].values
        if f"{key}_pval_steady" in adata.var.keys()
        else None
    )
    steady_u = (
        adata.var[f"{key}_steady_u"].values
        if f"{key}_steady_u" in adata.var.keys()
        else None
    )
    steady_s = (
        adata.var[f"{key}_steady_s"].values
        if f"{key}_steady_s" in adata.var.keys()
        else None
    )
    u0_scaled = (u0 / scaling if scaled else u0,)
    return std_u, std_s, u0_scaled, s0, pval_steady, steady_u, steady_s


# TODO: Add docstrings
def get_divergence(
    adata, mode="soft", use_latent_time=None, use_connectivities=None, **kwargs
):
    """TODO."""
    vdata = adata[:, ~np.isnan(adata.var["fit_alpha"].values)].copy()
    alpha, beta, gamma, scaling, t_ = get_vars(vdata)
    std_u, std_s, u0, s0, pval_steady, steady_u, steady_s = get_latent_vars(vdata)

    kwargs_ = {
        "kernel_width": None,
        "normalized": True,
        "var_scale": True,
        "reg_par": None,
        "min_confidence": 1e-2,
        "constraint_time_increments": False,
        "fit_steady_states": True,
        "fit_basal_transcription": None,
        "std_u": std_u,
        "std_s": std_s,
        "pval_steady": pval_steady,
        "steady_u": steady_u,
        "steady_s": steady_s,
    }
    kwargs_.update(adata.uns["recover_dynamics"])
    kwargs_.update(**kwargs)

    reg_time = None
    if use_latent_time is True:
        use_latent_time = "latent_time"
    if isinstance(use_latent_time, str) and use_latent_time in adata.obs.keys():
        reg_time = adata.obs[use_latent_time].values
    u, s = get_reads(vdata, use_raw=kwargs_["use_raw"])
    if kwargs_["fit_basal_transcription"]:
        u, s = u - u0, s - s0
    tau = (
        np.array(vdata.layers["fit_tau"]) if "fit_tau" in vdata.layers.keys() else None
    )
    tau_ = (
        np.array(vdata.layers["fit_tau_"])
        if "fit_tau_" in vdata.layers.keys()
        else None
    )

    kwargs_.update(
        {"t_": t_, "tau": tau, "tau_": tau_, "reg_time": reg_time, "mode": mode}
    )
    conn = get_connectivities(adata) if use_connectivities else None

    res = compute_divergence(
        u, s, alpha, beta, gamma, scaling, connectivities=conn, **kwargs_
    )
    return res
