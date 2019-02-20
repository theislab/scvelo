import warnings
import numpy as np
exp = np.exp


"""Dynamics delineation and simulation"""


def unspliced(tau, u0, alpha, beta):
    expu = exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


def spliced(tau, s0, u0, alpha, beta, gamma):
    c = (alpha - u0 * beta) / (gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


def tau_u(u, u0, alpha, beta):
    def log(x):  # to avoid invalid values for log.
        return np.log(np.clip(x, 1e-6, 1 - 1e-6))
    u_ratio = (u - alpha / beta) / (u0 - alpha / beta)
    return - 1 / beta * log(u_ratio)


def tau_s(s, s0, u0, alpha, beta, gamma, u=None, tau=None, eps=1e-2):
    if tau is None:
        tau = tau_u(u, u0, alpha, beta) if u is not None else 1
    tau_prev, loss, n_iter, max_iter, mixed_states = 1e6, 1e6, 0, 10, np.any(alpha == 0)
    b0 = (alpha - beta * u0) / (gamma - beta)
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


def find_swichting_time(u, s, tau, o, alpha, beta, gamma):
    off, on = o == 0, o == 1
    if off.sum() > 0:
        u_, s_, tau_ = u[off], s[off], tau[off]
        c = (alpha / (gamma - beta) - alpha / gamma) * exp(-gamma * tau_)
        c_obs = (s_ - beta / (gamma - beta) * u_ + c)
        exp_t0_ = (c_obs * c).sum() / (c ** 2).sum()
        t0_ = -1 / gamma * np.log(exp_t0_) if exp_t0_ > 0 else np.max(tau[on])
    else:
        t0_ = np.max(tau)
    return t0_


def assign_timepoints(u, s, alpha, beta, gamma, t0_=None, u0_=None, s0_=None, n_timepoints=300):
    if t0_ is None:
        t0_ = tau_u(u0_, 0, alpha, beta)
    if u0_ is None or s0_ is None:
        u0_ = unspliced(t0_, 0, alpha, beta)
        s0_ = spliced(t0_, 0, 0, alpha, beta, gamma)

    tpoints = np.linspace(0, t0_, num=n_timepoints)
    tpoints_ = np.linspace(0, tau_u(np.min(u[s > 0]), u0_, 0, beta), num=n_timepoints)[1:]

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


def fit_alpha(u, s, beta, gamma, tau, o, fit_scaling=False):
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
    c_gamma = (1 - exps) / gamma + (exps - expu) / (gamma - beta)
    c_gamma_ = ((1 - exps0_) / gamma + (exps0_ - expu0_) / (gamma - beta)) * exps_ - (1 - expu0_) * (exps_ - expu_) / (gamma - beta)

    # concatenating together
    c = np.concatenate([c_beta, c_gamma, c_beta_, c_gamma_]).T
    x = np.concatenate([u[on], s[on], u[off], s[off]]).T
    alpha = (c * x).sum() / (c ** 2).sum()

    if fit_scaling:  # alternatively compute alpha and scaling simultaneously
        c = np.concatenate([c_gamma, c_gamma_]).T
        x = np.concatenate([s[on], s[off]]).T
        alpha_alt = (c * x).sum() / (c ** 2).sum()

        c = np.concatenate([c_beta, c_beta_]).T
        x = np.concatenate([u[on], u[off]]).T
        scaling = (c * x).sum() / (c ** 2).sum() / alpha_alt  # ~ alpha * z / alpha

    return alpha


def fit_scaling(u, t, t_, alpha, beta):
    tau, alpha, u0, _ = vectorize(t, t_, alpha, beta)
    ut = unspliced(tau, u0, alpha, beta)
    return (u * ut).sum() / (ut ** 2).sum()


def vectorize(t, t_, alpha, beta, gamma=None, alpha_=0, u0=0, s0=0):
    o = np.array(t < t_, dtype=bool)
    tau = t * o + (t - t_) * (1 - o)

    u0_ = unspliced(t_, u0, alpha, beta)
    s0_ = spliced(t_, s0, u0, alpha, beta, gamma if gamma is not None else beta / 2)

    # vectorize u0, s0 and alpha
    u0 = u0 * o + u0_ * (1 - o)
    s0 = s0 * o + s0_ * (1 - o)
    alpha = alpha * o + alpha_ * (1 - o)

    return tau, alpha, u0, s0


"""State-independent derivatives"""


# def du_du0(beta, tau):
#     return exp(-beta * tau)

# def ds_ds0(gamma, tau):
#     return exp(-gamma * tau)

# def ds_du0(beta, gamma, tau):
#     return - beta / (gamma - beta) * (exp(-gamma * tau) - exp(-beta * tau))

def dus_u0s0(tau, beta, gamma):
    du_u0 = exp(-beta * tau)
    ds_s0 = exp(-gamma * tau)
    ds_u0 = - beta / (gamma - beta) * (ds_s0 - du_u0)
    return du_u0, ds_s0, ds_u0


def dus_tau(tau, alpha, beta, gamma, u0=0, s0=0, du0_t0=0, ds0_t0=0):
    expu, exps, cb, cc = exp(-beta * tau), exp(-gamma * tau), alpha - beta * u0, alpha - gamma * s0
    du_tau = (cb - du0_t0) * expu
    ds_tau = (cc - ds0_t0) * exps - cb / (gamma - beta) * (gamma * exps - beta * expu) + du0_t0 * beta / (gamma - beta) * (exps - expu)
    return du_tau, ds_tau


def du(tau, alpha, beta, u0=0, du0=[0, 0, 0]):
    # du0 is the derivative du0 / d(alpha, beta, tau)
    expu, cb = exp(-beta * tau), alpha / beta
    du_a = du0[0] * expu + 1. / beta * (1 - expu)
    du_b = du0[1] * expu - cb / beta * (1 - expu) + (cb - u0) * tau * expu
    du_tau = (alpha - beta * u0 - du0[2]) * expu
    return du_a, du_b, du_tau


def ds(tau, alpha, beta, gamma, u0=0, s0=0, du0=[0, 0, 0], ds0=[0, 0, 0, 0]):
    # ds0 is the derivative ds0 / d(alpha, beta, gamma, tau)
    expu, exps, = exp(-beta * tau), exp(-gamma * tau)
    expus, cb, cc = exps - expu, alpha / beta, alpha / gamma

    cbu = (alpha - beta * u0) / (gamma - beta)
    ccu = (alpha - gamma * u0) / (gamma - beta)
    ccs = alpha / gamma - s0 - cbu

    ds_a = ds0[0] * exps + 1. / gamma * (1 - exps) + 1 / (gamma - beta) * (1 - beta * du0[0]) * expus
    ds_b = ds0[1] * exps + cbu * tau * expu + 1 / (gamma - beta) * (ccu - beta * du0[1]) * expus
    ds_c = ds0[2] * exps + ccs * tau * exps - cc / gamma * (1 - exps) - cbu / (gamma - beta) * expus

    ds_dtau = (alpha - gamma * s0 - ds0[3]) * exps - cbu * (gamma * exps - beta * expu) + du0[2] * beta / (gamma - beta) * (exps - expu)

    return ds_a, ds_b, ds_c, ds_dtau


def derivatives(u, s, t, t0_, alpha, beta, gamma, scaling=1, alpha_=0, u0=0, s0=0, update_switch=False, weights=None):
    o = np.array(t < t0_, dtype=int)

    du0 = np.array(du(t0_, alpha, beta, u0))[:, None] * (1 - o)[None, :]
    ds0 = np.array(ds(t0_, alpha, beta, gamma, u0, s0))[:, None] * (1 - o)[None, :]

    tau, alpha, u0, s0 = vectorize(t, t0_, alpha, beta, gamma, alpha_, u0, s0)
    udiff = np.array(unspliced(tau, u0, alpha, beta) * scaling - u)
    sdiff = np.array(spliced(tau, s0, u0, alpha, beta, gamma) - s)

    if weights is not None:
        udiff = np.multiply(udiff, weights)
        sdiff = np.multiply(sdiff, weights)

    # state-dependent derivatives:
    du_a, du_b, du_tau = du(tau, alpha, beta, u0, du0)
    du_a, du_b, du_tau = du_a * scaling, du_b * scaling, du_tau * scaling
    ds_a, ds_b, ds_c, ds_tau = ds(tau, alpha, beta, gamma, u0, s0, du0, ds0)

    dl_a = (du_a * (1 - o)).dot(udiff) + (ds_a * (1 - o)).dot(sdiff)
    dl_a_ = (du_a * o).dot(udiff) + (ds_a * o).dot(sdiff)

    dl_b = du_b.dot(udiff) + ds_b.dot(sdiff)
    dl_c = ds_c.dot(sdiff)

    if update_switch:
        dl_tau = du_tau * udiff + ds_tau * sdiff
        dl_t0_ = - du_tau.dot(udiff) - ds_tau.dot(sdiff)
    else:
        dl_tau, dl_t0_ = None, None

    return dl_a, dl_b, dl_c, dl_a_, dl_tau, dl_t0_
