from .dynamical_model_utils import tau_u, inv, unspliced, spliced, vectorize
import warnings
import numpy as np
exp = np.exp


def tau_s(s, s0, u0, alpha, beta, gamma, u=None, tau=None, eps=1e-2):
    if tau is None:
        tau = tau_u(u, u0, alpha, beta) if u is not None else 1
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
        t0_ = tau_u(u0_, 0, alpha, beta)
    if u0_ is None or s0_ is None:
        u0_, s0_ = (unspliced(t0_, 0, alpha, beta), spliced(t0_, 0, 0, alpha, beta, gamma))

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
    # du_tau = (alpha - beta * u0 - du0[2]) * expuå
    return du_a, du_b


def ds(tau, alpha, beta, gamma, u0=0, s0=0, du0=[0, 0, 0], ds0=[0, 0, 0, 0], dtau=[0, 0, 0]):
    # ds0 is the derivative ds0 / d(alpha, beta, gamma, tau)
    expu, exps, = exp(-beta * tau), exp(-gamma * tau)
    expus = exps - expu

    cbu = (alpha - beta * u0) * inv(gamma - beta)
    ccu = (alpha - gamma * u0) * inv(gamma - beta)
    ccs = alpha / gamma - s0 - cbu

    # dsu0 = ds0 * exps - beta / (gamma - beta) * du0

    ds_a = ds0[0] * exps + 1. / gamma * (1 - exps) + 1 * inv(gamma - beta) * (1 - beta * du0[0]) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[0]
    ds_b = ds0[1] * exps + cbu * tau * expu + 1 * inv(gamma - beta) * (ccu - beta * du0[1]) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[1]
    ds_c = ds0[2] * exps + ccs * tau * exps - alpha / gamma**2 * (1 - exps) - cbu * inv(gamma - beta) * expus + (ccs * gamma * exps + cbu * beta * expu) * dtau[2]

    # ds_dtau = (alpha - gamma * s0 - ds0[3]) * exps - cbu * (gamma * exps - beta * expu) + du0[2] * beta * inv(gamma - beta) * (exps - expu)

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

    # dl_tau = du_tau * udiff + ds_tau * sdiff
    # dl_t0_ = - du_tau.dot(udiff) - ds_tau.dot(sdiff)
    dl_tau, dl_t0_ = None, None

    return dl_a, dl_b, dl_c, dl_a_, dl_tau, dl_t0_
