import numpy as np
exp = np.exp


"""Dynamics description and simulation"""


def unspliced(tau, u0, alpha, beta):
    expu = exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


def spliced(tau, s0, u0, alpha, beta, gamma):
    c = (alpha - u0 * beta) / (gamma - beta)
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


def vectorize(t, t_, alpha, beta, gamma, alpha_=0, u0=0, s0=0):
    o = np.array(t < t_, dtype=bool)
    tau = t * o + (t - t_) * (1 - o)

    u0_ = unspliced(t_, u0, alpha, beta)
    s0_ = spliced(t_, s0, u0, alpha, beta, gamma)

    # vectorize u0, s0 and alpha
    u0 = u0 * o + u0_ * (1 - o)
    s0 = s0 * o + s0_ * (1 - o)
    alpha = alpha * o + alpha_ * (1 - o)

    return tau, alpha, u0, s0


"""State-independent derivatives"""


def du_du0(beta, tau):
    return exp(-beta * tau)


def ds_ds0(gamma, tau):
    return exp(-gamma * tau)


def ds_du0(beta, gamma, tau):
    return - beta / (gamma - beta) * (exp(-gamma * tau) - exp(-beta * tau))


def du_dalpha(beta, tau):
    return 1 / beta * (1 - exp(-beta * tau))


def ds_dalpha(beta, gamma, tau):
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return 1 / gamma * (1 - exps) + 1 / (gamma - beta) * (exps - expu)


def du_dbeta(alpha, beta, tau, u0=0):
    expu = exp(-beta * tau)
    return (alpha / beta - u0) * tau * expu - alpha / beta ** 2 * (1 - expu)


def ds_dbeta(alpha, beta, gamma, tau, u0=0):
    cu = (alpha - beta * u0) / (gamma - beta)
    cs = (alpha - gamma * u0) / (gamma - beta) ** 2
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return cu * tau * expu + cs * (exps - expu)


def ds_dgamma(alpha, beta, gamma, tau, u0=0, s0=0):
    cu = (alpha - beta * u0) / (gamma - beta)
    cs = alpha / gamma - s0
    expu, exps = exp(-beta * tau), exp(-gamma * tau)
    return (cs - cu) * tau * exps - alpha / gamma ** 2 * (1 - exps) - cu / (gamma - beta) * (exps - expu)


def derivatives(alpha, beta, gamma, tau, o, diffU, diffS, alpha_=0, u0=0, s0=0, u0_=None, s0_=None, weights=None):
    # state-dependent derivatives:
    du_u = du_du0(beta, tau)
    ds_u = ds_du0(beta, gamma, tau)
    dL_du0 = (du_u * o).dot(diffU) + (ds_u * o).dot(diffS)
    dL_du0_ = (du_u * (1 - o)).dot(diffU) + (ds_u * (1 - o)).dot(diffS)

    ds_s = ds_ds0(gamma, tau)
    dL_ds0 = (ds_s * o).dot(diffS)
    dL_ds0_ = (ds_s * (1 - o)).dot(diffS)

    du_a = du_dalpha(beta, tau)
    ds_a = ds_dalpha(beta, gamma, tau)
    dL_dalpha = (du_a * o).dot(diffU) + (ds_a * o).dot(diffS)
    dL_dalpha_ = (du_a * (1 - o)).dot(diffU) + (ds_a * (1 - o)).dot(diffS)

    # independent derivatives: beta, gamma
    if u0_ is None: u0_ = alpha / beta
    if s0_ is None: s0_ = alpha / gamma

    u0 = u0 * o + u0_ * (1 - o)
    s0 = s0 * o + s0_ * (1 - o)

    alpha = alpha * o + alpha_ * (1 - o)

    dL_dbeta = du_dbeta(alpha, beta, tau, u0).dot(diffU) + ds_dbeta(alpha, beta, gamma, tau, u0).dot(diffS)
    dL_dgamma = ds_dgamma(alpha, beta, gamma, tau, u0, s0).dot(diffS)

    return dL_dalpha, dL_dbeta, dL_dgamma, dL_dalpha_, dL_du0, dL_ds0, dL_du0_, dL_ds0_


"""State-independent derivatives"""


def du_dalpha_o(beta, tau, state=1):
    expU = np.exp(-beta * tau)
    _expU = (1 - expU) * (state == 1) + expU * (state == 0)
    return 1 / beta * _expU


def du_dbeta_o(alpha, beta, tau, state=1, u0=0):
    expU = np.exp(-beta * tau)
    _expU = (1 - expU) * (state == 1) + expU * (state == 0)
    c = (alpha / beta - u0) * (state == 1) - alpha / beta * (state == 0)
    return c * tau * expU - alpha / beta ** 2 * _expU


def ds_dalpha_o(beta, gamma, tau, state=1):
    expU, expS = np.exp(-beta * tau), np.exp(-gamma * tau)
    _expU = (1 - expU) * (state == 1) + expU * (state == 0)
    _expS = (1 - expS) * (state == 1) + expS * (state == 0)
    return 1 / gamma * _expS - 1 / (gamma - beta) * (_expS - _expU)


def ds_dbeta_o(alpha, beta, gamma, tau, state=1, u0=0):
    _cS = (alpha - gamma * u0) / (gamma - beta) ** 2 * (state == 1) - alpha / (gamma - beta) ** 2 * (state == 0)
    _cU = (alpha - beta * u0) / (gamma - beta) * (state == 1) - alpha / (gamma - beta) * (state == 0)
    expU, expS = np.exp(-beta * tau), np.exp(-gamma * tau)
    return _cU * tau * expU + _cS * (expS - expU)


def ds_dgamma_o(alpha, beta, gamma, tau, state=1, u0=0, s0=0):
    _cU = (alpha - beta * u0) / (gamma - beta) * (state == 1) - alpha / (gamma - beta) * (state == 0)
    _cS = alpha / gamma - s0 * (state == 1) - alpha / gamma * (state == 0)
    expU, expS = np.exp(-beta * tau), np.exp(-gamma * tau)
    _expS = (1 - expS) * (state == 1) + expS * (state == 0)
    return (_cS - _cU) * tau * expS - alpha / gamma ** 2 * _expS - _cU / (gamma - beta) * (expS - expU)


def derivatives_o(alpha, beta, gamma, tau, o, diffU, diffS, alpha_=0, u0=0, s0=0, u0_=None, s0_=None):
    dL_dalpha = du_dalpha_o(beta, tau, o).dot(diffU) + ds_dalpha_o(beta, gamma, tau, o).dot(diffS)
    dL_dbeta = du_dbeta_o(alpha, beta, tau, o).dot(diffU) + ds_dbeta_o(alpha, beta, gamma, tau, o).dot(diffS)
    dL_dgamma = ds_dgamma_o(alpha, beta, gamma, tau, o).dot(diffS)

    return dL_dalpha, dL_dbeta, dL_dgamma, 0, 0, 0, 0, 0
