from .utils import prod_sum_obs
from scipy.optimize import minimize
import numpy as np
import warnings


def solve_cov(x, y, fit_offset=False):
    """Solution to least squares: gamma = cov(X,Y) / var(X)
    """
    n_obs, n_var = x.shape
    if fit_offset:
        cov_xy, cov_xx = prod_sum_obs(x, y) / n_obs, prod_sum_obs(x, x) / n_obs
        mean_x, mean_y = x.mean(0), y.mean(0)
        numerator_offset = cov_xx * mean_y - cov_xy * mean_x
        numerator_gamma = cov_xy - mean_x * mean_y
        offset, gamma = (numerator_offset, numerator_gamma) / (cov_xx - mean_x * mean_x)
    else:
        offset, gamma = np.zeros(n_var), prod_sum_obs(x, y) / prod_sum_obs(x, x)
    return offset, gamma


def solve2_inv(x, y, x2, y2, res_std=None, res2_std=None, fit_offset=False, fit_offset2=False):
    """Solution to the 2-dim generalized least squares: gamma = inv(X'QX)X'QY
    """
    n_obs, n_var = x.shape
    offset, offset_ss = np.zeros(n_var, dtype="float32"), np.zeros(n_var, dtype="float32")
    gamma = np.ones(n_var, dtype="float32")

    if (res_std is None) or (res2_std is None): res_std, res2_std = np.ones(n_var), np.ones(n_var)
    ones, zeros = np.ones(n_obs), np.zeros(n_obs)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x, y = np.vstack((x/res_std, x2/res2_std)), np.vstack((y/res_std, y2/res2_std))

    if fit_offset and fit_offset2:
        for i in range(n_var):
            A = np.c_[np.vstack((np.c_[ones/res_std[i], zeros], np.c_[zeros, ones/res2_std[i]])), x[:, i]]
            offset[i], offset_ss[i], gamma[i] = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))
    elif fit_offset:
        for i in range(n_var):
            A = np.c_[np.hstack((ones/res_std[i], zeros)), x[:, i]]
            offset[i], gamma[i] = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))
    elif fit_offset2:
        for i in range(n_var):
            A = np.c_[np.hstack((zeros, ones/res2_std[i])), x[:, i]]
            offset_ss[i], gamma[i] = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))
    else:
        for i in range(n_var):
            A = np.c_[x[:, i]]
            gamma[i] = np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(y[:, i]))

    return offset, offset_ss, gamma


def solve2_mle(Ms, Mu, Mus, Mss, fit_offset=False, fit_offset2=False):
    """Maximizing the log likelihood using weights according to empirical bayes
    """
    n_obs, n_var = Ms.shape
    offset, offset_ss = np.zeros(n_var, dtype="float32"), np.zeros(n_var, dtype="float32")
    gamma = np.ones(n_var, dtype="float32")

    def sse(A, data, b):
        sigma = (A.dot(data) - b).std(1)
        return np.log(sigma).sum()  # np.log(np.sqrt(2*np.pi) * sigma).sum() + (.5 * (res/sigma[:, None])**2).sum()

    if fit_offset and fit_offset2:
        for i in range(n_var):
            data = np.vstack((Mu[:, i], Ms[:, i], Mus[:, i], Mss[:, i]))
            offset[i], offset_ss[i], gamma[i] = \
                minimize(lambda m: sse(np.array([[1, -m[2], 0, 0], [1, m[2], 2, -2 * m[2]]]),
                                       data, b=np.array(m[0], m[1])), x0=(1e-4, 1e-4, 1), method="L-BFGS-B").x
    elif fit_offset:
        for i in range(n_var):
            data = np.vstack((Mu[:, i], Ms[:, i], Mus[:, i], Mss[:, i]))
            offset[i], gamma[i] = \
                minimize(lambda m: sse(np.array([[1, -m[1], 0, 0], [1, m[1], 2, -2 * m[1]]]),
                                       data, b=np.array(m[0], 0)), x0=(1e-4, 1), method="L-BFGS-B").x
    elif fit_offset2:
        for i in range(n_var):
            data = np.vstack((Mu[:, i], Ms[:, i], Mus[:, i], Mss[:, i]))
            offset_ss[i], gamma[i] = \
                minimize(lambda m: sse(np.array([[1, -m[1], 0, 0], [1, m[1], 2, -2 * m[1]]]),
                                       data, b=np.array(0, m[0])), x0=(1e-4, 1), method="L-BFGS-B").x
    else:
        for i in range(n_var):
            data = np.vstack((Mu[:, i], Ms[:, i], Mus[:, i], Mss[:, i]))
            gamma[i] = \
                minimize(lambda m: sse(np.array([[1, -m, 0, 0], [1, m, 2, -2 * m]]), data, b=0),
                         x0=gamma[i], method="L-BFGS-B").x
    return offset, offset_ss, gamma
