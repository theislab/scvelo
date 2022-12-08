import warnings

import numpy as np

from anndata import AnnData

from scvelo.core import invert, SplicingDynamics


# TODO: Add docstrings
# TODO Use `SplicingDynamics`
def unspliced(tau, u0, alpha, beta):
    """TODO."""
    expu = np.exp(-beta * tau)
    return u0 * expu + alpha / beta * (1 - expu)


# TODO: Add docstrings
def spliced(tau, s0, u0, alpha, beta, gamma):
    """TODO."""
    c = (alpha - u0 * beta) * invert(gamma - beta)
    expu, exps = np.exp(-beta * tau), np.exp(-gamma * tau)
    return s0 * exps + alpha / gamma * (1 - exps) + c * (exps - expu)


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


def simulation(
    n_obs=300,
    n_vars=None,
    alpha=None,
    beta=None,
    gamma=None,
    alpha_=None,
    t_max=None,
    noise_model="normal",
    noise_level=1,
    switches=None,
    random_seed=0,
):
    """Simulation of mRNA splicing kinetics.

    Simulated mRNA metabolism with transcription, splicing and degradation.
    The parameters for each reaction are randomly sampled from a log-normal distribution
    and time events follow the Poisson law. The total time spent in a transcriptional
    state is varied between two and ten hours.

    .. image:: https://user-images.githubusercontent.com/31883718/79432471-16c0a000-7fcc-11ea-8d62-6971bcf4181a.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    np.random.seed(random_seed)

    def draw_poisson(n):
        from random import seed, uniform  # draw from poisson

        seed(random_seed)
        t = np.cumsum([-0.1 * np.log(uniform(0, 1)) for _ in range(n - 1)])
        return np.insert(t, 0, 0)  # prepend t0=0

    def simulate_dynamics(tau, alpha, beta, gamma, u0, s0, noise_model, noise_level):
        ut, st = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=[u0, s0]
        ).get_solution(tau, stacked=False)
        if noise_model == "normal":  # add noise
            ut += np.random.normal(
                scale=noise_level * np.percentile(ut, 99) / 10, size=len(ut)
            )
            st += np.random.normal(
                scale=noise_level * np.percentile(st, 99) / 10, size=len(st)
            )
        ut, st = np.clip(ut, 0, None), np.clip(st, 0, None)
        return ut, st

    def simulate_gillespie(alpha, beta, gamma):
        # update rules:
        # transcription (u+1,s), splicing (u-1,s+1), degradation (u,s-1), nothing (u,s)
        update_rule = np.array([[1, 0], [-1, 1], [0, -1], [0, 0]])

        def update(props):
            if np.sum(props) > 0:
                props /= np.sum(props)
            p_cumsum = props.cumsum()
            p = np.random.rand()
            i = 0
            while p > p_cumsum[i]:
                i += 1
            return update_rule[i]

        u, s = np.zeros(len(alpha)), np.zeros(len(alpha))
        for i, alpha_i in enumerate(alpha):
            u_, s_ = (u[i - 1], s[i - 1]) if i > 0 else (0, 0)

            if (alpha_i == 0) and (u_ == 0) and (s_ == 0):
                du, ds = 0, 0
            else:
                du, ds = update(props=np.array([alpha_i, beta * u_, gamma * s_]))

            u[i], s[i] = (u_ + du, s_ + ds)
        return u, s

    alpha = 5 if alpha is None else alpha
    beta = 0.5 if beta is None else beta
    gamma = 0.3 if gamma is None else gamma
    alpha_ = 0 if alpha_ is None else alpha_

    t = draw_poisson(n_obs)
    if t_max is not None:
        t *= t_max / np.max(t)
    t_max = np.max(t)

    def cycle(array, n_vars=None):
        if isinstance(array, (np.ndarray, list, tuple)):
            return (
                array if n_vars is None else array * int(np.ceil(n_vars / len(array)))
            )
        else:
            return [array] if n_vars is None else [array] * n_vars

    # switching time point obtained as fraction of t_max rounded down
    switches = (
        cycle([0.4, 0.7, 1, 0.1], n_vars)
        if switches is None
        else cycle(switches, n_vars)
    )
    t_ = np.array([np.max(t[t < t_i * t_max]) for t_i in switches])

    noise_level = cycle(noise_level, len(switches) if n_vars is None else n_vars)

    n_vars = min(len(switches), len(noise_level)) if n_vars is None else n_vars
    U = np.zeros(shape=(len(t), n_vars))
    S = np.zeros(shape=(len(t), n_vars))

    def is_list(x):
        return isinstance(x, (tuple, list, np.ndarray))

    for i in range(n_vars):
        alpha_i = alpha[i] if is_list(alpha) and len(alpha) != n_obs else alpha
        beta_i = beta[i] if is_list(beta) and len(beta) != n_obs else beta
        gamma_i = gamma[i] if is_list(gamma) and len(gamma) != n_obs else gamma
        tau, alpha_vec, u0_vec, s0_vec = vectorize(
            t, t_[i], alpha_i, beta_i, gamma_i, alpha_=alpha_, u0=0, s0=0
        )

        if noise_model == "gillespie":
            U[:, i], S[:, i] = simulate_gillespie(alpha_vec, beta, gamma)
        else:
            U[:, i], S[:, i] = simulate_dynamics(
                tau,
                alpha_vec,
                beta_i,
                gamma_i,
                u0_vec,
                s0_vec,
                noise_model,
                noise_level[i],
            )

    if is_list(alpha) and len(alpha) == n_obs:
        alpha = np.nan
    if is_list(beta) and len(beta) == n_obs:
        beta = np.nan
    if is_list(gamma) and len(gamma) == n_obs:
        gamma = np.nan

    obs = {"true_t": t.round(2)}
    var = {
        "true_t_": t_[:n_vars],
        "true_alpha": np.ones(n_vars) * alpha,
        "true_beta": np.ones(n_vars) * beta,
        "true_gamma": np.ones(n_vars) * gamma,
        "true_scaling": np.ones(n_vars),
    }
    layers = {"unspliced": U, "spliced": S}

    return AnnData(S, obs, var, layers=layers)
