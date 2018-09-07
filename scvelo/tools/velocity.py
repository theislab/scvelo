from ..logging import logg, settings
from ..preprocessing.moments import moments, second_order_moments
from .solver import solve_cov, solve2_inv, solve2_mle
from .utils import prod_sum_obs
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def velocity(adata, vkey='velocity', mode='stochastic', fit_offset=False, fit_offset2=False, filter_genes=False, copy=False):
    """Estimates velocities in a gene-specific manner

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name under which to refer to the computed velocities for `velocity_graph` and `velocity_embedding`.
    mode: `'deterministic'`, `'stochastic'` or `'bayes'` (default: `'stochastic'`)
        Whether to run the estimation using the deterministic or stochastic model of transcriptional dynamics.
        `'bayes'` solves the stochastic model and accounts for heteroscedasticity, but is slower than `'stochastic'`.
    fit_offset: `bool` (default: `False`)
        Whether to fit with offset for first order moment dynamics.
    fit_offset2: `bool`, (default: `False`)
        Whether to fit with offset for second order moment dynamics.
    filter_genes: `bool` (default: `True`)
        Whether to remove genes that are not used for further velocity analysis.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity: `.layers`
        velocity vectors for each individual cell
    variance_velocity: `.layers`
        velocity vectors for the cell variances
    velocity_offset, velocity_beta, velocity_gamma, velocity_r2: `.var`
        parameters
    """
    if 'Ms' not in adata.layers.keys(): moments(adata)
    Ms, Mu = adata.layers['Ms'], adata.layers['Mu']

    logg.info('computing velocities', r=True)

    beta = np.ones(adata.n_vars, dtype="float32")  # estimate all rates in units of the splicing rate
    offset, gamma = solve_cov(Ms, Mu, fit_offset)

    stochastic = any([mode in item for item in ['stochastic', 'bayes', 'alpha']])
    if stochastic:
        Mss, Mus = second_order_moments(adata)
        offset2, gamma2 = solve_cov(2 * Mss - Ms, 2 * Mus + Mu, fit_offset2)

        # initialize covariance matrix
        res_std = (Mu - gamma[None, :] * Ms - offset[None, :]).std(0)
        res2_std = (2 * Mus + Mu - gamma2[None, :] * (2 * Mss - Ms) - offset2[None, :]).std(0)

        # solve multiple regression
        offset, offset2, gamma = solve2_mle(Ms, Mu, Mus, Mss, fit_offset, fit_offset2) if mode == 'bayes' \
            else solve2_inv(Ms, Mu, 2 * Mss - Ms, 2 * Mus + Mu, res_std, res2_std, fit_offset, fit_offset2)

        adata.layers['variance_velocity'] = beta[None, :] * (2 * Mus - 2 * Ms * Mu + Mu) - \
                                    gamma[None, :] * (2 * Mss - 2 * Ms ** 2 - Ms) - \
                                    offset2[None, :] + 2 * offset[None, :] * Ms

        # if mode == 'alpha':
        #     Muu = second_order_moments_u(adata)
        #     offset2u, alpha = solve_cov(np.ones(Mu.shape) + 2 * Mu, 2 * Muu - Mu)
        #     pars.extend([offset2u, alpha])
        #     pars_str.extend(['_offset2u', '_alpha'])

    res, tot = Mu - gamma[None, :] * Ms - offset[None, :], Mu - Mu.mean(0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = np.ones(adata.n_vars) - prod_sum_obs(res, res) / prod_sum_obs(tot, tot)

    gamma[np.isnan(gamma)], r2[np.isnan(r2)], res[np.isnan(res)] = 0, 0, 0

    pars = [offset, offset2, beta, gamma, r2] if stochastic else [offset, beta, gamma, r2]
    pars_str = ['_offset', '_offset2', '_beta', '_gamma', '_r2'] if stochastic else ['_offset', '_beta', '_gamma', '_r2']
    for i, par in enumerate(pars_str): adata.var[vkey + par] = pars[i]

    adata.layers[vkey] = Mu - gamma[None, :] * Ms - offset[None, :]

    adata.var['velocity_genes'] = (r2 > .01) & (gamma > .01)
    if filter_genes: adata._inplace_subset_var(adata.var['velocity_genes'])

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.layers`\n'
        '    \'' + vkey + '\', velocity vectors for each individual cell')

    return adata if copy else None
