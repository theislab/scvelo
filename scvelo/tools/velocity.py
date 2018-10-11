from ..logging import logg, settings
from ..preprocessing.moments import moments, second_order_moments
from .solver import solve_cov, solve2_inv, solve2_mle
from .utils import R_squared
from .velocity_confidence import velocity_confidence

from scipy.sparse import issparse
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

    n_counts = (adata.layers['unspliced'] > 0).sum(0)
    velocity_genes = n_counts.A1 > 10 if issparse(adata.layers['unspliced']) else n_counts > 10

    #lb_s, ub_s = np.percentile(Ms, [2, 98], axis=0)
    #lb_u, ub_u = np.percentile(Mu, [2, 98], axis=0)
    #idx = ((Ms <= lb_s) & (Mu < lb_u)) | ((Ms >= ub_s) & (Mu > ub_u))

    offset, gamma = solve_cov(Ms, Mu, fit_offset)

    stochastic = any([mode in item for item in ['stochastic', 'bayes', 'alpha']])
    if stochastic:
        res, tot = Mu - gamma[None, :] * Ms - offset[None, :], Mu - Mu.mean(0)
        r2 = R_squared(res, tot)
        velocity_genes = (r2 > .01) & (gamma > .01) & velocity_genes

        _Ms, _Mu = Ms[:, velocity_genes], Mu[:, velocity_genes]
        Mss, Mus = second_order_moments(adata[:, velocity_genes])
        offset2, gamma2 = solve_cov(2 * Mss - _Ms, 2 * Mus + _Mu, fit_offset2)
        #offset2, gamma2 = offset[velocity_genes], gamma[velocity_genes]

        # initialize covariance matrix
        res_std = res[:, velocity_genes].std(0)
        res2_std = (2 * Mus + _Mu - gamma2[None, :] * (2 * Mss - _Ms) - offset2[None, :]).std(0)

        # solve multiple regression
        offset[velocity_genes], offset2, gamma[velocity_genes] = \
            solve2_mle(_Ms, _Mu, Mus, Mss, fit_offset, fit_offset2) if mode == 'bayes' \
                else solve2_inv(_Ms, _Mu, 2 * Mss - _Ms, 2 * Mus + _Mu, res_std, res2_std, fit_offset, fit_offset2)

        _offset2 = np.zeros(adata.n_vars)
        _offset2[velocity_genes] = offset2
        adata.layers['variance_velocity'] = np.zeros(adata.shape)
        adata.layers['variance_velocity'][:, velocity_genes]\
            = beta[None, velocity_genes] * (2 * Mus - 2 * _Ms * _Mu + _Mu) - \
              gamma[None, velocity_genes] * (2 * Mss - 2 * _Ms ** 2 - _Ms) - \
              offset2 + 2 * offset[None, velocity_genes] * _Ms

        # if mode == 'alpha':
        #     Muu = second_order_moments_u(adata)
        #     offset2u, alpha = solve_cov(np.ones(Mu.shape) + 2 * Mu, 2 * Muu - Mu)
        #     pars.extend([offset2u, alpha])
        #     pars_str.extend(['_offset2u', '_alpha'])

    residual = Mu - gamma[None, :] * Ms - offset[None, :]
    residual[np.isnan(residual)] = 0
    adata.layers[vkey] = residual

    r2 = R_squared(residual=residual, total=Mu - Mu.mean(0))
    gamma[np.isnan(gamma)] = 0
    velocity_genes = velocity_genes & (r2 > .01) & (gamma > .01)

    pars = [offset, _offset2, beta, gamma, r2.round(2)] if stochastic else [offset, beta, gamma, r2.round(2)]
    pars_str = ['_offset', '_offset2', '_beta', '_gamma', '_r2'] if stochastic else ['_offset', '_beta', '_gamma', '_r2']
    for i, par in enumerate(pars_str): adata.var[vkey + par] = pars[i]

    if False:
        u, s = adata.layers['unspliced'].A, adata.layers['spliced'].A
        r2_raw = R_squared(residual=u - gamma[None, :] * s - offset[None, :], total=u - u.mean(0))
        adata.var[vkey + '_r2_raw'] = r2_raw
        velocity_genes = velocity_genes & (r2_raw > .01)

    adata.var['velocity_genes'] = velocity_genes
    if filter_genes: adata._inplace_subset_var(velocity_genes)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.layers`\n'
        '    \'' + vkey + '\', velocity vectors for each individual cell')

    velocity_confidence(adata, vkey)

    return adata if copy else None
