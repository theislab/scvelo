from ..logging import logg, settings
from ..preprocessing.moments import moments, second_order_moments
from .solver import solve_cov, solve2_inv, solve2_mle
from .utils import R_squared
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class Velocity:
    def __init__(self, adata=None, Ms=None, Mu=None, residual=None):
        self._adata = adata
        self._Ms = adata.layers['Ms'] if Ms is None else Ms
        self._Mu = adata.layers['Mu'] if Mu is None else Mu

        n_obs, n_vars = self._Ms.shape
        self._residual, self._residual2 = residual, None
        self._offset, self._offset2 = np.zeros(n_vars, dtype=np.float32), np.zeros(n_vars, dtype=np.float32)
        self._gamma, self._r2 = np.zeros(n_vars, dtype=np.float32), np.zeros(n_vars, dtype=np.float32)
        self._beta, self._velocity_genes = np.ones(n_vars, dtype=np.float32), np.ones(n_vars, dtype=bool)

    def compute_deterministic(self, fit_offset=False):
        self._offset, self._gamma = solve_cov(self._Ms, self._Mu, fit_offset)
        self._residual = self._Mu - self._gamma * self._Ms
        if fit_offset: self._residual -= self._offset

        self._r2 = R_squared(self._residual, total=self._Mu - self._Mu.mean(0))
        self._velocity_genes = (self._r2 > .01) & (self._gamma > .01)

    def compute_stochastic(self, fit_offset=False, fit_offset2=False, mode=None):
        if self._residual is None: self.compute_deterministic(fit_offset)
        idx = self._velocity_genes
        is_subset = True if len(set(idx)) > 1 else False

        _adata = self._adata[:, idx] if is_subset else self._adata
        _Ms = self._Ms[:, idx] if is_subset else self._Ms
        _Mu = self._Mu[:, idx] if is_subset else self._Mu
        _residual = self._residual[:, idx] if is_subset else self._residual

        _Mss, _Mus = second_order_moments(_adata)

        var_ss = 2 * _Mss - _Ms
        cov_us = 2 * _Mus + _Mu

        _offset2, _gamma2 = solve_cov(var_ss, cov_us, fit_offset2)

        # initialize covariance matrix
        res_std = _residual.std(0)
        res2_std = (cov_us - _gamma2 * var_ss - _offset2).std(0)

        # solve multiple regression
        self._offset[idx], self._offset2[idx], self._gamma[idx] = \
            solve2_mle(_Ms, _Mu, _Mus, _Mss, fit_offset, fit_offset2) if mode == 'bayes' \
                else solve2_inv(_Ms, _Mu, var_ss, cov_us, res_std, res2_std, fit_offset, fit_offset2)

        self._residual = self._Mu - self._gamma * self._Ms
        if fit_offset: self._residual -= self._offset

        _residual2 = (cov_us - 2 * _Ms * _Mu) - self._gamma[idx] * (var_ss - 2 * _Ms ** 2)
        if fit_offset: _residual2 += 2 * self._offset[idx] * _Ms
        if fit_offset2: _residual2 -= self._offset2[idx]
        if is_subset:
            self._residual2 = np.zeros(self._Ms.shape, dtype=np.float32)
            self._residual2[:, idx] = _residual2
        else:
            self._residual2 = _residual2

        # if mode == 'alpha':
        #     Muu = second_order_moments_u(adata)
        #     offset2u, alpha = solve_cov(np.ones(Mu.shape) + 2 * Mu, 2 * Muu - Mu)
        #     pars.extend([offset2u, alpha])
        #     pars_str.extend(['_offset2u', '_alpha'])

    def get_residuals(self):
        return self._residual, self._residual2

    def get_pars(self):
        return self._offset, self._offset2, self._beta, self._gamma, self._r2, self._velocity_genes

    def get_pars_names(self):
        return ['_offset', '_offset2', '_beta', '_gamma', '_r2', '_genes']


def velocity(data, vkey='velocity', mode=None, fit_offset=False, fit_offset2=False, filter_genes=False, copy=False):
    """Estimates velocities in a gene-specific manner

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
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
    adata = data.copy() if copy else data
    if 'Ms' not in adata.layers.keys(): moments(adata)

    logg.info('computing velocities', r=True)
    velo = Velocity(adata)
    velo.compute_deterministic(fit_offset)

    stochastic = any([mode is not None and mode in item for item in ['stochastic', 'bayes', 'alpha']])
    if stochastic:
        if filter_genes and len(set(velo._velocity_genes)) > 1:
            idx = velo._velocity_genes
            adata._inplace_subset_var(idx)
            velo = Velocity(adata, residual=velo._residual[:, idx])
        velo.compute_stochastic(fit_offset, fit_offset2, mode)
        adata.layers[vkey], adata.layers['variance_' + vkey] = velo.get_residuals()
    else:
        adata.layers[vkey], _ = velo.get_residuals()

    pars = velo.get_pars()
    for i, key in enumerate(velo.get_pars_names()):
        if len(set(pars[i])) > 1: adata.var[vkey + key] = pars[i]

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + vkey + '\', velocity vectors for each individual cell (adata.layers)')

    if filter_genes and len(set(velo._velocity_genes)) > 1:  # re-initialize after filtering
        adata._inplace_subset_var(velo._velocity_genes)

    return adata if copy else None
