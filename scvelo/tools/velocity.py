from .. import settings
from .. import logging as logg
from ..preprocessing.moments import moments, second_order_moments
from .optimization import leastsq_NxN, leastsq_generalized, maximum_likelihood
from .utils import R_squared, groups_to_bool, make_dense, strings_to_categoricals

import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class Velocity:
    def __init__(self, adata=None, Ms=None, Mu=None, groups_for_fit=None, groupby=None, residual=None, use_raw=False,):
        self._adata = adata
        self._Ms = adata.layers['spliced'] if use_raw else adata.layers['Ms'] if Ms is None else Ms
        self._Mu = adata.layers['unspliced'] if use_raw else adata.layers['Mu'] if Mu is None else Mu
        self._Ms, self._Mu = make_dense(self._Ms), make_dense(self._Mu)

        n_obs, n_vars = self._Ms.shape
        self._residual, self._residual2 = residual, None
        self._offset, self._offset2 = np.zeros(n_vars, dtype=np.float32), np.zeros(n_vars, dtype=np.float32)
        self._gamma, self._r2 = np.zeros(n_vars, dtype=np.float32), np.zeros(n_vars, dtype=np.float32)
        self._beta, self._velocity_genes = np.ones(n_vars, dtype=np.float32), np.ones(n_vars, dtype=bool)
        self._groups_for_fit = groups_to_bool(adata, groups_for_fit, groupby)

    def compute_deterministic(self, fit_offset=False, perc=None):
        Ms = self._Ms if self._groups_for_fit is None else self._Ms[self._groups_for_fit]
        Mu = self._Mu if self._groups_for_fit is None else self._Mu[self._groups_for_fit]

        self._offset, self._gamma = leastsq_NxN(Ms, Mu, fit_offset, perc)
        self._residual = self._Mu - self._gamma * self._Ms
        if fit_offset: self._residual -= self._offset

        self._r2 = R_squared(self._residual, total=self._Mu - self._Mu.mean(0))
        self._velocity_genes = (self._r2 > .01) & (self._gamma > .01)

    def compute_stochastic(self, fit_offset=False, fit_offset2=False, mode=None, perc=None):
        if self._residual is None: self.compute_deterministic(fit_offset=fit_offset, perc=perc)
        idx = self._velocity_genes
        is_subset = True if len(set(idx)) > 1 else False

        _adata = self._adata[:, idx] if is_subset else self._adata
        _Ms = self._Ms[:, idx] if is_subset else self._Ms
        _Mu = self._Mu[:, idx] if is_subset else self._Mu
        _residual = self._residual[:, idx] if is_subset else self._residual

        _Mss, _Mus = second_order_moments(_adata)

        var_ss = 2 * _Mss - _Ms
        cov_us = 2 * _Mus + _Mu

        _offset2, _gamma2 = leastsq_NxN(var_ss, cov_us, fit_offset2)

        # initialize covariance matrix
        res_std = _residual.std(0)
        res2_std = (cov_us - _gamma2 * var_ss - _offset2).std(0)

        # solve multiple regression
        self._offset[idx], self._offset2[idx], self._gamma[idx] = \
            maximum_likelihood(_Ms, _Mu, _Mus, _Mss, fit_offset, fit_offset2) if mode == 'bayes' \
                else leastsq_generalized(_Ms, _Mu, var_ss, cov_us, res_std, res2_std, fit_offset, fit_offset2, perc)

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
        #     offset2u, alpha = leastsq_NxN(np.ones(Mu.shape) + 2 * Mu, 2 * Muu - Mu)
        #     pars.extend([offset2u, alpha])
        #     pars_str.extend(['_offset2u', '_alpha'])

    def get_pars(self):
        return self._offset, self._offset2, self._beta, self._gamma, self._r2, self._velocity_genes

    def get_pars_names(self):
        return ['_offset', '_offset2', '_beta', '_gamma', '_r2', '_genes']


def write_residuals(adata, vkey, residual=None, cell_subset=None):
    if residual is not None:
        if cell_subset is None:
            adata.layers[vkey] = residual
        else:
            if vkey not in adata.layers.keys(): adata.layers[vkey] = np.zeros(adata.shape, dtype=np.float32)
            adata.layers[vkey][cell_subset] = residual


def write_pars(adata, vkey, pars, pars_names, add_key=None):
    for i, key in enumerate(pars_names):
        key = vkey + key + '_' + add_key if add_key is not None else vkey + key
        if len(set(pars[i])) > 1: adata.var[key] = pars[i]


def velocity(data, vkey='velocity', mode=None, fit_offset=False, fit_offset2=False, filter_genes=False,
             groups=None, groupby=None, groups_for_fit=None, use_raw=False, perc=None, copy=False):
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
    groups: `str`, `list` (default: `None`)
        Subset of groups, e.g. [‘g1’, ‘g2’, ‘g3’], to which velocity analysis shall be restricted.
    groupby: `str`, `list` or `np.ndarray` (default: `None`)
        Key of observations grouping to consider.
    groups_for_fit: `str`, `list` or `np.ndarray` (default: `None`)
        Subset of groups, e.g. [‘g1’, ‘g2’, ‘g3’], to which steady-state fitting shall be restricted.
    use_raw: `bool` (default: `False)
        Whether to use raw data for estimation.
    perc: `float` (default: `None`)
        Percentile, e.g. 98, upon which velocity estimation is performed (for robustness).
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
    if not use_raw and 'Ms' not in adata.layers.keys(): moments(adata)

    logg.info('computing velocities', r=True)

    strings_to_categoricals(adata)
    categories = adata.obs[groupby].cat.categories \
        if groupby is not None and groups is None and groups_for_fit is None else [None]

    for cat in categories:
        groups = cat if cat is not None else groups

        cell_subset = groups_to_bool(adata, groups, groupby)
        _adata = adata if groups is None else adata[cell_subset]

        velo = Velocity(_adata, groups_for_fit=groups_for_fit, groupby=groupby, use_raw=use_raw)
        velo.compute_deterministic(fit_offset, perc=perc)

        if any([mode is not None and mode in item for item in ['stochastic', 'bayes', 'alpha']]):
            if filter_genes and len(set(velo._velocity_genes)) > 1:
                adata._inplace_subset_var(velo._velocity_genes)
                residual = velo._residual[:, velo._velocity_genes]
                _adata = adata if groups is None else adata[cell_subset]
                velo = Velocity(_adata, residual=residual, groups_for_fit=groups_for_fit, groupby=groupby)
            velo.compute_stochastic(fit_offset, fit_offset2, mode, perc=perc)

        write_residuals(adata, vkey, velo._residual, cell_subset)
        write_residuals(adata, 'variance_' + vkey, velo._residual2, cell_subset)
        write_pars(adata, vkey, velo.get_pars(), velo.get_pars_names(), add_key=cat)

        if filter_genes and len(set(velo._velocity_genes)) > 1:
            adata._inplace_subset_var(velo._velocity_genes)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added \n' 
              '    \'' + vkey + '\', velocity vectors for each individual cell (adata.layers)')

    return adata if copy else None
