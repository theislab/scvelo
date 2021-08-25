import warnings

import numpy as np

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import LinearRegression
from scvelo.preprocessing.moments import moments, second_order_moments
from .optimization import leastsq_generalized, maximum_likelihood
from .utils import groups_to_bool, make_dense, R_squared, strings_to_categoricals

warnings.simplefilter(action="ignore", category=FutureWarning)


class Velocity:
    def __init__(
        self,
        adata=None,
        Ms=None,
        Mu=None,
        groups_for_fit=None,
        groupby=None,
        residual=None,
        constrain_ratio=None,
        min_r2=0.01,
        min_ratio=0.01,
        use_highly_variable=True,
        r2_adjusted=True,
        use_raw=False,
    ):
        self._adata = adata
        self._Ms, self._Mu = Ms, Mu
        if Ms is None:
            self._Ms = adata.layers["spliced"] if use_raw else adata.layers["Ms"]
        if Mu is None:
            self._Mu = adata.layers["unspliced"] if use_raw else adata.layers["Mu"]
        self._Ms, self._Mu = make_dense(self._Ms), make_dense(self._Mu)

        n_obs, n_vars = self._Ms.shape
        self._residual, self._residual2 = residual, None
        self._offset = np.zeros(n_vars, dtype=np.float32)
        self._offset2 = np.zeros(n_vars, dtype=np.float32)
        self._gamma = np.zeros(n_vars, dtype=np.float32)
        self._qreg_ratio = np.zeros(n_vars, dtype=np.float32)
        self._r2 = np.zeros(n_vars, dtype=np.float32)
        self._beta = np.ones(n_vars, dtype=np.float32)
        self._velocity_genes = np.ones(n_vars, dtype=bool)
        self._groups_for_fit = groups_to_bool(adata, groups_for_fit, groupby)
        self._constrain_ratio = constrain_ratio
        self._r2_adjusted = r2_adjusted
        self._min_r2 = min_r2
        self._min_ratio = min_ratio
        self._highly_variable = None
        if use_highly_variable is not None and adata is not None:
            if "highly_variable" in adata.var.keys():
                self._highly_variable = adata.var["highly_variable"]

    def compute_deterministic(self, fit_offset=False, perc=None):
        subset = self._groups_for_fit
        Ms = self._Ms if subset is None else self._Ms[subset]
        Mu = self._Mu if subset is None else self._Mu[subset]

        lr = LinearRegression(fit_intercept=fit_offset, percentile=perc)
        lr.fit(Ms, Mu)
        self._offset = lr.intercept_
        self._gamma = lr.coef_

        if self._constrain_ratio is not None:
            if np.size(self._constrain_ratio) < 2:
                self._constrain_ratio = [None, self._constrain_ratio]
            cr = self._constrain_ratio
            self._gamma = np.clip(self._gamma, cr[0], cr[1])

        self._residual = self._Mu - self._gamma * self._Ms
        if fit_offset:
            self._residual -= self._offset
        _residual = self._residual

        # velocity genes
        if self._r2_adjusted:
            lr = LinearRegression(fit_intercept=fit_offset)
            lr.fit(Ms, Mu)
            _offset = lr.intercept_
            _gamma = lr.coef_

            _residual = self._Mu - _gamma * self._Ms
            if fit_offset:
                _residual -= _offset

        self._qreg_ratio = np.array(self._gamma)  # quantile regression ratio

        self._r2 = R_squared(_residual, total=self._Mu - self._Mu.mean(0))
        self._velocity_genes = (
            (self._r2 > self._min_r2)
            & (self._gamma > self._min_ratio)
            & (np.max(self._Ms > 0, 0) > 0)
            & (np.max(self._Mu > 0, 0) > 0)
        )

        if self._highly_variable is not None:
            self._velocity_genes &= self._highly_variable

        if np.sum(self._velocity_genes) < 2:
            min_r2 = np.percentile(self._r2, 80)
            self._velocity_genes = self._r2 > min_r2
            min_r2 = np.round(min_r2, 4)
            logg.warn(
                f"You seem to have very low signal in splicing dynamics.\n"
                f"The correlation threshold has been reduced to {min_r2}.\n"
                f"Please be cautious when interpreting results."
            )

    def compute_stochastic(
        self, fit_offset=False, fit_offset2=False, mode=None, perc=None
    ):
        if self._residual is None:
            self.compute_deterministic(fit_offset=fit_offset, perc=perc)

        idx = np.ones(self._velocity_genes.shape, dtype=bool)
        if np.any(self._velocity_genes):
            idx = self._velocity_genes
        is_subset = len(set(idx)) > 1

        _adata = self._adata[:, idx] if is_subset else self._adata
        _Ms = self._Ms[:, idx] if is_subset else self._Ms
        _Mu = self._Mu[:, idx] if is_subset else self._Mu
        _residual = self._residual[:, idx] if is_subset else self._residual

        _Mss, _Mus = second_order_moments(_adata)

        var_ss = 2 * _Mss - _Ms
        cov_us = 2 * _Mus + _Mu

        lr = LinearRegression(fit_intercept=fit_offset2)
        lr.fit(var_ss, cov_us)
        _offset2 = lr.intercept_
        _gamma2 = lr.coef_

        # initialize covariance matrix
        res_std = _residual.std(0)
        res2_std = (cov_us - _gamma2 * var_ss - _offset2).std(0)

        # solve multiple regression
        self._offset[idx], self._offset2[idx], self._gamma[idx] = (
            maximum_likelihood(_Ms, _Mu, _Mus, _Mss, fit_offset, fit_offset2)
            if mode == "bayes"
            else leastsq_generalized(
                _Ms,
                _Mu,
                var_ss,
                cov_us,
                res_std,
                res2_std,
                fit_offset,
                fit_offset2,
                perc,
            )
        )

        self._residual = self._Mu - self._gamma * self._Ms
        if fit_offset:
            self._residual -= self._offset

        _residual2 = (cov_us - 2 * _Ms * _Mu) - self._gamma[idx] * (
            var_ss - 2 * _Ms ** 2
        )
        if fit_offset:
            _residual2 += 2 * self._offset[idx] * _Ms
        if fit_offset2:
            _residual2 -= self._offset2[idx]
        if is_subset:
            self._residual2 = np.zeros(self._Ms.shape, dtype=np.float32)
            self._residual2[:, idx] = _residual2
        else:
            self._residual2 = _residual2

    def get_pars(self):
        return (
            self._offset,
            self._offset2,
            self._beta,
            self._gamma,
            self._qreg_ratio,
            self._r2,
            self._velocity_genes,
        )

    def get_pars_names(self):
        return [
            "_offset",
            "_offset2",
            "_beta",
            "_gamma",
            "_qreg_ratio",
            "_r2",
            "_genes",
        ]


def write_residuals(adata, vkey, residual=None, cell_subset=None):
    if residual is not None:
        if cell_subset is None:
            adata.layers[vkey] = residual
        else:
            if vkey not in adata.layers.keys():
                adata.layers[vkey] = np.zeros(adata.shape, dtype=np.float32)
            adata.layers[vkey][cell_subset] = residual


def write_pars(adata, vkey, pars, pars_names, add_key=None):
    for i, key in enumerate(pars_names):
        key = f"{vkey}{key}_{add_key}" if add_key is not None else f"{vkey}{key}"
        if len(set(pars[i])) > 1:
            adata.var[key] = pars[i]
        elif key in adata.var.keys():
            del adata.var[key]


def velocity(
    data,
    vkey="velocity",
    mode="stochastic",
    fit_offset=False,
    fit_offset2=False,
    filter_genes=False,
    groups=None,
    groupby=None,
    groups_for_fit=None,
    constrain_ratio=None,
    use_raw=False,
    use_latent_time=None,
    perc=[5, 95],
    min_r2=1e-2,
    min_likelihood=1e-3,
    r2_adjusted=None,
    use_highly_variable=True,
    diff_kinetics=None,
    copy=False,
    **kwargs,
):
    """Estimates velocities in a gene-specific manner.

    The steady-state model [Manno18]_ determines velocities by quantifying how
    observations deviate from a presumed steady-state equilibrium ratio of unspliced to
    spliced mRNA levels. This steady-state ratio is obtained by performing a linear
    regression restricting the input data to the extreme quantiles. By including
    second-order moments, the stochastic model [Bergen19]_ exploits not only the balance
    of unspliced to spliced mRNA levels but also their covariation. By contrast, the
    likelihood-based dynamical model [Bergen19]_ solves the full splicing kinetics and
    generalizes RNA velocity estimation to transient systems. It is also capable of
    capturing non-observed steady states.

    .. image:: https://user-images.githubusercontent.com/31883718/69636491-ff057100-1056-11ea-90b7-d04098112ce1.png

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name under which to refer to the computed velocities
        for `velocity_graph` and `velocity_embedding`.
    mode: `'deterministic'`, `'stochastic'` or `'dynamical'` (default: `'stochastic'`)
        Whether to run the estimation using the steady-state/deterministic,
        stochastic or dynamical model of transcriptional dynamics.
        The dynamical model requires to run `tl.recover_dynamics` first.
    fit_offset: `bool` (default: `False`)
        Whether to fit with offset for first order moment dynamics.
    fit_offset2: `bool`, (default: `False`)
        Whether to fit with offset for second order moment dynamics.
    filter_genes: `bool` (default: `True`)
        Whether to remove genes that are not used for further velocity analysis.
    groups: `str`, `list` (default: `None`)
        Subset of groups, e.g. [‘g1’, ‘g2’, ‘g3’],
        to which velocity analysis shall be restricted.
    groupby: `str`, `list` or `np.ndarray` (default: `None`)
        Key of observations grouping to consider.
    groups_for_fit: `str`, `list` or `np.ndarray` (default: `None`)
        Subset of groups, e.g. [‘g1’, ‘g2’, ‘g3’],
        to which steady-state fitting shall be restricted.
    constrain_ratio: `float` or tuple of type `float` or None: (default: `None`)
        Bounds for the steady-state ratio.
    use_raw: `bool` (default: `False`)
        Whether to use raw data for estimation.
    use_latent_time: `bool`or `None` (default: `None`)
        Whether to use latent time as a regularization for velocity estimation.
    perc: `float` (default: `[5, 95]`)
        Percentile, e.g. 98, for extreme quantile fit.
    min_r2: `float` (default: 0.01)
        Minimum threshold for coefficient of determination
    min_likelihood: `float` (default: `None`)
        Minimal likelihood for velocity genes to fit the model on.
    r2_adjusted: `bool` (default: `None`)
        Whether to compute coefficient of determination
        on full data fit (adjusted) or extreme quantile fit (None)
    use_highly_variable: `bool` (default: True)
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    velocity: `.layers`
        velocity vectors for each individual cell
    velocity_genes, velocity_beta, velocity_gamma, velocity_r2: `.var`
        parameters
    """  # noqa E501

    adata = data.copy() if copy else data
    if not use_raw and "Ms" not in adata.layers.keys():
        moments(adata)

    logg.info("computing velocities", r=True)

    strings_to_categoricals(adata)

    if mode is None or (mode == "dynamical" and "fit_alpha" not in adata.var.keys()):
        mode = "stochastic"
        logg.warn(
            "Falling back to stochastic model. "
            "For the dynamical model run tl.recover_dynamics first."
        )

    if mode in {"dynamical", "dynamical_residuals"}:
        from .dynamical_model_utils import get_divergence, get_reads, get_vars

        gene_subset = ~np.isnan(adata.var["fit_alpha"].values)
        vdata = adata[:, gene_subset]
        alpha, beta, gamma, scaling, t_ = get_vars(vdata)

        connect = not adata.uns["recover_dynamics"]["use_raw"]
        kwargs_ = {
            "kernel_width": None,
            "normalized": True,
            "var_scale": True,
            "reg_par": None,
            "min_confidence": 1e-2,
            "constraint_time_increments": False,
            "fit_steady_states": True,
            "fit_basal_transcription": None,
            "use_connectivities": connect,
            "time_connectivities": connect,
            "use_latent_time": use_latent_time,
        }
        kwargs_.update(adata.uns["recover_dynamics"])
        kwargs_.update(**kwargs)

        if "residuals" in mode:
            u, s = get_reads(vdata, use_raw=adata.uns["recover_dynamics"]["use_raw"])
            if kwargs_["fit_basal_transcription"]:
                u, s = u - adata.var["fit_u0"], s - adata.var["fit_s0"]
            o = vdata.layers["fit_t"] < t_
            vt = u * beta - s * gamma  # ds/dt
            wt = (alpha * o - beta * u) * scaling  # du/dt
        else:
            vt, wt = get_divergence(vdata, mode="velocity", **kwargs_)

        vgenes = adata.var.fit_likelihood > min_likelihood
        if min_r2 is not None:
            if "fit_r2" not in adata.var.keys():
                velo = Velocity(
                    adata,
                    groups_for_fit=groups_for_fit,
                    groupby=groupby,
                    constrain_ratio=constrain_ratio,
                    min_r2=min_r2,
                    use_highly_variable=use_highly_variable,
                    use_raw=use_raw,
                )
                velo.compute_deterministic(fit_offset=fit_offset, perc=perc)
                adata.var["fit_r2"] = velo._r2
            vgenes &= adata.var.fit_r2 > min_r2

        lb, ub = np.nanpercentile(adata.var.fit_scaling, [10, 90])
        vgenes = (
            vgenes
            & (adata.var.fit_scaling > np.min([lb, 0.03]))
            & (adata.var.fit_scaling < np.max([ub, 3]))
        )

        adata.var[f"{vkey}_genes"] = vgenes

        adata.layers[vkey] = np.ones(adata.shape) * np.nan
        adata.layers[vkey][:, gene_subset] = vt

        adata.layers[f"{vkey}_u"] = np.ones(adata.shape) * np.nan
        adata.layers[f"{vkey}_u"][:, gene_subset] = wt

        if filter_genes and len(set(vgenes)) > 1:
            adata._inplace_subset_var(vgenes)

    elif mode in {"steady_state", "deterministic", "stochastic"}:
        categories = (
            adata.obs[groupby].cat.categories
            if groupby is not None and groups is None and groups_for_fit is None
            else [None]
        )

        for cat in categories:
            groups = cat if cat is not None else groups

            cell_subset = groups_to_bool(adata, groups, groupby)
            _adata = adata if groups is None else adata[cell_subset]
            velo = Velocity(
                _adata,
                groups_for_fit=groups_for_fit,
                groupby=groupby,
                constrain_ratio=constrain_ratio,
                min_r2=min_r2,
                r2_adjusted=r2_adjusted,
                use_highly_variable=use_highly_variable,
                use_raw=use_raw,
            )
            velo.compute_deterministic(fit_offset=fit_offset, perc=perc)

            if mode == "stochastic":
                if filter_genes and len(set(velo._velocity_genes)) > 1:
                    adata._inplace_subset_var(velo._velocity_genes)
                    residual = velo._residual[:, velo._velocity_genes]
                    _adata = adata if groups is None else adata[cell_subset]
                    velo = Velocity(
                        _adata,
                        residual=residual,
                        groups_for_fit=groups_for_fit,
                        groupby=groupby,
                        constrain_ratio=constrain_ratio,
                        use_highly_variable=use_highly_variable,
                    )
                velo.compute_stochastic(fit_offset, fit_offset2, mode, perc=perc)

            write_residuals(adata, vkey, velo._residual, cell_subset)
            write_residuals(adata, f"variance_{vkey}", velo._residual2, cell_subset)
            write_pars(adata, vkey, velo.get_pars(), velo.get_pars_names(), add_key=cat)

            if filter_genes and len(set(velo._velocity_genes)) > 1:
                adata._inplace_subset_var(velo._velocity_genes)

    else:
        raise ValueError(
            "Mode can only be one of these: deterministic, stochastic or dynamical."
        )

    if f"{vkey}_genes" in adata.var.keys() and np.sum(adata.var[f"{vkey}_genes"]) < 10:
        logg.warn(
            "Too few genes are selected as velocity genes. "
            "Consider setting a lower threshold for min_r2 or min_likelihood."
        )

    if diff_kinetics:
        if not isinstance(diff_kinetics, str):
            diff_kinetics = "fit_diff_kinetics"
        if diff_kinetics in adata.var.keys():
            if diff_kinetics in adata.uns["recover_dynamics"]:
                groupby = adata.uns["recover_dynamics"]["fit_diff_kinetics"]
            else:
                groupby = "clusters"
            clusters = adata.obs[groupby]
            for i, v in enumerate(np.array(adata.var[diff_kinetics].values, dtype=str)):
                if len(v) > 0 and v != "nan":
                    idx = 1 - clusters.isin([a.strip() for a in v.split(",")])
                    adata.layers[vkey][:, i] *= idx
                    if mode == "dynamical":
                        adata.layers[f"{vkey}_u"][:, i] *= idx

    adata.uns[f"{vkey}_params"] = {"mode": mode, "fit_offset": fit_offset, "perc": perc}

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        f"    '{vkey}', velocity vectors for each individual cell (adata.layers)"
    )

    return adata if copy else None


def velocity_genes(
    data,
    vkey="velocity",
    min_r2=0.01,
    min_ratio=0.01,
    use_highly_variable=True,
    copy=False,
):
    """Estimates velocities in a gene-specific manner

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name under which to refer to the computed velocities.
    min_r2: `float` (default: 0.01)
        Minimum threshold for coefficient of determination
    min_ratio: `float` (default: 0.01)
        Minimum threshold for quantile regression un/spliced ratio.
    use_highly_variable: `bool` (default: True)
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Updates `adata` attributes
    velocity_genes: `.var`
        genes to be used for further velocity analysis (velocity graph and embedding)
    """

    adata = data.copy() if copy else data
    if f"{vkey}_genes" not in adata.var.keys():
        velocity(adata, vkey)
    vgenes = np.ones(adata.n_vars, dtype=bool)

    if "Ms" in adata.layers.keys() and "Mu" in adata.layers.keys():
        vgenes &= np.max(adata.layers["Ms"] > 0, 0) > 0
        vgenes &= np.max(adata.layers["Mu"] > 0, 0) > 0

    if min_r2 is not None and f"{vkey}_r2" in adata.var.keys():
        vgenes &= adata.var[f"{vkey}_r2"] > min_r2

    if min_ratio is not None and f"{vkey}_qreg_ratio" in adata.var.keys():
        vgenes &= adata.var[f"{vkey}_qreg_ratio"] > min_ratio

    if use_highly_variable and "highly_variable" in adata.var.keys():
        vgenes &= adata.var["highly_variable"].values

    if np.sum(vgenes) < 2:
        logg.warn(
            "You seem to have very low signal in splicing dynamics.\n"
            "Consider reducing the thresholds and be cautious with interpretations.\n"
        )

    adata.var[f"{vkey}_genes"] = vgenes

    logg.info("Number of obtained velocity_genes:", np.sum(adata.var[f"{vkey}_genes"]))

    return adata if copy else None
