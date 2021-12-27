import random
from typing import List, Optional, Union

import hypothesis.strategies as st
from hypothesis import given
from hypothesis.extra.numpy import arrays

import numpy as np
from scipy.sparse import csr_matrix, issparse

from anndata import AnnData


# TODO: Add possibility to generate adata object with floats as counts
# TODO: Add possibility to generate different columns with different data types in
#       adata.obs and adata.var
@st.composite
def get_adata(
    draw,
    n_obs: Optional[int] = None,
    n_vars: Optional[int] = None,
    min_obs: Optional[int] = 1,
    max_obs: Optional[int] = 100,
    min_vars: Optional[int] = 1,
    max_vars: Optional[int] = 100,
    layer_keys: Optional[Union[List, str]] = None,
    min_layers: Optional[int] = 2,
    max_layers: int = 2,
    obsm_keys: Optional[Union[List, str]] = None,
    min_obsm: Optional[int] = 2,
    max_obsm: Optional[int] = 2,
    varm_keys: Optional[Union[List, str]] = None,
    min_varm: Optional[int] = 2,
    max_varm: Optional[int] = 2,
    obs_col_names=None,
    min_obs_cols=2,
    max_obs_cols=2,
    var_col_names=None,
    min_var_cols=2,
    max_var_cols=2,
    sparse_entries: bool = False,
) -> AnnData:
    """Generate an AnnData object.

    The largest possible value of a numerical entry is `1e5`.

    Arguments
    ---------
        n_obs
            Number of observations. If set to `None`, a random integer between `1` and
            `max_obs` will be drawn. Defaults to `None`.
        n_vars
            Number of variables. If set to `None`, a random integer between `1` and
            `max_vars` will be drawn. Defaults to `None`.
        min_obs
            Minimum number of observations. If set to `None`, there is no lower limit.
            Defaults to `1`.
        max_obs
            Maximum number of observations. If set to `None`, there is no upper limit.
            Defaults to `100`.
        min_vars
            Minimum number of variables. If set to `None`, there is no lower limit.
            Defaults to `1`.
        max_vars
            Maximum number of variables. If set to `None`, there is no upper limit.
            Defaults to `100`.
        layer_keys
            Names of layers. If set to `None`, layers will be named at random. Defaults
            to `None`.
        min_layers
            Minimum number of layers. Is set to the number of provided layer names if
            `layer_keys` is not `None`. Defaults to `2`.
        max_layers
            Maximum number of layers. Is set to the number of provided layer
            names if `layer_keys` is not `None`. Defaults to `2`.
        obsm_keys
            Names of multi-dimensional observations annotation. If set to `None`, names
            will be generated at random. Defaults to `None`.
        min_obsm
            Minimum number of multi-dimensional observations annotation. Is set to the
            number of keys if `obsm_keys` is not `None`. Defaults to `2`.
        max_obsm
            Maximum number of multi-dimensional observations annotation. Is set to the
            number of keys if `obsm_keys` is not `None`. Defaults to `2`.
        varm_keys
            Names of multi-dimensional variables annotation. If set to `None`, names
            will be generated at random. Defaults to `None`.
        min_varm
            Minimum number of multi-dimensional variables annotation. Is set to the
            number of keys if `varm_keys` is not `None`. Defaults to `2`.
        max_varm
            Maximum number of multi-dimensional variables annotation. Is set to the
            number of keys if `varm_keys` is not `None`. Defaults to `2`.
        obs_col_names
            Names of columns in `adata.obs`. If set to `None`, colums will be named at
            random. Defaults to `None`.
        min_obs_cols
            Minimum number of columns in `adata.obs`. Is set to the number of provided
            column names if `obs_col_names` is not `None`. Defaults to `2`.
        max_obs_cols
            Maximum number of columns in `adata.obs`. Is set to the number of provided
            column names if `obs_col_names` is not `None`. Defaults to `2`.
        var_col_names
            Names of columns in `adata.var`. If set to `None`, colums will be named at
            random. Defaults to `None`.
        min_var_cols
            Minimum number of columns in `adata.var`. Is set to the number of provided
            column names if `var_col_names` is not `None`. Defaults to `2`.
        max_var_cols
            Maximum number of columns in `adata.var`. Is set to the number of provided
            column names if `var_col_names` is not `None`. Defaults to `2`.
        sparse_entries
            Whether or not to make AnnData entries sparse.

    Returns
    -------
    AnnData
        Generated :class:`~anndata.AnnData` object.
    """

    if n_obs is None:
        n_obs = draw(st.integers(min_value=min_obs, max_value=max_obs))
    if n_vars is None:
        n_vars = draw(st.integers(min_value=min_vars, max_value=max_vars))

    if isinstance(layer_keys, str):
        layer_keys = [layer_keys]
    if isinstance(obsm_keys, str):
        obsm_keys = [obsm_keys]
    if isinstance(obs_col_names, str):
        obs_col_names = [obs_col_names]
    if isinstance(var_col_names, str):
        var_col_names = [var_col_names]

    if layer_keys is not None:
        min_layers = len(layer_keys)
        max_layers = len(layer_keys)
    if obsm_keys is not None:
        min_obsm = len(obsm_keys)
        max_obsm = len(obsm_keys)
    if varm_keys is not None:
        min_varm = len(varm_keys)
        max_varm = len(varm_keys)
    if obs_col_names is not None:
        min_obs_cols = len(obs_col_names)
        max_obs_cols = len(obs_col_names)
    if var_col_names is not None:
        min_var_cols = len(var_col_names)
        max_var_cols = len(var_col_names)

    X = draw(
        arrays(
            dtype=int,
            elements=st.integers(min_value=0, max_value=1e2),
            shape=(n_obs, n_vars),
        )
    )

    layers = draw(
        st.dictionaries(
            st.text(
                st.characters(
                    blacklist_categories=("Cs",),
                    blacklist_characters=("X"),
                ),
                min_size=1,
            )
            if layer_keys is None
            else st.sampled_from(layer_keys),
            arrays(
                dtype=int,
                elements=st.integers(min_value=0, max_value=1e2),
                shape=(n_obs, n_vars),
            ),
            min_size=min_layers,
            max_size=max_layers,
        )
    )

    obsm = draw(
        st.dictionaries(
            st.text(
                st.characters(
                    blacklist_categories=("Cs",),
                    blacklist_characters=("X"),
                ),
                min_size=1,
            )
            if obsm_keys is None
            else st.sampled_from(obsm_keys),
            arrays(
                dtype=int,
                elements=st.integers(min_value=0, max_value=1e2),
                shape=st.tuples(
                    st.integers(min_value=n_obs, max_value=n_obs),
                    st.integers(min_value=min_vars, max_value=max_vars),
                ),
            ),
            min_size=min_obsm,
            max_size=max_obsm,
        )
    )

    varm = draw(
        st.dictionaries(
            st.text(
                st.characters(
                    blacklist_categories=("Cs",),
                    blacklist_characters=("X"),
                ),
                min_size=1,
            )
            if varm_keys is None
            else st.sampled_from(varm_keys),
            arrays(
                dtype=int,
                elements=st.integers(min_value=0, max_value=1e2),
                shape=st.tuples(
                    st.integers(min_value=n_vars, max_value=n_vars),
                    st.integers(min_value=min_obs, max_value=max_obs),
                ),
            ),
            min_size=min_varm,
            max_size=max_varm,
        )
    )

    obs = draw(
        st.dictionaries(
            st.text(min_size=1)
            if obs_col_names is None
            else st.sampled_from(obs_col_names),
            st.lists(
                elements=st.integers(min_value=0, max_value=1e2),
                min_size=n_obs,
                max_size=n_obs,
            ),
            min_size=min_obs_cols,
            max_size=max_obs_cols,
        )
    )

    var = draw(
        st.dictionaries(
            st.text(min_size=1)
            if var_col_names is None
            else st.sampled_from(var_col_names),
            st.lists(
                elements=st.integers(min_value=0, max_value=1e2),
                min_size=n_vars,
                max_size=n_vars,
            ),
            min_size=min_var_cols,
            max_size=max_var_cols,
        )
    )

    # Make keys for layers and obsm unique
    for key in set(layers.keys()).intersection(obsm.keys()):
        layers[f"{key}_"] = layers.pop(key)

    if sparse_entries:
        layers = {key: csr_matrix(val) for key, val in layers.items()}
        obsm = {key: csr_matrix(val) for key, val in obsm.items()}
        varm = {key: csr_matrix(val) for key, val in varm.items()}
        return AnnData(
            X=csr_matrix(X), layers=layers, obsm=obsm, varm=varm, obs=obs, var=var
        )
    else:
        return AnnData(X=X, layers=layers, obsm=obsm, varm=varm, obs=obs, var=var)


class TestAdataGeneration:
    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_default_adata_generation(self, adata: AnnData):
        assert type(adata) is AnnData
        assert "X" not in adata.layers
        assert "X" not in adata.obsm
        assert "X" not in adata.varm

    @given(adata=get_adata(max_obs=5, max_vars=5, sparse_entries=True))
    def test_sparse_adata_generation(self, adata: AnnData):
        assert type(adata) is AnnData
        assert issparse(adata.X)
        assert np.all([issparse(adata.layers[layer]) for layer in adata.layers])
        assert np.all([issparse(adata.obsm[name]) for name in adata.obsm])
        assert np.all([issparse(adata.varm[name]) for name in adata.varm])

    @given(
        adata=get_adata(
            n_obs=2,
            n_vars=2,
            layer_keys=["unspliced", "spliced"],
            obsm_keys="X_umap",
            varm_keys=["varm_entry_1", "varm_entry_2"],
            obs_col_names=["louvain", "donor", "day"],
            var_col_names=["alpha", "beta", "gamma"],
        )
    )
    def test_custom_adata_generation(self, adata: AnnData):
        assert adata.X.shape == (2, 2)
        assert len(adata.layers) == 2
        assert len(adata.obsm) == 1
        assert len(adata.varm) == 2
        assert set(adata.layers.keys()) == {"unspliced", "spliced"}
        assert set(adata.obsm.keys()) == {"X_umap"}
        assert set(adata.varm.keys()) == {"varm_entry_1", "varm_entry_2"}
        assert set(adata.obs.columns) == {"louvain", "donor", "day"}
        assert set(adata.var.columns) == {"alpha", "beta", "gamma"}

    @given(adata=get_adata(max_obs=5, max_vars=5, min_obs_cols=0, max_obs_cols=10))
    def test_setting_number_obs_columns(self, adata):
        assert len(adata.obs.columns) >= 0
        assert len(adata.obs.columns) <= 10

    @given(adata=get_adata(max_obs=5, max_vars=5, min_var_cols=0, max_var_cols=10))
    def test_setting_number_var_columns(self, adata):
        assert len(adata.var.columns) >= 0
        assert len(adata.var.columns) <= 10


class TestBase:
    def _subset_modalities(
        self,
        adata: AnnData,
        n_modalities: int,
        from_layers: bool = True,
        from_obsm: bool = True,
    ):
        """Subset modalities of an AnnData object."""

        modalities = ["X"]
        if from_layers:
            modalities += list(adata.layers.keys())
        if from_obsm:
            modalities += list(adata.obsm.keys())
        return random.sample(modalities, min(len(modalities), n_modalities))

    def _subset_columns(
        self,
        adata: AnnData,
        n_cols: int,
        from_obs: bool = True,
        from_var: bool = True,
    ):
        """Subset columns of an AnnData object in `obs` and `var` slots."""

        columns = []
        if from_obs:
            columns += list(adata.obs.columns)
        if from_var:
            columns += list(adata.var.columns)
        return random.sample(columns, min(len(columns), n_cols))

    def _convert_to_float(self, adata: AnnData):
        """Convert AnnData entries in `layer` and `obsm` into floats."""

        for layer in adata.layers:
            adata.layers[layer] = adata.layers[layer].astype(float)
        for obs in adata.obsm:
            adata.obsm[obs] = adata.obsm[obs].astype(float)
