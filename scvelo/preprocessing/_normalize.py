from typing import Optional, Set

import numpy as np
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs

from anndata import AnnData

from scvelo import logging as logg
from scvelo.core import get_modality, sum
from ._preprocessing_base import NormalizerBase


# TODO: Extend functionality to achieve same behavior as scv.pp.normalize_per_cell
class CountNormalizer(NormalizerBase):
    def __init__(
        self,
        inplace: bool = True,
        enforce: bool = False,
        use_initial_size: bool = True,
        target_count: Optional[float] = None,
    ):
        """Class for normalizing count data.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotated data should
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce normalizing modalities.
        use_initial_size
            Boolean flag to indicate whether to scale w.r.t. initial size per
            observation or not.
        target_count
            Number of counts per observation after normalization. If set to `None`, the
            counts per observation are normalized to the median of all counts per
            observation prior to normalization. Defaults to `None`.
        """

        super().__init__(inplace=inplace, enforce=enforce)
        self.use_initial_size = use_initial_size
        self.target_count = target_count

    def _transform(self, adata: AnnData, modalities: Set) -> Optional[AnnData]:
        """Normalize counts per observation.

        Arguments
        ---------
        adata
            Annotated data matrix of shape `n_obs` × `n_vars`.
        modalities
            Set of modalities to normalize.

        Returns
        -------
        Optional[AnnData]
            Returns annotated data if `self.inplace` is set to `False`, `None`
            otherwise.
        """

        for modality in modalities:
            count_data = get_modality(adata=adata, modality=modality)

            if self.use_initial_size:
                counts = self._get_initial_counts(adata=adata, modality=modality)
            else:
                counts = sum(count_data, axis=1)

            if self.target_count is None:
                target_count = np.median(counts)
            else:
                target_count = self.target_count

            if target_count != 0:
                counts = counts / target_count

            counts[counts == 0] = 1

            if issparse(count_data):
                sparsefuncs.inplace_row_scale(count_data, 1 / counts)
            else:
                count_data /= np.array(counts[:, None])

            logg.info(f"Normalized modality {modality}.")

        if "n_counts" not in adata.obs.columns:
            adata.obs["n_counts"] = sum(adata.X, axis=1)

        return adata if not self.inplace else None
