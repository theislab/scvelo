from abc import abstractmethod
from functools import reduce
from typing import Dict, List, Optional, Set, Union

import numpy as np
from numpy import ndarray
from pandas import DataFrame
from scipy.sparse import issparse, spmatrix

from anndata import AnnData

import scvelo.logging as logg
from scvelo.core import get_modality, ScveloBase, sum


# TODO: Add `fit` and `fit_transform`
class PreprocessingBase(ScveloBase):
    def __init__(self, inplace: bool = True, enforce: bool = False):
        """Base class for preprocessing of data to find protein velocity.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotaded data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce preprocessing step.
        """

        super().__init__(inplace=inplace)
        self.enforce = enforce

    def __call__(self, data, modalities, **kwargs):
        return self.transform(data, modalities, **kwargs)

    def _is_normalized(
        self,
        adata: AnnData,
        modality: str,
        n_cells: Optional[int] = 5,
        atol: float = 1e-3,
    ) -> Optional[bool]:
        """Check if a modality is not normalized.

        Arguments
        ---------
        adata
            Annotated data to check for normalization.
        modality
            Name of modality to check.
        n_cells
            Number of cells to check normalization on.
        atol
            Absolute tolerance for deciding whether modality is normalized or not.

        Returns
        -------
        Optional[bool]
            `True` if modality looks normalized, `False` otherwise.
        """

        data_matrix = get_modality(adata=adata, modality=modality)
        if issparse(data_matrix):
            data_matrix = data_matrix[:n_cells].data
        else:
            data_matrix = data_matrix[:n_cells]

        if np.allclose(data_matrix % 1, 0, atol=atol):
            return False
        else:
            return True

    def _remove_normalized_modalities(
        self, adata: AnnData, modalities: Set, task: str = None
    ):
        """Remove normalized modalities.

        Arguments
        ---------
        adata
            Annotated data object.
        modalities
            Set of modalities.
        task
            Preprocessing task.
        """

        normalized_modalities = set()

        for modality in modalities:
            if self._is_normalized(adata, modality):
                if task is None:
                    logg.warn(
                        f"Modality `{modality}` is excluded as it looks already "
                        "normalized."
                    )
                else:
                    logg.warn(
                        f"Modality `{modality}` is excluded from {task} as it looks "
                        "already normalized."
                    )
                normalized_modalities.add(modality)

        modalities -= normalized_modalities

    def _remove_not_normalized_modalities(
        self, adata: AnnData, modalities: Set, task: str = None
    ):
        """Remove not normalized modalities.

        Arguments
        ---------
        adata
            Annotated data object.
        modalities
            Set of modalities.
        task
            Preprocessing task.
        """

        not_normalized_modalities = set()
        for modality in modalities:
            if not self._is_normalized(adata, modality):
                if task is None:
                    logg.warn(
                        f"Modality `{modality}` is excluded as it seems not yet "
                        "normalized."
                    )
                else:
                    logg.warn(
                        f"Modality `{modality}` is excluded from {task} as it seems "
                        "not yet normalized."
                    )
                not_normalized_modalities.add(modality)

        modalities -= not_normalized_modalities

    def _set_initial_counts(self, adata: AnnData):
        """Set initial counts of modalities.

        Arguments
        ---------
        adata
            Annotated data object.
        modalities
            Set of modalities.
        """

        for modality in {"X"} | {*adata.layers} | {*adata.obsm}:
            if modality == "X":
                if "initial_size" not in adata.obs.columns:
                    adata.obs["initial_size"] = sum(
                        get_modality(adata=adata, modality=modality), axis=1
                    )
            elif f"initial_size_{modality}" not in adata.obs.columns:
                adata.obs[f"initial_size_{modality}"] = sum(
                    get_modality(adata=adata, modality=modality), axis=1
                )

    # TODO: Handle case where initial counts have not been set
    def _get_initial_counts(self, adata: AnnData, modality: str) -> ndarray:
        """Retrieve initial counts per observation.

        Arguments
        ---------
        adata
            Annotated data object.
        modality
            Name of modality for which to retrieve initial counts.

        Returns
        -------
        ndarray
            Initial counts of specified modality.
        """

        if modality == "X":
            return adata.obs["initial_size"].values
        else:
            return adata.obs[f"initial_size_{modality}"].values

    @abstractmethod
    def _transform(
        self, adata: AnnData, modalities: Set, **kwargs
    ) -> Optional[AnnData]:
        pass

    def _update_modalities(self, adata: AnnData, modalities: Set):
        pass

    # TODO: Allow more general input type for `data`
    def transform(
        self, data: AnnData, modalities: Optional[Union[List, str]] = None, **kwargs
    ) -> Optional[AnnData]:
        """Transform data

        Arguments
        ---------
        adata
            Annotated data object.
        modality
            Name of modality for which to retrieve initial counts.

        Returns
        -------
        Optional[AnnData]
            Transformed annotated data is `inplace` was set to `False`, `None`
            otherwise.
        """

        adata = data.copy() if not self.inplace else data
        if isinstance(modalities, str):
            modalities = [modalities]

        if modalities is None:
            modalities = {"X"}
        else:
            # TODO: Add warning if provided modalities do not exist
            modalities = set(modalities) & ({"X"} | {*adata.layers} | {*adata.obsm})

        if not self.enforce:
            self._update_modalities(adata=adata, modalities=modalities)

        self._set_initial_counts(adata=adata)

        adata = self._transform(adata=adata, modalities=modalities, **kwargs)

        return adata if not self.inplace else None


class FilterBase(PreprocessingBase):
    def __init__(
        self,
        inplace: bool = True,
        enforce: bool = False,
        keep: Optional[Dict[str, List[Union[str, int]]]] = None,
    ):
        """Base class for filtering of annotated data.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotaded data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce preprocessing step.
        """

        super().__init__(inplace=inplace, enforce=enforce)

    def _subset_var(self, adata: AnnData, modality: str, filter_mask: ndarray):
        """Subset variables using a filter mask.

        Arguments
        ---------
        data
            Annotated data to subset from.
        modality
            Modality to subset.
        filter_mask
            Boolean ndarray indicating which variables to subset.
        """

        if modality in adata.obsm.keys():
            if isinstance(adata.obsm[modality], DataFrame):
                adata.obsm[modality] = adata.obsm[modality].loc[:, filter_mask]
            else:
                adata.obsm[modality] = adata.obsm[modality][:, filter_mask]
        else:
            adata._inplace_subset_var(filter_mask)


class CountFilterBase(FilterBase):
    def __init__(
        self,
        inplace: bool = True,
        enforce: bool = False,
        lower_bound: Optional[int] = None,
        upper_bound: Optional[int] = None,
        shared_counts: bool = False,
    ):
        """Base class for filtering variables based on their counts.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotaded data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce filtering modalities.
        lower_bound
            Minimum number of occurances to not be filtered out.
        upper_bound
            Maximum number of occurances to not be filtered out.
        shared_counts
            Boolearn flag to filter by shared counts.
        """

        super().__init__(inplace=inplace, enforce=enforce)

        self.lower_bound = lower_bound if lower_bound is not None else -np.inf
        self.upper_bound = upper_bound if upper_bound is not None else np.inf
        self.shared_counts = shared_counts

    def _update_modalities(self, adata: AnnData, modalities: Set):
        """Remove already normalized modalities

        Arguments
        ---------
        adata
            Annotated data to applt count filter to.
        modalities
            Candidate modalities for filtering.
        """

        self._remove_normalized_modalities(
            adata=adata, modalities=modalities, task="filtering"
        )

    def _sum_shared_counts(
        self, adata: AnnData, modalities: Set
    ) -> Union[ndarray, spmatrix]:
        """Sum counts across modalities.

        If a count entry is trivial, i.e. `0`, the shared counts will be set to `0`.

        Arguments
        ---------
        adata
            Annotated data from which shared counts are calculated.
        modalities
            Set of modalities for which to calculate shared counts.

        Returns
        -------
        Union[ndarray, spmatrix]
            Summed counts across modalities.
        """

        count_matrices = [
            get_modality(adata=adata, modality=modality) for modality in modalities
        ]

        nonzeros = reduce(
            lambda x, y: (x > 0).multiply(y > 0) if issparse(x) else (x > 0) * (y > 0),
            count_matrices,
        )

        if issparse(nonzeros):
            return nonzeros.multiply(np.array(count_matrices).sum(axis=0))
        else:
            return nonzeros * np.array(count_matrices).sum(axis=0)

    def get_filter_mask(
        self, adata: AnnData, modality: ndarray, vars_to_keep: List[Union[int, str]]
    ) -> ndarray:
        """Create filter for genes of a given modality.

        Arguments
        ---------
        adata
            Annotated data to apply count filter to.
        modality
            Modality to find filter mask for.
        vars_to_keep
            List of variables to exclude from filtering.

        Returns
        -------
        ndarray
            Mask to filter out genes.
        """

        summed_counts = sum(modality, axis=0)
        filter_mask = (self.lower_bound <= summed_counts) & (
            summed_counts <= self.upper_bound
        )

        # `len(vars_to_keep) > 0` assures that `vars_to_keep` is not empty and that
        # `isinstance(vars_to_keep[0], str)`, thus, does not fail
        if len(vars_to_keep) > 0 and isinstance(vars_to_keep[0], str):
            vars_to_keep = [
                adata.var.index.get_loc(variable) for variable in vars_to_keep
            ]

        filter_mask[vars_to_keep] = True

        return filter_mask
