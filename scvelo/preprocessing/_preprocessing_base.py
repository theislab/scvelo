from abc import abstractmethod
from typing import List, Optional, Set, Union

import numpy as np
from numpy import ndarray
from scipy.sparse import issparse

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
            Boolean flag to indicate whether operations on annotated data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce preprocessing step.
        """

        super().__init__(inplace=inplace)
        self.enforce = enforce

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
