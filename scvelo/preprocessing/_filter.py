from typing import Dict, List, Optional, Set

from typing_extensions import Literal

import numpy as np
from numpy import ndarray

from anndata import AnnData

import scvelo.logging as logg
from scvelo.core import get_modality
from ._preprocessing_base import CountFilterBase


def _log_progress(
    kind: Literal["cell_count", "gene_count"],
    filter_mask: ndarray,
    modality: str,
    lower_bound: float,
    upper_bound: float,
) -> None:
    """Log filter progress.

    Arguments
    ---------
    kind
        Kind of filter applied.
    filter_mask
        Filter mask applied.
    modality
        Modality being filtered.
    lower_bound
        Lower bound for variable to be filtered.
    upper_bound
        Upper bound for variable to be filtered.
    """

    n_filtered_genes = np.invert(filter_mask).sum()

    # Final log: Filtered out X genes CONNECTING_PHRASE less / more than Y FILTER_BASE
    if kind == "cell_count":
        connecting_phrase = "detected in"
        filter_base = "cells"
    elif kind == "gene_count":
        connecting_phrase = "with"
        filter_base = "counts"

    if n_filtered_genes > 0:
        logg.info(f"Filtered out {n_filtered_genes} genes ", end=" ")

        if not np.isinf(lower_bound) and not np.isinf(upper_bound):
            logg.info(
                f"{connecting_phrase} less than {lower_bound} and more than "
                f"{upper_bound} {filter_base} ({modality}).",
                no_indent=True,
            )
        elif np.isinf(upper_bound):
            logg.info(
                f"{connecting_phrase} less than {lower_bound} {filter_base} "
                "({modality}).",
                no_indent=True,
            )
        elif np.isinf(lower_bound):
            logg.info(
                f"{connecting_phrase} more than {lower_bound} {filter_base} "
                "({modality}).",
                no_indent=True,
            )


class CellCountFilter(CountFilterBase):
    def __init__(
        self,
        inplace: bool = True,
        enforce: bool = False,
        min_n_obs: Optional[int] = None,
        max_n_obs: Optional[int] = None,
        shared_counts: bool = False,
    ):
        """Class to filter variables based on number of observations expressing them.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotaded data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce filtering modalities.
        min_n_obs
            Minimum number of observations in which a variable is expressed to not be
            filtered out.
        max_n_obs
            Maximum number of observations in which a variable is expressed to not be
            filtered out.
        shared_counts
            Boolearn flag to filter by shared counts
        """

        super().__init__(
            inplace=inplace,
            enforce=enforce,
            lower_bound=min_n_obs,
            upper_bound=max_n_obs,
            shared_counts=shared_counts,
        )

    def _transform(
        self,
        adata: AnnData,
        modalities: Set,
        vars_to_keep: Optional[Dict[str, List[int]]] = None,
    ) -> Optional[AnnData]:
        """Filter out genes in modalities.

        Arguments
        ---------
        adata
            Annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells,
            columns to genes.
        modalities
            Set of modalities to filter.
        vars_to_keep
            Dictionary of variables to keep for a modality. If `None`, no variables are
            excluded from filtering.

        Returns
        -------
        Optional[AnnData]
            Returns annotated data if `self.inplace` is set to `False`, `None`
            otherwise.
        """

        if vars_to_keep is None:
            vars_to_keep = {}

        if self.shared_counts:
            shared_counts = self._sum_shared_counts(adata=adata, modalities=modalities)

            filter_mask = self.get_filter_mask(
                adata=adata,
                modality=shared_counts > 0,
                vars_to_keep=list(set().union(*vars_to_keep.values())),
            )

            # Filtering by shared counts requires the count matrices of the modalities
            # to be of the same size. As the filter mask is the same across all
            # modalities, it suffices to filter one modality, which automatically
            # filters all.
            self._subset_var(
                adata=adata, modality=modalities.pop(), filter_mask=filter_mask
            )

            _log_progress(
                kind="cell_count",
                filter_mask=filter_mask,
                modality="shared",
                lower_bound=self.lower_bound,
                upper_bound=self.upper_bound,
            )
        else:
            for modality in modalities:
                filter_mask = self.get_filter_mask(
                    adata=adata,
                    modality=get_modality(adata=adata, modality=modality) > 0,
                    vars_to_keep=vars_to_keep.get(modality, []),
                )
                self._subset_var(
                    adata=adata, modality=modality, filter_mask=filter_mask
                )

                _log_progress(
                    kind="cell_count",
                    filter_mask=filter_mask,
                    modality=modality,
                    lower_bound=self.lower_bound,
                    upper_bound=self.upper_bound,
                )

        return adata if not self.inplace else None


class GeneCountFilter(CountFilterBase):
    def __init__(
        self,
        inplace: bool = True,
        enforce: bool = False,
        min_counts: Optional[int] = None,
        max_counts: Optional[int] = None,
        shared_counts: bool = False,
    ):
        """Class for filtering variables based on counts across all observations.

        Arguments
        ---------
        inplace
            Boolean flag to indicate whether operations on annotaded data should be
            performed inplace or not. Otherwise, it is ignored.
        enforce
            Boolean flag to enforce filtering modalities.
        min_counts
            Minimum number of counts of a variable to not be filtered out.
        max_counts
            Maximum number of counts of a variable to not be filtered out.
        shared_counts
            Boolearn flag to filter by shared counts.
        """

        super().__init__(
            inplace=inplace,
            enforce=enforce,
            lower_bound=min_counts,
            upper_bound=max_counts,
            shared_counts=shared_counts,
        )

    def _transform(
        self,
        adata: AnnData,
        modalities: Set,
        vars_to_keep: Optional[Dict[str, List[int]]] = None,
    ) -> Optional[AnnData]:
        """Filter out genes in modalities.

        Arguments
        ---------
        adata
            Annotated data matrix of shape `n_obs` × `n_vars`.
        modalities
            Set of modalities to filter.
        vars_to_keep
            Dictionary of variables to keep for a modality. If `None`, no variables are
            excluded from filtering.

        Returns
        -------
        Optional[AnnData]
            Returns annotated data if `self.inplace` is set to `False`, `None`
            otherwise.
        """

        if vars_to_keep is None:
            vars_to_keep = {}

        if self.shared_counts:
            shared_counts = self._sum_shared_counts(adata=adata, modalities=modalities)

            filter_mask = self.get_filter_mask(
                adata=adata,
                modality=shared_counts,
                vars_to_keep=list(set().union(*vars_to_keep.values())),
            )

            # Filtering by shared counts requires the count matrices of the modalities
            # to be of the same size. As the filter mask is the same across all
            # modalities, it suffices to filter one modality, which automatically
            # filters all.
            self._subset_var(
                adata=adata, modality=modalities.pop(), filter_mask=filter_mask
            )

            _log_progress(
                kind="gene_count",
                filter_mask=filter_mask,
                modality="shared",
                lower_bound=self.lower_bound,
                upper_bound=self.upper_bound,
            )
        else:
            for modality in modalities:
                filter_mask = self.get_filter_mask(
                    adata=adata,
                    modality=get_modality(adata=adata, modality=modality),
                    vars_to_keep=vars_to_keep.get(modality, []),
                )
                self._subset_var(
                    adata=adata, modality=modality, filter_mask=filter_mask
                )

                _log_progress(
                    kind="gene_count",
                    filter_mask=filter_mask,
                    modality=modality,
                    lower_bound=self.lower_bound,
                    upper_bound=self.upper_bound,
                )

        return adata if not self.inplace else None
