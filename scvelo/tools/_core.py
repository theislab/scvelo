from abc import abstractmethod
from typing import NamedTuple

import torch

from anndata import AnnData


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    U_KEY: str = "U"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()

DEFAULT_ACTIVATION_FUNCTION = torch.nn.Softplus()


class BaseInference:
    """Base Inference class for all velocity methods."""

    def __init__(self, adata: AnnData):
        self._adata = adata
        self._state_dict = None
        super().__init__()

    @abstractmethod
    def fit(self, *args, **kwargs):
        """Fit the model."""

    @abstractmethod
    def state_dict(self):
        """Return the state of the model."""
        return self._state_dict

    @abstractmethod
    def export_results_adata(self):
        """Export the results to the AnnData object."""

    def get_velocity(self, *args, **kwargs):
        """Return the velocity."""
