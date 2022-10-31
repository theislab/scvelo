from abc import ABCMeta, abstractmethod

from anndata import AnnData


class BaseInference(ABCMeta):
    """Base Inference class for all velocity methods."""

    def __init__(self, adata: AnnData):
        self._adata = adata
        self._state_dict = None
        super().__init__()

    @abstractmethod
    def fit(self, *args, **kwargs):
        """Fit the model."""

    def state_dict(self):
        """Return the state of the model."""
        return self._state_dict

    @abstractmethod
    def load_state_dict(self, state_dict):
        """Load the state of the model."""

    @abstractmethod
    def export_results_adata(self):
        """Export the results to the AnnData object."""

    def get_velocity(self, *args, **kwargs):
        """Return the velocity."""
