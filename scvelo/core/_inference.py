from abc import ABCMeta, abstractmethod

from anndata import AnnData


class BaseInference(ABCMeta):
    """Base Inference class for all velocity methods."""

    def __init__(self, adata: AnnData):
        self.adata = adata
        super().__init__()

    @abstractmethod
    def fit(self, *args, **kwargs):
        """Fit the model."""

    @abstractmethod
    def state_dict(self):
        """Return the state of the model."""

    @abstractmethod
    def load_state_dict(self, state_dict):
        """Load the state of the model."""

    def export_results_adata(self):
        """Export the results to the AnnData object."""
