from .scatter import scatter
from .velocity_embedding import velocity_embedding
from .velocity_embedding_grid import velocity_embedding_grid
from .velocity_embedding_stream import velocity_embedding_stream
from .velocity_graph import velocity_graph
from .docs import doc_scatter, doc_params

import matplotlib.pyplot as pl
from matplotlib import rcParams


def _wraps_plot_scatter(wrapper):
    annots_orig = {k: v for k, v in wrapper.__annotations__.items() if k not in {'self', 'adata', 'kwargs'}}
    annots = {k: v for k, v in scatter.__annotations__.items()}
    wrapper.__annotations__ = {**annots, **annots_orig}
    wrapper.__wrapped__ = scatter
    return wrapper


def _wraps_plot_velocity_embedding(wrapper):
    annots_orig = {k: v for k, v in wrapper.__annotations__.items() if k not in {'self', 'adata', 'kwargs'}}
    annots = {k: v for k, v in velocity_embedding.__annotations__.items()}
    wrapper.__annotations__ = {**annots, **annots_orig}
    wrapper.__wrapped__ = velocity_embedding
    return wrapper


def gridspec(ncols=4, nrows=1, figsize=None, dpi=None):
    if figsize is None: figsize = rcParams['figure.figsize']
    gs = pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))
    return gs


class GridSpec:
    def __init__(self, ncols=4, nrows=1, figsize=None, dpi=None, **scatter_kwargs):
        """\
        GridSpec

        Arguments
        ---------
        ncols: `int` (default: 4)
            Number of panels per row.
        nrows: `int` (default: 1)
            Number of panels per column.
        figsize: tuple (default: `None`)
            Figure size.
        dpi: `int` (default: `None`)
            Figure dpi.
        scatter_kwargs:
            Arguments to be passed to all scatter plots, e.g. `frameon=False`.

        Returns
        -------
            If `show==False` a `matplotlib.Axis`
        """
        self.ncols, self.nrows, self.figsize, self.dpi = ncols, nrows, figsize, dpi
        self.scatter_kwargs = scatter_kwargs
        self.scatter_kwargs.update({'show': False})
        self.initialize_grid()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pl.show()

    def initialize_grid(self):
        self.gs = gridspec(self.ncols, self.nrows, self.figsize, self.dpi)
        geo = self.gs[0].get_geometry()
        self.max_count, self.count = geo[0] * geo[1], 0

    def get_ax(self):
        self.count += 1
        if self.count > self.max_count:
            self.initialize_grid()
        return pl.subplot(self.gs[self.count - 1])

    @_wraps_plot_scatter
    @doc_params(scatter=doc_scatter)
    def scatter(self, adata, **kwargs):
        """\
        Scatter plot along observations or variables axes.

        Parameters
        ---------
        adata: :class:`~anndata.AnnData`
            Annotated data matrix.
        {scatter}

        Returns
        -------
        If `show==False` a `matplotlib.Axis`
        """
        return scatter(adata, ax=self.get_ax(), **kwargs, **self.scatter_kwargs)

    @_wraps_plot_velocity_embedding
    def velocity_embedding(self, adata, **kwargs):
        return velocity_embedding(adata, ax=self.get_ax(), **kwargs, **self.scatter_kwargs)
