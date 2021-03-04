from functools import partial

import matplotlib.pyplot as pl

# todo: auto-complete and docs wrapper
from .scatter import scatter
from .utils import get_figure_params, hist
from .velocity_embedding import velocity_embedding
from .velocity_embedding_grid import velocity_embedding_grid
from .velocity_embedding_stream import velocity_embedding_stream
from .velocity_graph import velocity_graph


def _wraps_plot(wrapper, func):
    args = {"self", "kwargs"}
    annots_orig = {k: v for k, v in wrapper.__annotations__.items() if k not in args}
    annots = {k: v for k, v in func.__annotations__.items()}
    wrapper.__annotations__ = {**annots, **annots_orig}
    wrapper.__wrapped__ = func
    return wrapper


_wraps_plot_scatter = partial(_wraps_plot, func=scatter)
_wraps_plot_hist = partial(_wraps_plot, func=hist)
_wraps_plot_velocity_graph = partial(_wraps_plot, func=velocity_graph)
_wraps_plot_velocity_embedding = partial(_wraps_plot, func=velocity_embedding)
_wraps_plot_velocity_embedding_grid = partial(_wraps_plot, func=velocity_embedding_grid)
_wraps_plot_velocity_embedding_stream = partial(
    _wraps_plot, func=velocity_embedding_stream
)


def gridspec(ncols=4, nrows=1, figsize=None, dpi=None):
    figsize, dpi = get_figure_params(figsize, dpi, ncols)
    gs = pl.GridSpec(
        nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi)
    )
    return gs


class GridSpec:
    def __init__(self, ncols=4, nrows=1, figsize=None, dpi=None, **scatter_kwargs):
        """Specifies the geometry of the grid that a subplots can be placed in

        Example

        .. code:: python

            with scv.GridSpec() as pl:
                pl.scatter(adata, basis='pca')
                pl.scatter(adata, basis='umap')
                pl.hist(adata.obs.initial_size)

        Parameters
        ----------
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
        """
        self.ncols, self.nrows, self.figsize, self.dpi = ncols, nrows, figsize, dpi
        self.scatter_kwargs = scatter_kwargs
        self.scatter_kwargs.update({"show": False})
        self.get_new_grid()
        self.new_row = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.new_row and self.count < self.max_count:
            ax = pl.subplot(self.gs[self.max_count - 1])
            ax.axis("off")
        pl.show()

    def get_new_grid(self):
        self.gs = gridspec(self.ncols, self.nrows, self.figsize, self.dpi)
        geo = self.gs[0].get_geometry()
        self.max_count, self.count, self.new_row = geo[0] * geo[1], 0, True

    def get_ax(self):
        if self.count >= self.max_count:
            self.get_new_grid()
        self.count += 1
        return pl.subplot(self.gs[self.count - 1])

    def get_kwargs(self, kwargs=None):
        _kwargs = self.scatter_kwargs.copy()
        if kwargs is not None:
            _kwargs.update(kwargs)
        _kwargs.update({"ax": self.get_ax(), "show": False})
        return _kwargs

    @_wraps_plot_scatter
    def scatter(self, adata, **kwargs):
        return scatter(adata, **self.get_kwargs(kwargs))

    @_wraps_plot_velocity_embedding
    def velocity_embedding(self, adata, **kwargs):
        return velocity_embedding(adata, **self.get_kwargs(kwargs))

    @_wraps_plot_velocity_embedding_grid
    def velocity_embedding_grid(self, adata, **kwargs):
        return velocity_embedding_grid(adata, **self.get_kwargs(kwargs))

    @_wraps_plot_velocity_embedding_stream
    def velocity_embedding_stream(self, adata, **kwargs):
        return velocity_embedding_stream(adata, **self.get_kwargs(kwargs))

    @_wraps_plot_velocity_graph
    def velocity_graph(self, adata, **kwargs):
        return velocity_graph(adata, **self.get_kwargs(kwargs))

    @_wraps_plot_hist
    def hist(self, array, **kwargs):
        return hist(array, **self.get_kwargs(kwargs))
