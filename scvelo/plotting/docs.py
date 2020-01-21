"""Shared docstrings for plotting function parameters.
"""
from textwrap import dedent


def doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """
    def dec(obj):
        obj.__doc__ = dedent(obj.__doc__).format(**kwds)
        return obj
    return dec


doc_scatter = """\
basis: `str` (default='umap')
    Key for embedding.
color: `str`,  list of `str` or `None` (default: `None`)
    Key for annotations of observations/cells or variables/genes
use_raw : `bool` (default: `None`)
    Use `raw` attribute of `adata` if present.
layer: `str`,  list of `str` or `None` (default: `None`)
    Specify the layer for `color`.
color_map: `str` (default: `matplotlib.rcParams['image.cmap']`)
    String denoting matplotlib color map.
colorbar: `bool` (default: `False`)
    Whether to show colorbar.
palette: list of `str` (default: `None`)
    Colors to use for plotting groups (categorical annotation).
size: `float` (default: 5)
    Point size.
alpha: `float` (default: 1)
    Set blending - 0 transparent to 1 opaque.
linewidth: `float` (default: 1)
    Scaling factor for the width of occurring lines.
linecolor: `str` ir list of `str` (default: 'k')
    Color of lines from velocity fits, linear fits and polynomial fits
perc: tuple, e.g. [2,98] (default: `None`)
    Specify percentile for continuous coloring.
groups: `str` (default: `all groups`)
    Restrict to a few categories in categorical observation annotation.
sort_order: `bool` (default: `True`)
    For continuous annotations used as color parameter, plot data points with higher values on top of others.
components: `str` or list of `str` (default: '1,2')
    For instance, ['1,2', '2,3'].
projection: {'2d', '3d'} (default: '2d')
    Projection of plot.
legend_loc: str (default: 'none')
    Location of legend, either 'on data', 'right margin' or valid keywords
    for matplotlib.legend.
legend_fontsize: `int` (default: `None`)
    Legend font size.
legend_fontweight: {'normal', 'bold', ...} (default: `None`)
    Legend font weight. Defaults to 'bold' if `legend_loc = 'on data'`, otherwise to 'normal'. 
    Available are `['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']`.
right_margin: `float` or list of `float` (default: `None`)
    Adjust the width of the space right of each plotting panel.
left_margin: `float` or list of `float` (default: `None`)
    Adjust the width of the space left of each plotting panel.
xlabel: `str` (default: `None`)
    Label of x-axis.
ylabel: `str` (default: `None`)
    Label of y-axis.
title: `str` (default: `None`)
    Provide title for panels either as, e.g. `["title1", "title2", ...]`.
fontsize: `float` (default: `None`)
    Label font size.
figsize: tuple (default: `(7,5)`)
    Figure size.
xlim: tuple, e.g. [0,1] or `None` (default: `None`)
    Restrict x-limits of the axis.
ylim: tuple, e.g. [0,1] or `None` (default: `None`)
    Restrict y-limits of the axis.
show_density: `bool` or `str` or `None` (default: `None`)
    Whether to show density of counts attached to the x and y axes. Color of the plot can also be given as `str`.
show_assignments: `bool` or `str` or `None` (default: `None`)
    Whether to show assignments to the model curve. Color of the assignments can also be given as `str`.
show_linear_fit: `bool` or `str` or `None` (default: `None`)
    Whether to show linear fit to the data. Color of the plot can also be given as `str`.
show_polyfit: `bool` or `str` or `int` or `None` (default: `None`)
    Whether to show polynomial fit to ???. Color of the polyfit plot can also be given as `str`. If `int` is given,
    determines the degree of the polynomial fit (else the degree defaults to 2).
rug: `str` or `None` (default: `None`)
    If categorical observation annotation (e.g. 'clusters') is given, a rugplot is attached to the x-axis showing the 
    distribution of data membership to each of the categories.
add_outline: `bool` (default: `False`)
    Whether to show an outline around scatter plot dots
outline_width: tuple type `scalar` or `None` (default: `(0.3, 0.05)`)
    Width of the inner and outer outline
outline_color: tuple of type `str` or `None` (default: `('black', 'white')`)
    Inner and outer matplotlib color of the outline
n_convolve: `int` or `None` (default: `None`)
    If `int` is given, data is smoothed by convolution along the x-axis with kernel size n_convolve.
smooth: `bool` or `int` (default: `None`)
    Whether to convolve/average the color values over the nearest neighbors. If `int`, it specifies number of neighbors.
dpi: `int` (default: 80)
    Figure dpi.
frameon: `bool` (default: `True`)
    Draw a frame around the scatter plot.
ncols : `int` (default: `None`)
    Number of panels per row.
wspace : `float` (default: None)
    Adjust the width of the space between multiple panels.
hspace : `float` (default: None)
    Adjust the height of the space between multiple panels.
show: `bool`, optional (default: `None`)
    Show the plot, do not return axis.
save: `bool` or `str`, optional (default: `None`)
    If `True` or a `str`, save the figure. A string is appended to the default filename. 
    Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
ax: `matplotlib.Axes`, optional (default: `None`)
    A matplotlib axes object. Only works if plotting a single component.\
"""