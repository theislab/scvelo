"""Settings
"""

# set global verbosity level to show errors(0), warnings(1), info(2) and hints(3)
verbosity = 3


plot_suffix = ''
"""Global suffix that is appended to figure filenames.
"""

file_format_data = 'h5ad'
"""File format for saving AnnData objects.
Allowed are 'txt', 'csv' (comma separated value file) for exporting and 'h5ad'
(hdf5) for lossless saving.
"""

file_format_figs = 'pdf'
"""File format for saving figures.
For example 'png', 'pdf' or 'svg'. Many other formats work as well (see
`matplotlib.pyplot.savefig`).
"""

autosave = False
"""Save plots/figures as files in directory 'figs'.
Do not show plots/figures interactively.
"""

autoshow = True
"""Show all plots/figures automatically if autosave == False.
There is no need to call the matplotlib pl.show() in this case.
"""

writedir = './write/'
"""Directory where the function scanpy.write writes to by default.
"""

cachedir = './cache/'
"""Default cache directory.
"""

figdir = './figures/'
"""Directory where plots are saved.
"""

max_memory = 15
"""Maximal memory usage in Gigabyte.
Is currently not well respected....
"""

n_jobs = 1
"""Default number of jobs/ CPUs to use for parallel computing.
"""

logfile = ''
"""Name of logfile. By default is set to '' and writes to standard output."""

categories_to_ignore = ['N/A', 'dontknow', 'no_gate', '?']
"""Categories that are omitted in plotting etc.
"""

_frameon = False
"""See set_figure_params.
"""

_rcParams_style = None
"""See set_figure_params.
"""


# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------

from matplotlib import rcParams
from scanpy.plotting.rcmod import set_rcParams_scanpy
from scanpy.plotting.utils import default_palette


def set_rcParams_scvelo(fontsize=8, color_map=None, frameon=None):
    """Set matplotlib.rcParams to scvelo defaults."""

    # dpi options (mpl default: 100, 100)
    rcParams['figure.dpi'] = 100
    rcParams['savefig.dpi'] = 150

    # figure (mpl default: 0.125, 0.96, 0.15, 0.91)
    rcParams['figure.figsize'] = (7, 5)
    rcParams['figure.subplot.left'] = 0.18
    rcParams['figure.subplot.right'] = 0.96
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.top'] = 0.91

    # lines (defaults:  1.5, 6, 1)
    rcParams['lines.linewidth'] = 1.5  # the line width of the frame
    rcParams['lines.markersize'] = 6
    rcParams['lines.markeredgewidth'] = 1

    # font
    rcParams['font.sans-serif'] = \
        ['Arial', 'Helvetica', 'DejaVu Sans',
         'Bitstream Vera Sans', 'sans-serif']

    fontsize = fontsize
    labelsize = 0.9 * fontsize

    # fonsizes (mpl default: 10, medium, large, medium)
    rcParams['font.size'] = fontsize
    rcParams['legend.fontsize'] = labelsize
    rcParams['axes.titlesize'] = fontsize
    rcParams['axes.labelsize'] = labelsize

    # legend (mpl default: 1, 1, 2, 0.8)
    rcParams['legend.numpoints'] = 1
    rcParams['legend.scatterpoints'] = 1
    rcParams['legend.handlelength'] = 0.5
    rcParams['legend.handletextpad'] = 0.4

    # color cycle
    rcParams['axes.prop_cycle'] = default_palette()

    # axes
    rcParams['axes.linewidth'] = 0.8
    rcParams['axes.edgecolor'] = 'black'
    rcParams['axes.facecolor'] = 'white'

    # ticks (mpl default: k, k, medium, medium)
    rcParams['xtick.color'] = 'k'
    rcParams['ytick.color'] = 'k'
    rcParams['xtick.labelsize'] = labelsize
    rcParams['ytick.labelsize'] = labelsize

    # axes grid (mpl default: False, #b0b0b0)
    rcParams['axes.grid'] = False
    rcParams['grid.color'] = '.8'

    # color map
    rcParams['image.cmap'] = 'RdBu_r' if color_map is None else color_map

    # frame (mpl default: True)
    frameon = False if frameon is None else frameon
    global _frameon
    _frameon = frameon


def set_figure_params(style='scvelo', figsize=None, dpi=None, dpi_save=None, frameon=None, vector_friendly=True,
                      color_map=None, format='pdf', transparent=False, ipython_format='png2x'):
    """Set resolution/size, styling and format of figures.

    Arguments
    ---------
    style : `str` (default: `None`)
        Init default values for ``matplotlib.rcParams`` suited for `scvelo` or `scanpy`.
        Use `None` for the default matplotlib values.
    figsize: `[float, float]` (default: `None`)
        Width and height for default figure size.
    dpi : `int` (default: `None`)
        Resolution of rendered figures - this influences the size of figures in notebooks.
    dpi_save : `int` (default: `None`)
        Resolution of saved figures. This should typically be higher to achieve
        publication quality.
    frameon : `bool` (default: `None`)
        Add frames and axes labels to scatter plots.
    vector_friendly : `bool` (default: `True`)
        Plot scatter plots using `png` backend even when exporting as `pdf` or `svg`.
    color_map : `str` (default: `None`)
        Convenience method for setting the default color map.
    format : {'png', 'pdf', 'svg', etc.} (default: 'pdf')
        This sets the default format for saving figures: `file_format_figs`.
    transparent : `bool` (default: `True`)
        Save figures with transparent back ground. Sets
        `rcParams['savefig.transparent']`.
    ipython_format : list of `str` (default: 'png2x')
        Only concerns the notebook/IPython environment; see
        `IPython.core.display.set_matplotlib_formats` for more details.
    """
    try:
        import IPython
        IPython.core.display.set_matplotlib_formats(ipython_format)
    except:
        pass
    from matplotlib import rcParams
    global _rcParams_style
    _rcParams_style = style
    global _vector_friendly
    _vector_friendly = vector_friendly
    global file_format_figs
    file_format_figs = format
    if transparent is not None:
        rcParams['savefig.transparent'] = transparent
    if style is 'scvelo':
        set_rcParams_scvelo(color_map=color_map, frameon=frameon)
    elif style is 'scanpy':
        # dpi is not specified by scanpy directly in the defaults
        if dpi is None:
            rcParams['figure.dpi'] = 80
        if dpi_save is None:
            rcParams['savefig.dpi'] = 150
        frameon = True if frameon is None else frameon
        global _frameon
        _frameon = frameon
        set_rcParams_scanpy(color_map=color_map)
    # Overwrite style options if given
    if figsize is not None:
        rcParams['figure.figsize'] = figsize
    if dpi is not None:
        rcParams['figure.dpi'] = dpi
    if dpi_save is not None:
        rcParams['savefig.dpi'] = dpi_save


def set_rcParams_defaults():
    """Reset `matplotlib.rcParams` to defaults."""
    from matplotlib import rcParamsDefault
    rcParams.update(rcParamsDefault)


# ------------------------------------------------------------------------------
# Private global variables & functions
# ------------------------------------------------------------------------------

_vector_friendly = False
"""Set to true if you want to include pngs in svgs and pdfs.
"""

_low_resolution_warning = True
"""Print warning when saving a figure with low resolution."""

def _set_start_time():
    from time import time
    return time()

_start = _set_start_time()
"""Time when the settings module is first imported."""

_previous_time = _start
"""Variable for timing program parts."""

_previous_memory_usage = -1
"""Stores the previous memory usage."""
