import os
from collections import abc

from cycler import Cycler, cycler

import numpy as np
import pandas as pd
from pandas import Index
from scipy import stats
from scipy.sparse import issparse

import matplotlib.pyplot as pl
import matplotlib.transforms as tx
from matplotlib import patheffects, rcParams
from matplotlib.collections import LineCollection
from matplotlib.colors import cnames, is_color_like, ListedColormap, to_rgb
from matplotlib.gridspec import SubplotSpec
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scvelo import logging as logg
from scvelo import settings
from scvelo.tools.utils import strings_to_categoricals
from . import palettes

try:
    from numpy.core._exceptions import UFuncTypeError
except ModuleNotFoundError:
    from numpy._core._exceptions import UFuncTypeError


"""helper functions"""


# TODO: Add docstrings
def make_dense(X):
    """TODO."""
    if issparse(X):
        XA = X.toarray() if X.ndim == 2 else X.A1
    else:
        XA = X.A1 if isinstance(X, np.matrix) else X
    return np.array(XA)


# TODO: Add docstrings
def is_view(adata):
    """TODO."""
    return (
        adata.is_view
        if hasattr(adata, "is_view")
        else adata.isview
        if hasattr(adata, "isview")
        else adata._isview
        if hasattr(adata, "_isview")
        else True
    )


# TODO: Add docstrings
def is_categorical(data, c=None):
    """TODO."""
    if c is None:
        return isinstance(
            data.dtype, pd.CategoricalDtype
        )  # if data is categorical/array
    if not is_view(data):  # if data is anndata view
        strings_to_categoricals(data)
    return (
        isinstance(c, str)
        and c in data.obs.keys()
        and isinstance(data.obs[c].dtype, pd.CategoricalDtype)
    )


# TODO: Add docstrings
def is_int(key):
    """TODO."""
    return isinstance(key, (int, np.integer))


# TODO: Add docstrings
def is_list(key):
    """TODO."""
    return isinstance(key, (list, tuple, np.record))


# TODO: Add docstrings
def is_list_or_array(key):
    """TODO."""
    return isinstance(key, (list, tuple, np.record, np.ndarray))


# TODO: Add docstrings
def is_list_of_str(key, max_len=None):
    """TODO."""
    if max_len is not None:
        return (
            is_list_or_array(key)
            and len(key) < max_len
            and all(isinstance(item, str) for item in key)
        )
    else:
        return is_list(key) and all(isinstance(item, str) for item in key)


# TODO: Add docstrings
def is_list_of_list(lst):
    """TODO."""
    return lst is not None and any(
        isinstance(list_element, list) for list_element in lst
    )


# TODO: Add docstrings
def is_list_of_int(lst):
    """TODO."""
    return is_list_or_array(lst) and all(is_int(item) for item in lst)


# TODO: Add docstrings
def to_list(key, max_len=20):
    """TODO."""
    if isinstance(key, Index) or is_list_of_str(key, max_len):
        key = list(key)
    return key if is_list(key) and (max_len is None or len(key) < max_len) else [key]


# TODO: Add docstrings
def to_val(key):
    """TODO."""
    return key[0] if isinstance(key, (list, tuple)) and len(key) == 1 else key


# TODO: Add docstrings
def get_figure_params(figsize, dpi=None, ncols=1):
    """TODO."""
    figsize = rcParams["figure.figsize"] if figsize is None else figsize
    dpi = rcParams["figure.dpi"] if dpi is None else dpi
    if settings.presenter_view and figsize[0] * ncols * (dpi / 80) > 12:
        figscale = 12 / (figsize[0] * ncols)
        figsize = (figsize[0] * figscale, figsize[1] * figscale)
        dpi = min(80, dpi)
    return figsize, dpi


# TODO: Add docstrings
def get_ax(ax, show=None, figsize=None, dpi=None, projection=None):
    """TODO."""
    figsize, _ = get_figure_params(figsize)
    if ax is None:
        projection = "3d" if projection == "3d" else None
        _, ax = pl.subplots(
            figsize=figsize, dpi=dpi, subplot_kw={"projection": projection}
        )
    elif isinstance(ax, SubplotSpec):
        geo = ax.get_geometry()
        if show is None:
            show = geo[-1] + 1 == geo[0] * geo[1]
        ax = pl.subplot(ax)
    return ax, show


# TODO: Add docstrings
def get_kwargs(kwargs, dict_new_kwargs):
    """TODO."""
    kwargs = kwargs.copy()
    kwargs.update(dict_new_kwargs)
    return kwargs


# TODO: Add docstrings
def check_basis(adata, basis):
    """TODO."""
    if basis in adata.obsm.keys() and f"X_{basis}" not in adata.obsm.keys():
        adata.obsm[f"X_{basis}"] = adata.obsm[basis]
        logg.info(f"Renamed '{basis}' to convention 'X_{basis}' (adata.obsm).")


# TODO: Add docstrings
def get_basis(adata, basis):
    """TODO."""
    if isinstance(basis, str) and basis.startswith("X_"):
        basis = basis[2:]
    check_basis(adata, basis)
    return basis


# TODO: Add docstrings
def to_valid_bases_list(adata, keys):
    """TODO."""
    if isinstance(keys, pd.DataFrame):
        keys = keys.index
    if not isinstance(keys, str):
        keys = list(np.ravel(keys))
    keys = to_list(keys, max_len=np.inf)
    if all(isinstance(item, str) for item in keys):
        for i, key in enumerate(keys):
            if key.startswith("X_"):
                keys[i] = key = key[2:]
            check_basis(adata, key)
        valid_keys = np.hstack(
            [
                adata.obs.keys(),
                adata.var.keys(),
                adata.varm.keys(),
                adata.obsm.keys(),
                [key[2:] for key in adata.obsm.keys()],
                list(adata.layers.keys()),
            ]
        )
        keys_ = keys
        keys = [key for key in keys if key in valid_keys or key in adata.var_names]
        keys_ = [key for key in keys_ if key not in keys]
        if len(keys_) > 0:
            msg_embedding = ""
            if len(keys_) == 1 and keys_[0] in {"diffmap", "umap", "tsne"}:
                msg_embedding = f"You need to run `scv.tl.{keys_[0]}` first."
            logg.warn(", ".join(keys_), "not found.", msg_embedding)
    return keys


# TODO: Add docstrings
def get_components(components=None, basis=None, projection=None):
    """TODO."""
    if components is None:
        components = "1,2,3" if projection == "3d" else "1,2"
    if isinstance(components, str):
        components = components.split(",")
    components = np.array(components).astype(int) - 1
    if "diffmap" in basis or "vmap" in basis:
        components += 1
    return components


# TODO: Add docstrings
def get_obs_vector(adata, basis, layer=None, use_raw=None):
    """TODO."""
    return (
        adata.obs_vector(basis, layer=layer)
        if layer in adata.layers.keys()
        else adata.raw.obs_vector(basis)
        if use_raw
        else adata.obs_vector(basis)
    )


# TODO: Add docstrings
def get_value_counts(adata, color):
    """TODO."""
    value_counts = adata.obs[color].value_counts()
    probs = np.array(adata.obs[color])
    for cat in value_counts.index:
        probs[probs == cat] = value_counts[cat]
    return np.array(probs, dtype=np.float32)


# TODO: Add docstrings
def get_groups(adata, groups, groupby=None):
    """TODO."""
    if not isinstance(groupby, str) or groupby not in adata.obs.keys():
        groupby = (
            "clusters"
            if "clusters" in adata.obs.keys()
            else "louvain"
            if "louvain" in adata.obs.keys()
            else None
        )
    if groups is True:
        return None, groupby
    if groups is not None and not isinstance(groups, str) and len(groups) == 1:
        groups = groups[0]
    if isinstance(groups, str):
        cats = [""]
        if is_categorical(adata, groupby):
            cats = adata.obs[groupby].cat.categories
        if ":" in groups and not np.any([":" in cat for cat in cats]):
            groupby, groups = groups.split(":")
            groups = groups.strip()
        if "," in groups and not np.any(["," in cat for cat in cats]):
            groups = [a.strip() for a in groups.split(",")]
    if isinstance(groups, str):
        groups = [groups]
    return groups, groupby


# TODO: Add docstrings
def groups_to_bool(adata, groups, groupby=None):
    """TODO."""
    groups, groupby = get_groups(adata, groups, groupby)
    if isinstance(groups, (list, tuple, np.ndarray, np.record)):
        if groupby is not None and isinstance(groups[0], str):
            groups = np.array([key in groups for key in adata.obs[groupby]])

    if groupby is not None and groupby in adata.obs.keys():
        c = adata.obs[groupby]
        if np.any(pd.isnull(c)):
            valid = np.array(~pd.isnull(c))
            groups = (
                valid if groups is None or len(groups) != len(c) else groups & valid
            )

    groups = (
        np.ravel(groups)
        if isinstance(groups, (list, tuple, np.ndarray, np.record))
        else None
    )
    return groups


# TODO: Add docstrings
def gets_vals_from_color_gradients(adata, color=None, **scatter_kwargs):
    """TODO."""
    color_gradients = scatter_kwargs.pop("color_gradients")
    scatter_kwargs.update({"color_gradients": None})
    if "colorbar" not in scatter_kwargs or scatter_kwargs["colorbar"] is None:
        scatter_kwargs.update({"colorbar": False})
    if "s" not in scatter_kwargs:
        size = scatter_kwargs["size"]
        scatter_kwargs["s"] = default_size(adata) if size is None else size
    if not any([v in scatter_kwargs for v in ["vmin", "vmax", "vmid"]]):
        scatter_kwargs["vmid"] = 0
    if isinstance(color_gradients, str) and color_gradients in adata.obsm.keys():
        if color is None:
            color = color_gradients
        color_gradients = adata.obsm[color_gradients]
    elif (
        isinstance(color_gradients, (list, tuple))
        and color_gradients[0] in adata.obs.keys()
    ):
        color_gradients = pd.DataFrame(
            np.stack([adata.obs[c] for c in color_gradients]).T, columns=color_gradients
        )
    if color is None:
        color = "clusters_gradients"
    palette = scatter_kwargs.pop("palette")
    if palette is None and hasattr(color_gradients, "colors"):
        palette = list(color_gradients.colors)

    pd_colgrad = pd.DataFrame(color_gradients)
    vals = np.clip(pd_colgrad.values, 0, None)
    names = (
        color_gradients.names
        if hasattr(color_gradients, "names")
        else pd_colgrad.columns
    )

    adata.obs[color] = pd.Categorical(
        [f"{names[i]}" for i in np.argmax(vals, 1)], categories=names
    )
    set_colors_for_categorical_obs(adata, color, palette)

    return vals, names, color, scatter_kwargs


"""get default parameters"""


# TODO: Add docstrings
def default_basis(adata, **kwargs):
    """TODO."""
    if "x" in kwargs and "y" in kwargs:
        keys, x, y = ["embedding"], kwargs.pop("x"), kwargs.pop("y")
        adata.obsm["X_embedding"] = np.stack([x, y]).T
        if "velocity_embedding" in adata.obsm.keys():
            del adata.obsm["velocity_embedding"]
    else:
        keys = [
            key for key in ["pca", "tsne", "umap"] if f"X_{key}" in adata.obsm.keys()
        ]
    if not keys:
        raise ValueError("No basis specified.")
    return keys[-1] if len(keys) > 0 else None


# TODO: Add docstrings
def default_size(adata):
    """TODO."""
    adjusted, classic = 1.2e5 / adata.n_obs, 20
    return (
        np.mean([adjusted, classic])
        if settings._rcParams_style == "scvelo"
        else adjusted
    )


# TODO: Add docstrings
def default_color(adata, add_outline=None):
    """TODO."""
    if (
        isinstance(add_outline, str)
        and add_outline in adata.var.keys()
        and "recover_dynamics" in adata.uns.keys()
        and add_outline in adata.uns["recover_dynamics"]
    ):
        return adata.uns["recover_dynamics"][add_outline]
    return (
        "clusters"
        if "clusters" in adata.obs.keys()
        else "louvain"
        if "louvain" in adata.obs.keys()
        else "grey"
    )


# TODO: Add docstrings
def default_color_map(adata, c):
    """TODO."""
    cmap = None
    if isinstance(c, str) and c in adata.obs.keys() and not is_categorical(adata, c):
        c = adata.obs[c]
    elif isinstance(c, int):
        cmap = "viridis_r"
    if len(np.array(c).flatten()) == adata.n_obs:
        try:
            if np.min(c) in [-1, 0, False] and np.max(c) in [1, True]:
                cmap = "viridis_r"
        except UFuncTypeError as e:
            logg.warn(f"Setting `cmap` to `None`: {e}")
            cmap = None
    return cmap


# TODO: Add docstrings
def default_legend_loc(adata, color, legend_loc):
    """TODO."""
    n_categories = 0
    if is_categorical(adata, color):
        n_categories = len(adata.obs[color].cat.categories)
    if legend_loc is False:
        legend_loc = "none"
    elif legend_loc is None:
        legend_loc = "upper right" if n_categories <= 4 else "on data"
    return legend_loc


# TODO: Add docstrings
def default_xkey(adata, use_raw):
    """TODO."""
    use_raw = "spliced" in adata.layers.keys() and (
        use_raw or "Ms" not in adata.layers.keys()
    )
    return "spliced" if use_raw else "Ms" if "Ms" in adata.layers.keys() else "X"


# TODO: Add docstrings
def default_ykey(adata, use_raw):
    """TODO."""
    use_raw = "unspliced" in adata.layers.keys() and (
        use_raw or "Mu" not in adata.layers.keys()
    )
    return "unspliced" if use_raw else "Mu" if "Mu" in adata.layers.keys() else None


# TODO: Add docstrings
def default_arrow(size):
    """TODO."""
    if isinstance(size, (list, tuple)) and len(size) == 3:
        head_l, head_w, ax_l = size
    elif isinstance(size, (int, np.integer, float)):
        head_l, head_w, ax_l = 12 * size, 10 * size, 8 * size
    else:
        head_l, head_w, ax_l = 12, 10, 8
    return head_l, head_w, ax_l


"""set axes parameters (ticks, frame, labels, title, """


# TODO: Add docstrings
def update_axes(
    ax,
    xlim=None,
    ylim=None,
    fontsize=None,
    is_embedding=False,
    frameon=None,
    figsize=None,
    aspect="auto",
):
    """TODO."""
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    frameon = settings._frameon if frameon is None else frameon
    if isinstance(frameon, str) and frameon == "artist":
        set_artist_frame(ax, figsize=figsize)
    elif frameon:
        if is_embedding:
            kwargs = {
                "bottom": False,
                "left": False,
                "labelbottom": False,
                "labelleft": False,
            }
            ax.tick_params(which="both", **kwargs)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            labelsize = int(fontsize * 0.75) if fontsize is not None else None
            ax.tick_params(axis="both", which="major", labelsize=labelsize)
        if isinstance(frameon, str) and frameon != "full":
            frameon = "bl" if frameon == "half" else frameon
            bf, lf, tf, rf = (f in frameon for f in ["bottom", "left", "top", "right"])
            if not np.any([bf, lf, tf, rf]):
                bf, lf, tf, rf = (f in frameon for f in ["b", "l", "t", "r"])
            ax.spines["top"].set_visible(tf)
            ax.spines["right"].set_visible(rf)
            if not bf:
                ax.set_xlabel("")
                ax.spines["bottom"].set_visible(False)
            if not lf:
                ax.set_ylabel("")
                ax.spines["left"].set_visible(False)
            kwargs = {"bottom": bf, "left": lf, "labelbottom": bf, "labelleft": lf}
            ax.tick_params(which="both", top=tf, right=rf, **kwargs)

    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)
        kwargs = {
            "bottom": False,
            "left": False,
            "labelbottom": False,
            "labelleft": False,
        }
        ax.tick_params(which="both", **kwargs)
        ax.set_frame_on(False)

    ax.set_aspect(aspect)

    if rcParams["savefig.transparent"]:
        ax.patch.set_alpha(0)


# TODO: Add docstrings
def set_artist_frame(ax, length=0.2, figsize=None):
    """TODO."""
    ax.tick_params(
        axis="both",
        which="both",
        labelbottom=False,
        labelleft=False,
        bottom=False,
        left=False,
        top=False,
        right=False,
    )
    for side in ["bottom", "right", "top", "left"]:
        ax.spines[side].set_visible(False)
    figsize = rcParams["figure.figsize"] if figsize is None else figsize
    aspect_ratio = figsize[0] / figsize[1]
    ax.xaxis.set_label_coords(length * 0.45, -0.035)
    ax.yaxis.set_label_coords(-0.025, length * aspect_ratio * 0.45)
    ax.xaxis.label.set_size(ax.xaxis.label.get_size() / 1.2)
    ax.yaxis.label.set_size(ax.yaxis.label.get_size() / 1.2)

    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDirectionArrows

    kwargs = {
        "loc": 3,
        "pad": -1,
        "back_length": 0,
        "fontsize": 0,
        "aspect_ratio": aspect_ratio,
    }
    kwargs.update({"text_props": {"ec": "k", "fc": "k", "lw": 0.1}})
    kwargs.update({"arrow_props": {"ec": "none", "fc": "k"}, "length": length})
    arr = AnchoredDirectionArrows(ax.transAxes, "none", "none", **kwargs)
    ax.add_artist(arr)


# TODO: Add docstrings
def set_label(xlabel, ylabel, fontsize=None, basis=None, ax=None, **kwargs):
    """TODO."""
    labels = np.array(["Ms", "Mu", "X"])
    labels_new = np.array(["spliced", "unspliced", "expression"])
    if xlabel in labels:
        xlabel = labels_new[xlabel == labels][0]
    if ylabel in labels:
        ylabel = labels_new[ylabel == labels][0]
    if ax is None:
        ax = pl.gca()
    kwargs.update({"fontsize": fontsize})
    if basis is not None:
        component_name = (
            "DC"
            if "diffmap" in basis
            else "tSNE"
            if basis == "tsne"
            else "UMAP"
            if basis == "umap"
            else "PC"
            if basis == "pca"
            else basis.replace("draw_graph_", "").upper()
            if "draw_graph" in basis
            else basis
        )
        ax.set_xlabel(f"{component_name}1", **kwargs)
        ax.set_ylabel(f"{component_name}2", **kwargs)
    if isinstance(xlabel, str):
        ax.set_xlabel(xlabel.replace("_", " "), **kwargs)
    if isinstance(ylabel, str):
        rotation = 0 if ylabel.startswith("$") or len(ylabel) == 1 else 90
        ax.set_ylabel(ylabel.replace("_", " "), rotation=rotation, **kwargs)


# TODO: Add docstrings
def set_title(title, layer=None, color=None, fontsize=None, ax=None):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    color = color if isinstance(color, str) and not is_color_like(color) else None
    if isinstance(title, str):
        title = title.replace("_", " ")
    elif isinstance(layer, str) and isinstance(color, str):
        title = f"{color}  {layer}".replace("_", " ")
    elif isinstance(color, str):
        title = color.replace("_", " ")
    else:
        title = ""
    ax.set_title(title, fontsize=fontsize)


# TODO: Add docstrings
def set_frame(ax, frameon):
    """TODO."""
    frameon = settings._frameon if frameon is None else frameon
    if not frameon:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_frame_on(False)


# TODO: Add docstrings
def set_legend(
    adata,
    ax,
    value_to_plot,
    legend_loc,
    scatter_array,
    legend_fontweight,
    legend_fontsize,
    legend_fontoutline,
    legend_align_text,
    groups,
):
    """Adds a legend to the given ax with categorial data."""
    # add legend
    if legend_fontoutline is None:
        legend_fontoutline = 1
    obs_vals = adata.obs[value_to_plot]
    str_cats = obs_vals.cat.categories.astype(str)
    obs_vals = obs_vals.cat.set_categories(str_cats, rename=True)
    color_keys = adata.uns[f"{value_to_plot}_colors"]
    if isinstance(color_keys, dict):
        color_keys = np.array([color_keys[c] for c in obs_vals.cat.categories])
    valid_cats = np.where(obs_vals.value_counts()[obs_vals.cat.categories] > 0)[0]
    categories = np.array(obs_vals.cat.categories)[valid_cats]
    colors = np.array(color_keys)[valid_cats]

    if groups is not None:
        groups, groupby = get_groups(adata, groups, value_to_plot)
        # only label groups with the respective color
        groups = [g for g in groups if g in categories]
        colors = [colors[list(categories).index(x)] for x in groups]
        categories = groups

    if legend_loc == "on data":
        legend_fontweight = "bold" if legend_fontweight is None else legend_fontweight
        # identify centroids to put labels
        texts = []
        for label in categories:
            x_pos, y_pos = np.nanmedian(scatter_array[obs_vals == label, :], axis=0)
            if isinstance(label, str):
                label = label.replace("_", " ")
            kwargs = {"verticalalignment": "center", "horizontalalignment": "center"}
            kwargs.update({"weight": legend_fontweight, "fontsize": legend_fontsize})
            pe = [patheffects.withStroke(linewidth=legend_fontoutline, foreground="w")]
            text = ax.text(x_pos, y_pos, label, path_effects=pe, **kwargs)
            texts.append(text)

        if legend_align_text:
            autoalign = "y" if legend_align_text is True else legend_align_text
            try:
                from adjustText import adjust_text as adj_text

                adj_text(texts, autoalign=autoalign, text_from_points=False, ax=ax)
            except ImportError:
                print("Please `pip install adjustText` for auto-aligning texts")

    else:
        for idx, label in enumerate(categories):
            if isinstance(label, str):
                label = label.replace("_", " ")
            ax.scatter([], [], c=[colors[idx]], label=label)
        ncol = 1 if len(categories) <= 14 else 2 if len(categories) <= 30 else 3
        kwargs = {"frameon": False, "fontsize": legend_fontsize, "ncol": ncol}
        if legend_loc == "upper right":
            ax.legend(loc="upper left", bbox_to_anchor=(1, 1), **kwargs)
        elif legend_loc == "lower right":
            ax.legend(loc="lower left", bbox_to_anchor=(1, 0), **kwargs)
        elif "right" in legend_loc:  # 'right', 'center right', 'right margin'
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), **kwargs)
        elif legend_loc != "none":
            ax.legend(loc=legend_loc, **kwargs)


# TODO: Add docstrings
def set_margin(ax, x, y, add_margin):
    """TODO."""
    add_margin = 0.1 if add_margin is True else add_margin
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    xmargin = (xmax - xmin) * add_margin
    ymargin = (ymax - ymin) * add_margin
    ax.set_xlim(xmin - xmargin, xmax + xmargin)
    ax.set_ylim(ymin - ymargin, ymax + ymargin)


"""get color values"""


# TODO: Add docstrings
def clip(c, perc):
    """TODO."""
    if np.size(perc) < 2:
        perc = [perc, 100] if perc < 50 else [0, perc]
    lb, ub = np.percentile(c, perc)
    return np.clip(c, lb, ub)


# TODO: Add docstrings
def get_colors(adata, c):
    """TODO."""
    if is_color_like(c):
        return c
    else:
        if f"{c}_colors" not in adata.uns.keys():
            palette = default_palette(None)
            palette = adjust_palette(palette, length=len(adata.obs[c].cat.categories))
            n_cats = len(adata.obs[c].cat.categories)
            adata.uns[f"{c}_colors"] = palette[:n_cats].by_key()["color"]
        if isinstance(adata.uns[f"{c}_colors"], dict):
            cluster_ix = adata.obs[c].values
        else:
            cluster_ix = adata.obs[c].cat.codes.values
        return np.array(
            [
                adata.uns[f"{c}_colors"][cluster_ix[i]]
                if cluster_ix[i] != -1
                else "lightgrey"
                for i in range(adata.n_obs)
            ]
        )


# TODO: Add docstrings
def interpret_colorkey(adata, c=None, layer=None, perc=None, use_raw=None):
    """TODO."""
    if c is None:
        c = default_color(adata)
    if issparse(c):
        c = make_dense(c).flatten()
    if is_categorical(adata, c):
        c = get_colors(adata, c)
    elif isinstance(c, str):
        if is_color_like(c) and c not in adata.var_names:
            pass
        elif c in adata.obs.keys():  # color by observation key
            c = adata.obs[c]
        elif c in adata.var_names or (
            use_raw and adata.raw is not None and c in adata.raw.var_names
        ):  # by gene
            if layer in adata.layers.keys():
                if perc is None and any(
                    layer_name in layer
                    for layer_name in ["spliced", "unspliced", "Ms", "Mu", "velocity"]
                ):
                    perc = [1, 99]  # to ignore outliers in non-logarithmized layers
                c = adata.obs_vector(c, layer=layer)
            elif layer is not None and np.any(
                [
                    layer_name in layer or "X" in layer
                    for layer_name in adata.layers.keys()
                ]
            ):
                l_array = np.hstack(
                    [
                        adata.obs_vector(c, layer=layer)[:, None]
                        for layer in adata.layers.keys()
                    ]
                )
                l_array = pd.DataFrame(l_array, columns=adata.layers.keys())
                l_array.insert(0, "X", adata.obs_vector(c))
                c = np.array(l_array.astype(np.float32).eval(layer))
            else:
                if layer is not None and layer != "X":
                    logg.warn(layer, "not found. Using .X instead.")
                if adata.raw is None and use_raw:
                    raise ValueError("AnnData object does not have `raw` counts.")
                c = adata.raw.obs_vector(c) if use_raw else adata.obs_vector(c)
            c = c.toarray().flatten() if issparse(c) else c
        elif c in adata.var.keys():  # color by observation key
            c = adata.var[c]
        elif np.any([var_key in c for var_key in adata.var.keys()]):
            var_keys = [
                k for k in adata.var.keys() if not isinstance(adata.var[k][0], str)
            ]
            var = adata.var[list(var_keys)]
            c = var.astype(np.float32).eval(c)
        elif np.any([obs_key in c for obs_key in adata.obs.keys()]):
            obs_keys = [
                k for k in adata.obs.keys() if not isinstance(adata.obs[k][0], str)
            ]
            obs = adata.obs[list(obs_keys)]
            c = obs.astype(np.float32).eval(c)
        elif not is_color_like(c):
            raise ValueError(
                "color key is invalid! pass valid observation annotation or a gene name"
            )
        if not isinstance(c, str) and perc is not None:
            c = clip(c, perc=perc)
    else:
        c = np.array(c).flatten()
        if perc is not None:
            c = clip(c, perc=perc)
    return c


# adapted from scanpy
def set_colors_for_categorical_obs(adata, value_to_plot, palette=None):
    """Sets adata.uns[f'{value_to_plot}_colors'] to given palette or default colors.

    Parameters
    ----------
    adata
        annData object
    value_to_plot
        name of a valid categorical observation
    palette
        Palette should be either a valid :func:`~matplotlib.pyplot.colormaps` string,
        a sequence of colors (in a format that can be understood by matplotlib,
        eg. RGB, RGBS, hex, or a cycler object with key='color'
    """
    from matplotlib.colors import to_hex

    from .palettes import additional_colors

    color_key = f"{value_to_plot}_colors"
    valid = True
    categories = adata.obs[value_to_plot].cat.categories
    length = len(categories)

    if isinstance(palette, str) and "default" in palette:
        palette = palettes.default_26 if length <= 28 else palettes.default_64
    if isinstance(palette, str) and palette in adata.uns:
        palette = (
            [adata.uns[palette][c] for c in categories]
            if isinstance(adata.uns[palette], dict)
            else adata.uns[palette]
        )
    if palette is None and color_key in adata.uns:
        color_keys = adata.uns[color_key]
        # Check if colors already exist in adata.uns and if they are a valid palette
        if isinstance(color_keys, np.ndarray) and isinstance(color_keys[0], dict):
            adata.uns[color_key] = adata.uns[color_key][0]
        # Flatten the dict to a list (mainly for anndata compatibilities)
        if isinstance(adata.uns[color_key], dict):
            adata.uns[color_key] = [adata.uns[color_key][c] for c in categories]
        color_keys = adata.uns[color_key]
        for color in color_keys:
            if not is_color_like(color):
                # check if valid color translate to a hex color value
                if color in additional_colors:
                    color = additional_colors[color]
                else:
                    logg.warn(
                        f"The following color value found in "
                        f"adata.uns['{value_to_plot}_colors'] is not valid: '{color}'. "
                        f"Default colors will be used instead."
                    )
                    valid = False
                    break
        if len(adata.uns[color_key]) < len(adata.obs[value_to_plot].cat.categories):
            valid = False
    elif palette is not None:
        # check is palette given is a valid matplotlib colormap
        if isinstance(palette, str) and palette in pl.colormaps():
            # this creates a palette from a colormap. E.g. 'Accent, Dark2, tab20'
            cmap = pl.get_cmap(palette)
            colors_list = [to_hex(x) for x in cmap(np.linspace(0, 1, length))]

        else:
            # check if palette is an array of length n_obs
            if isinstance(palette, (list, np.ndarray)) or is_categorical(palette):
                if len(adata.obs[value_to_plot]) == len(palette):
                    cats = pd.Categorical(adata.obs[value_to_plot])
                    colors = pd.Categorical(palette)
                    if len(cats) == len(colors):
                        palette = dict(zip(cats, colors))
            # check if palette is as dict and convert it to an ordered list
            if isinstance(palette, dict):
                palette = [palette[c] for c in categories]
            # check if palette is a list and convert it to a cycler
            if isinstance(palette, abc.Sequence):
                if len(palette) < length:
                    logg.warn(
                        "Length of palette colors is smaller than the number of "
                        f"categories (palette length: {len(palette)}, "
                        f"categories length: {length}. "
                        "Some categories will have the same color."
                    )
                # check that colors are valid
                _color_list = []
                for color in palette:
                    if not is_color_like(color):
                        # check if valid color and translate to a hex color value
                        if color in additional_colors:
                            color = additional_colors[color]
                        else:
                            logg.warn(
                                f"The following color value is not valid: '{color}'. "
                                f"Default colors will be used instead."
                            )
                            valid = False
                            break
                    _color_list.append(color)
                palette = cycler(color=_color_list)

            if not isinstance(palette, Cycler) or "color" not in palette.keys:
                logg.warn(
                    "Please check that the value of 'palette' is a valid "
                    "matplotlib colormap string (eg. Set2), a list of color names or "
                    "a cycler with a 'color' key. Default colors will be used instead."
                )
                valid = False

            if valid:
                cc = palette()
                colors_list = [to_hex(next(cc)["color"]) for x in range(length)]
        if valid:
            adata.uns[f"{value_to_plot}_colors"] = colors_list
    else:
        # No valid palette exists or was given
        valid = False

    # Set to defaults:
    if not valid:
        # check if default matplotlib palette has enough colors
        if len(rcParams["axes.prop_cycle"].by_key()["color"]) >= length:
            cc = rcParams["axes.prop_cycle"]()
            palette = [next(cc)["color"] for _ in range(length)]
        # Else fall back to default palettes
        else:
            if length <= 28:
                palette = palettes.default_26
            elif length <= len(palettes.default_64):  # 103 colors
                palette = palettes.default_64
            else:
                palette = ["grey" for _ in range(length)]
                logg.info(
                    f"the obs value {value_to_plot!r} has more than 103 categories. "
                    f"Uniform 'grey' color will be used for all categories."
                )

        adata.uns[f"{value_to_plot}_colors"] = palette[:length]


# TODO: Add docstrings
def set_colorbar(smp, ax, orientation="vertical", labelsize=None):
    """TODO."""
    cax = inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0)
    cb = pl.colorbar(smp, orientation=orientation, cax=cax)
    cb.set_alpha(1)
    cb.ax.tick_params(labelsize=labelsize)
    cb.locator = MaxNLocator(nbins=3, integer=True)
    cb.update_ticks()


# TODO: Add docstrings
def default_palette(palette=None):
    """TODO."""
    if palette is None:
        return rcParams["axes.prop_cycle"]
    elif not isinstance(palette, Cycler):
        return cycler(color=palette)
    else:
        return palette


# TODO: Add docstrings
def adjust_palette(palette, length):
    """TODO."""
    islist = False
    if isinstance(palette, list):
        islist = True
    if (islist and len(palette) < length) or (
        not isinstance(palette, list) and len(palette.by_key()["color"]) < length
    ):
        if length <= 28:
            palette = palettes.default_26
        elif length <= len(palettes.default_64):  # 103 colors
            palette = palettes.default_64
        else:
            palette = ["grey" for _ in range(length)]
            logg.info("more than 103 colors would be required, initializing as 'grey'")
        return palette if islist else cycler(color=palette)
    elif islist:
        return palette
    elif not isinstance(palette, Cycler):
        return cycler(color=palette)
    else:
        return palette


def rgb_custom_colormap(colors=None, alpha=None, N=256):
    """Creates a custom colormap. Colors can be given as names or rgb values.

    Parameters
    ----------
    colors: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    alpha: `list`, `np.ndarray` or `None` (default: `None`)
        Alpha of the colors. Must be same length as colors.
    N: `int` (default: `256`)
        y coordinate

    Returns
    -------
    :class:`~matplotlib.colors.ListedColormap`
    """
    if colors is None:
        colors = ["royalblue", "white", "forestgreen"]
    c = []
    if "transparent" in colors:
        if alpha is None:
            alpha = [1 if i != "transparent" else 0 for i in colors]
        colors = [i if i != "transparent" else "white" for i in colors]

    for color in colors:
        if isinstance(color, str):
            color = to_rgb(color if color.startswith("#") else cnames[color])
            c.append(color)
    if alpha is None:
        alpha = np.ones(len(c))

    vals = np.ones((N, 4))
    ints = len(c) - 1
    n = int(N / ints)

    for j in range(ints):
        start = n * j
        end = n * (j + 1)
        for i in range(3):
            vals[start:end, i] = np.linspace(c[j][i], c[j + 1][i], n)
        vals[start:end, -1] = np.linspace(alpha[j], alpha[j + 1], n)
    return ListedColormap(vals)


"""save figure"""


# TODO: Add docstrings
def savefig_or_show(writekey=None, show=None, dpi=None, ext=None, save=None):
    """TODO."""
    if isinstance(save, str):
        # check whether `save` contains a figure extension
        if ext is None:
            for try_ext in [".svg", ".pdf", ".png"]:
                if save.endswith(try_ext):
                    ext = try_ext[1:]
                    save = save.replace(try_ext, "")
                    break
        # append it
        if "/" in save:
            writekey = None
        writekey = (
            f"{writekey}_{save}" if writekey is not None and len(writekey) > 0 else save
        )
        save = True
    save = settings.autosave if save is None else save
    show = settings.autoshow if show is None else show

    if save:
        if dpi is None:
            # needed in nb b/c internal figures are also influenced by 'savefig.dpi'.
            if (
                not isinstance(rcParams["savefig.dpi"], str)
                and rcParams["savefig.dpi"] < 150
            ):
                if settings._low_resolution_warning:
                    logg.warn(
                        "You are using a low resolution (dpi<150) for saving figures.\n"
                        "Consider running `set_figure_params(dpi_save=...)`, which "
                        "will adjust `matplotlib.rcParams['savefig.dpi']`"
                    )
                    settings._low_resolution_warning = False
            else:
                dpi = rcParams["savefig.dpi"]
        if len(settings.figdir) > 0:
            if settings.figdir[-1] != "/":
                settings.figdir += "/"
            if not os.path.exists(settings.figdir):
                os.makedirs(settings.figdir)
        if ext is None:
            ext = settings.file_format_figs
        filepath = f"{settings.figdir}{settings.plot_prefix}{writekey}"
        if "/" in writekey:
            filepath = f"{writekey}"
        try:
            filename = filepath + f"{settings.plot_suffix}.{ext}"
            pl.savefig(filename, dpi=dpi, bbox_inches="tight")
        except ValueError as e:
            # save as .png if .pdf is not feasible (e.g. specific streamplots)
            filename = filepath + f"{settings.plot_suffix}.png"
            pl.savefig(filename, dpi=dpi, bbox_inches="tight")
            logg_message = (
                f"figure cannot be saved as {ext}, using png instead "
                f"({e.__str__().lower()})."
            )
            logg.msg(logg_message, v=1)
        logg.msg("saving figure to file", filename, v=1)
    if show:
        pl.show()
    if save:
        pl.close()  # clear figure


"""additional plots (linear fit, density, outline, rug)"""


# TODO: Add docstrings
def plot_linfit(
    x,
    y,
    add_linfit=True,
    add_legend=True,
    color=None,
    linewidth=None,
    fontsize=None,
    ax=None,
):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    idx_valid = ~np.isnan(x + y)
    x, y = x[idx_valid], y[idx_valid]

    if isinstance(add_linfit, str) and "no_intercept" in add_linfit:
        mu_x, mu_y = (0, 0)
    else:
        mu_x, mu_y = np.mean(x), np.mean(y)
    slope = (np.mean(x * y) - mu_x * mu_y) / (np.mean(x**2) - mu_x**2)
    offset = mu_y - slope * mu_x

    if isinstance(add_linfit, str) and "intercept" in add_linfit:
        str_list = [a.strip() for a in add_linfit.split(",")]
        if "," in add_linfit:
            add_linfit = [a for a in str_list if "intercept" not in a][0]
        else:
            add_linfit = None

    color = (
        add_linfit
        if isinstance(add_linfit, str)
        else color
        if isinstance(color, str)
        else "k"
    )
    xnew = np.linspace(np.min(x), np.max(x) * 1.02)
    ax.plot(xnew, offset + xnew * slope, linewidth=linewidth, color=color)
    if add_legend:
        kwargs = {"ha": "left", "va": "top", "fontsize": fontsize}
        bbox = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.2}
        txt = r"$\rho = $" + f"{np.round(stats.pearsonr(x, y)[0], 2)}"
        ax.text(0.05, 0.95, txt, transform=ax.transAxes, bbox=bbox, **kwargs)


# TODO: Add docstrings
def plot_polyfit(
    x,
    y,
    add_polyfit=True,
    add_legend=True,
    color=None,
    linewidth=None,
    fontsize=None,
    ax=None,
):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    idx_valid = ~np.isnan(x + y)
    x, y = x[idx_valid], y[idx_valid]
    lb, ub = np.percentile(x, [0.5, 99.5])
    idx = (x > lb) & (x < ub)
    x, y = x[idx], y[idx]
    if isinstance(add_polyfit, str) and "no_intercept" in add_polyfit:
        x, y = np.hstack([np.zeros(len(x)), x]), np.hstack([np.zeros(len(x)), y])
    deg = 2 if isinstance(add_polyfit, (str, bool)) else add_polyfit
    fit = np.polyfit(x, y, deg=deg)
    f = np.poly1d(fit)

    if isinstance(add_polyfit, str) and "intercept" in add_polyfit:
        str_list = [a.strip() for a in add_polyfit.split(",")]
        if "," in add_polyfit:
            add_polyfit = [a for a in str_list if "intercept" not in a][0]
        else:
            add_polyfit = None
    color = (
        add_polyfit
        if isinstance(add_polyfit, str)
        else color
        if isinstance(color, str)
        else "k"
    )

    xnew = np.linspace(np.min(x), np.max(x), num=100)
    ax.plot(xnew, f(xnew), color=color, linewidth=linewidth)

    if add_legend:
        R2 = np.sum((f(x) - np.mean(y)) ** 2) / np.sum((y - np.mean(y)) ** 2)
        kwargs = {"ha": "left", "va": "top", "fontsize": fontsize}
        bbox = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.2}
        txt = r"$R^2 = $" + f"{np.round(R2, 2)}"
        ax.text(0.05, 0.95, txt, transform=ax.transAxes, bbox=bbox, **kwargs)


# TODO: Add docstrings
def plot_vlines(adata, basis, vkey, xkey, linewidth=1, linecolor=None, ax=None):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    xnew = np.linspace(0, np.percentile(make_dense(adata[:, basis].layers[xkey]), 98))
    fits = [
        fit
        for fit in make_unique_list(vkey)
        if all(["velocity" in fit, f"{fit}_gamma" in adata.var.keys()])
    ]
    linecolor, lines = to_list(linecolor), []
    for i, fit in enumerate(fits):
        linestyle = "--" if f"variance_{fit}" in adata.layers.keys() else "-"
        gamma, beta, offset = 1, 1, 0
        if f"{fit}_gamma" in adata.var.keys():
            gamma = adata[:, basis].var[f"{fit}_gamma"].values
        if f"{fit}_beta" in adata.var.keys():
            beta = adata[:, basis].var[f"{fit}_beta"].values
        if f"{fit}_offset" in adata.var.keys():
            offset = adata[:, basis].var[f"{fit}_offset"].values
        if len(linecolor) > i and linecolor[i] is not None:
            c = linecolor[i]
        else:
            c = "k" if i == 0 else None
        kwargs = {"linestyle": linestyle, "linewidth": linewidth}
        (line,) = ax.plot(xnew, gamma / beta * xnew + offset / beta, c=c, **kwargs)
        lines.append(line)
        fits[i] = (
            f"steady-state ratio ({fit})" if len(fits) > 1 else "steady-state ratio"
        )
    return lines, fits


# TODO: Add docstrings
def plot_velocity_fits(
    adata,
    basis,
    vkey=None,
    use_raw=None,
    linewidth=None,
    linecolor=None,
    legend_loc=None,
    legend_fontsize=None,
    show_assignments=None,
    ax=None,
):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    if use_raw is None:
        use_raw = "Ms" not in adata.layers.keys()

    # linear fits
    if vkey is None:
        vkey = "dynamics" if "fit_alpha" in adata.var.keys() else "velocity"
    xkey = "spliced" if use_raw else "Ms"
    lines, fits = plot_vlines(adata, basis, vkey, xkey, linewidth, linecolor, ax=ax)

    # full dynamic fits
    from .simulation import show_full_dynamics

    if "true_alpha" in adata.var.keys() and (
        vkey is not None and "true_dynamics" in vkey
    ):
        line, fit = show_full_dynamics(
            adata,
            basis,
            key="true",
            use_raw=use_raw,
            linewidth=linewidth,
            linecolor=linecolor,
            ax=ax,
        )
        fits.append(fit)
        lines.append(line)
    if "fit_alpha" in adata.var.keys() and (vkey is None or "dynamics" in vkey):
        line, fit = show_full_dynamics(
            adata,
            basis,
            key="fit",
            use_raw=use_raw,
            linewidth=linewidth,
            linecolor=linecolor,
            show_assignments=show_assignments,
            ax=ax,
        )
        fits.append(fit)
        lines.append(line)

    if legend_loc in {True, None, "bottom right"}:
        legend_loc = "lower right"
    if len(fits) > 0 and legend_loc and legend_loc != "none":
        ax.legend(handles=lines, labels=fits, fontsize=legend_fontsize, loc=legend_loc)


# TODO: Add docstrings
def plot_density(
    x, y=None, add_density=True, eval_pts=50, scale=10, alpha=0.3, color="grey", ax=None
):
    """TODO."""
    from scipy.stats import gaussian_kde as kde

    if ax is None:
        ax = pl.gca()
    color = (
        add_density
        if isinstance(add_density, str)
        else color
        if isinstance(color, str)
        else "grey"
    )

    # density plot along x-coordinate
    xnew = np.linspace(min(x), max(x), eval_pts)
    vals = kde(x)(xnew)

    if y is not None:
        offset = max(y) / scale
        scale_s = offset / np.max(vals)
        offset *= 1.3
        vals = vals * scale_s - offset
    else:
        offset = 0

    ax.plot(xnew, vals, color=color)
    ax.fill_between(xnew, -offset, vals, alpha=alpha, color=color)
    ax.set_ylim(-offset)

    # density plot along y-coordinate
    if y is not None:
        ynew = np.linspace(min(y), max(y), eval_pts)
        vals = kde(y)(ynew)

        offset = max(x) / scale
        scale_u = offset / np.max(vals)
        offset *= 1.3
        vals = vals * scale_u - offset

        ax.plot(vals, ynew, color=color)
        ax.fill_betweenx(ynew, -offset, vals, alpha=alpha, color=color)
        ax.set_xlim(-offset)


# TODO: Add docstrings
def plot_outline(
    x, y, kwargs, outline_width=None, outline_color=None, zorder=None, ax=None
):
    """TODO."""
    # Adapted from scanpy. The default outline is a black edge
    # followed by a thin white edged added around connected clusters.
    if ax is None:
        ax = pl.gca()

    bg_width, gp_width = (0.3, 0.05) if outline_width is None else outline_width
    bg_color, gp_color = ("black", "white") if outline_color is None else outline_color

    s = kwargs.pop("s")
    point = np.sqrt(s)

    gp_size = (2 * (point * gp_width) + point) ** 2
    bg_size = (2 * (point * bg_width) + np.sqrt(gp_size)) ** 2

    kwargs["edgecolor"] = "none"
    zord = 0 if zorder is None else zorder
    if "rasterized" not in kwargs:
        kwargs["rasterized"] = settings._vector_friendly
    ax.scatter(x, y, s=bg_size, marker=".", c=bg_color, zorder=zord - 2, **kwargs)
    ax.scatter(x, y, s=gp_size, marker=".", c=gp_color, zorder=zord - 1, **kwargs)
    # restore size
    kwargs["s"] = s


# TODO: Add docstrings
def plot_rug(x, height=0.03, color=None, ax=None, **kwargs):
    """TODO."""
    if ax is None:
        ax = pl.gca()
    x = np.asarray(x)

    transform = tx.blended_transform_factory(ax.transData, ax.transAxes)
    line_segs = np.column_stack(
        [np.repeat(x, 2), np.tile([0, height], len(x))]
    ).reshape([len(x), 2, 2])
    kwargs.update({"transform": transform, "color": color})
    ax.add_collection(LineCollection(line_segs, **kwargs))
    ax.autoscale_view(scalex=True, scaley=False)


"""for velocity_embedding"""


# TODO: Add docstrings
def velocity_embedding_changed(adata, basis, vkey):
    """TODO."""
    if f"X_{basis}" not in adata.obsm.keys():
        changed = False
    else:
        changed = f"{vkey}_{basis}" not in adata.obsm_keys()
        if f"{vkey}_params" in adata.uns.keys():
            sett = adata.uns[f"{vkey}_params"]
            changed |= "embeddings" not in sett or basis not in sett["embeddings"]
    return changed


"""additional plots (linear fit, density, outline, rug)"""


def hist(
    arrays,
    alpha=0.5,
    bins=50,
    color=None,
    colors=None,
    labels=None,
    hist=None,
    kde=None,
    bw_method=None,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    cutoff=None,
    xscale=None,
    yscale=None,
    xticks=None,
    yticks=None,
    fontsize=None,
    legend_fontsize=None,
    figsize=None,
    normed=None,
    perc=None,
    exclude_zeros=None,
    axvline=None,
    axhline=None,
    pdf=None,
    ax=None,
    dpi=None,
    show=True,
    save=None,
    **kwargs,
):
    """Plot a histogram.

    Parameters
    ----------
    arrays: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    alpha: `list`, `np.ndarray` or `None` (default: `None`)
        Alpha of the colors. Must be same length as colors.
    bins : `int` or `sequence` (default: `50`)
        If an integer is given, ``bins + 1`` bin edges are calculated and
        returned, consistent with `numpy.histogram`.
        If `bins` is a sequence, gives bin edges, including left edge of
        first bin and right edge of last bin.  In this case, `bins` is
        returned unmodified.
    colors: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    labels : str or None (default: `None`)
            String, or sequence of strings to label clusters.
    hist: `bool` or `None` (default: `None`)
        Whether to show histogram.
    kde: `bool` or `None` (default: `None`)
        Whether to use kernel density estimation on data.
    bw_method : `str`, `scalar` or `callable`, (default: `None`)
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a
        scalar, this will be used directly as `kde.factor`.  If a callable,
        it should take a `gaussian_kde` instance as only parameter and
        return a scalar.  If None (default), nothing happens; the current
        `kde.covariance_factor` method is kept.
    xlabel: `str` (default: `None`)
        Label of x-axis.
    ylabel: `str` (default: `None`)
        Label of y-axis.
    xlim: tuple, e.g. [0,1] or `None` (default: `None`)
        Restrict x-limits of the axis.
    ylim: tuple, e.g. [0,1] or `None` (default: `None`)
        Restrict y-limits of the axis.
    cutoff: tuple, e.g. [0,1] or `float` or `None` (default: `None`)
        Bins will be cut off below and above the cutoff values.
    xscale: `log` or `None` (default: `None`)
        Scale of the x-axis.
    yscale: `log` or `None` (default: `None`)
        Scale of the y-axis.
    fontsize: `float` (default: `None`)
        Label font size.
    legend_fontsize: `int` (default: `None`)
        Legend font size.
    figsize: tuple (default: `(7,5)`)
        Figure size.
    normed: `bool` or `None` (default: `None`)
        Whether to normalize data.
    perc: tuple, e.g. [2,98] (default: `None`)
        Specify percentile for continuous coloring.
    exclude_zeros: `bool` or `None` (default: `None`)
        Whether to exclude zeros in data for the kde and hist plot.
    axvline: `float` or `None` (default: `None`)
        Plot a vertical line at the specified x-value.
    axhline `float` or `None` (default: `None`)
        Plot a horizontal line at the specified y-value.
    pdf: `str` or `None` (default: `None`)
        probability density function to be fitted,
        e.g., 'norm', 't', 'chi', 'beta', 'gamma', 'laplace' etc.
    ax: `matplotlib.Axes`, optional (default: `None`)
        A matplotlib axes object. Only works if plotting a single component.
    dpi: `int` (default: 80)
        Figure dpi.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the default filename.
        Infer the filetype if ending on {'.pdf', '.png', '.svg'}.

    Returns
    -------
    If `show==False` a `matplotlib.Axis`
    """
    if ax is None:
        fig, ax = pl.subplots(figsize=figsize, dpi=dpi)

    arrays = (
        [np.ravel(array) for array in arrays]
        if isinstance(arrays, (list, tuple)) or arrays.ndim > 1
        else [arrays]
    )
    if normed is None:
        normed = kde or pdf
    if hist is None:
        hist = not kde

    palette = default_palette(None).by_key()["color"][::-1]
    colors = to_list(colors) if color is None else to_list(color)
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    masked_arrays = np.ma.masked_invalid(np.hstack(arrays)).compressed()
    bmin, bmax = masked_arrays.min(), masked_arrays.max()
    if xlim is not None:
        bmin, bmax = max(bmin, xlim[0]), min(bmax, xlim[1])
    elif perc is not None:
        if np.size(perc) < 2:
            perc = [perc, 100] if perc < 50 else [0, perc]
        bmin, bmax = np.nanpercentile(masked_arrays, perc)
    bins = np.arange(bmin, bmax + (bmax - bmin) / bins, (bmax - bmin) / bins)

    if xscale == "log":
        bins = bins[bins > 0]
        bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

    if cutoff is not None:
        bins = (
            bins[(bins > cutoff[0]) & (bins < cutoff[1])]
            if isinstance(cutoff, list)
            else bins[bins < cutoff]
        )

    if isinstance(labels, str):
        labels = [labels]

    for i, x in enumerate(arrays):
        x_vals = np.array(x[np.isfinite(x)])
        if exclude_zeros:
            x_vals = np.array(x_vals[x_vals != 0])
        if kde:
            from scipy.stats import gaussian_kde

            kde_bins = gaussian_kde(x_vals, bw_method=bw_method)(bins)
            if not normed:
                kde_bins *= (bins[1] - bins[0]) * len(x_vals)
            ci, li = colors[i], labels[i] if labels is not None else None
            ax.plot(bins, kde_bins, color=ci)
            ax.fill_between(bins, 0, kde_bins, alpha=0.4, color=ci, label=li)
            ylim = np.min(kde_bins) if ylim is None else ylim
        if hist:
            ci, li = colors[i], labels[i] if labels is not None and not kde else None
            kwargs.update({"color": ci, "label": li})
            # TODO: Check if proper exception is used and potentially remove
            try:
                ax.hist(x_vals, bins=bins, alpha=alpha, density=normed, **kwargs)
            except TypeError:
                ax.hist(x_vals, bins=bins, alpha=alpha, **kwargs)
    if xlabel is None:
        xlabel = ""
    if ylabel is None:
        ylabel = ""
    set_label(xlabel, ylabel, fontsize=fontsize, ax=ax)

    if labels is not None:
        ax.legend(fontsize=legend_fontsize)

    if axvline:
        ax.axvline(axvline)
    if axhline:
        ax.axhline(axhline)

    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    update_axes(ax, xlim, ylim, fontsize, frameon=True)

    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    @pl.FuncFormatter
    def log_fmt(x, pos):
        return r"${%.2g}$" % (x)

    if xscale == "log":
        if xticks is None:
            lspace = np.linspace(-10, 10, 21)
            ticks = [a for a in [10**a for a in lspace] if bmin < a < bmax]
            ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(log_fmt)
        ax.minorticks_off()

    pdf = [pdf] if isinstance(pdf, str) else pdf
    if pdf is not None:
        fits = []
        for i, pdf_name in enumerate(pdf):
            xt = ax.get_xticks()
            xmin, xmax = min(xt), max(xt)
            lnspc = np.linspace(xmin, xmax, len(bins))

            if "(" in pdf_name:  # used passed parameters
                start = pdf_name.rfind("(")
                args, pdf_name = eval(pdf_name[start:]), pdf_name[:start]
            else:  # fit parameters
                args = getattr(stats, pdf_name).fit(x_vals)
            pd_vals = getattr(stats, pdf_name).pdf(lnspc, *args)
            logg.info("Fitting", pdf_name, np.round(args, 4), ".")
            fit = ax.plot(lnspc, pd_vals, label=pdf_name, color=colors[i])
            fits.extend(fit)
        ax.legend(handles=fits, labels=pdf, fontsize=legend_fontsize)

    if rcParams["savefig.transparent"]:
        ax.patch.set_alpha(0)

    savefig_or_show(dpi=dpi, save=save, show=show)
    if show is False:
        return ax


# TODO: Add docstrings
def plot(
    arrays,
    normalize=False,
    colors=None,
    labels=None,
    xlabel=None,
    ylabel=None,
    xscale=None,
    yscale=None,
    ax=None,
    figsize=None,
    dpi=None,
    show=True,
):
    """TODO."""
    ax, show = get_ax(ax, show, figsize, dpi)
    arrays = np.array(arrays)
    arrays = (
        arrays if isinstance(arrays, (list, tuple)) or arrays.ndim > 1 else [arrays]
    )

    palette = default_palette(None).by_key()["color"][::-1]
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    for i, array in enumerate(arrays):
        X = array[np.isfinite(array)]
        X = X / np.max(X) if normalize else X
        ax.plot(X, color=colors[i], label=labels[i] if labels is not None else None)

    ax.set_xlabel(xlabel if xlabel is not None else "")
    ax.set_ylabel(ylabel if xlabel is not None else "")
    if labels is not None:
        ax.legend()
    if xscale is not None:
        ax.xscale(xscale)
    if yscale is not None:
        ax.yscale(yscale)

    if not show:
        return ax


# TODO: Add docstrings
def fraction_timeseries(
    adata,
    xkey="clusters",
    tkey="dpt_pseudotime",
    bins=30,
    legend_loc="best",
    title=None,
    fontsize=None,
    ax=None,
    figsize=None,
    dpi=None,
    xlabel=None,
    ylabel=None,
    show=True,
):
    """TODO."""
    t = np.linspace(0, 1 + 1 / bins, bins)
    types = np.unique(adata.obs[xkey].values)

    y = []
    for i in range(bins - 1):
        mask = np.all(
            [adata.obs[tkey].values <= t[i + 1], adata.obs[tkey].values > t[i]], axis=0
        )
        x = list(adata[mask].obs[xkey].values)
        y.append([])
        for name in types:
            occur = x.count(name)
            y[-1].append(occur)
        y[-1] /= np.sum(y[-1])
    y = np.array(y).T

    ax = pl.figure(figsize=figsize, dpi=dpi) if ax is None else ax

    c = None
    if "clusters_colors" in adata.uns.keys():
        c = adata.uns["clusters_colors"]
    pl.stackplot(t[:-1], y, baseline="zero", labels=types, colors=c, edgecolor="white")

    pl.legend(types, loc=legend_loc)
    if title is not None:
        pl.title(title, fontsize=fontsize)
    pl.xlabel(tkey if xlabel is None else xlabel, fontsize=fontsize)
    pl.ylabel(f"{xkey} fractions" if ylabel is None else ylabel, fontsize=fontsize)
    pl.xlim(adata.obs[tkey].values.min(), adata.obs[tkey].values.max())
    pl.ylim(0, 1)

    if not show:
        return ax
    else:
        pl.show()


"""deprecated"""


# TODO: Add docstrings
def make_unique_list(key, allow_array=False):
    """TODO."""
    if isinstance(key, (Index, abc.KeysView)):
        key = list(key)
    is_list = (
        isinstance(key, (list, tuple, np.record))
        if allow_array
        else isinstance(key, (list, tuple, np.ndarray, np.record))
    )
    is_list_of_str = is_list and all(isinstance(item, str) for item in key)
    return key if is_list_of_str else key if is_list and len(key) < 20 else [key]


# TODO: Add docstrings
def make_unique_valid_list(adata, keys):
    """TODO."""
    keys = make_unique_list(keys)
    if all(isinstance(item, str) for item in keys):
        for i, key in enumerate(keys):
            if key.startswith("X_"):
                keys[i] = key = key[2:]
            check_basis(adata, key)
        valid_keys = np.hstack(
            [
                adata.obs.keys(),
                adata.var.keys(),
                adata.varm.keys(),
                adata.obsm.keys(),
                [key[2:] for key in adata.obsm.keys()],
                list(adata.layers.keys()),
            ]
        )
        keys_ = keys
        keys = [key for key in keys if key in valid_keys or key in adata.var_names]
        keys_ = [key for key in keys_ if key not in keys]
        if len(keys_) > 0:
            logg.warn(", ".join(keys_), "not found.")
    return keys


# TODO: Add docstrings
def get_temporal_connectivities(adata, tkey, n_convolve=30):
    """TODO."""
    from scvelo.tools.utils import normalize
    from scvelo.tools.velocity_graph import vals_to_csr

    # from scvelo.tools.utils import get_indices
    # c_idx = get_indices(get_connectivities(adata, recurse_neighbors=True))[0]
    # lspace = np.linspace(0, len(c_idx) - 1, len(c_idx), dtype=int)
    # c_idx = np.hstack([c_idx, lspace[:, None]])

    rows, cols, vals, n_obs, n_convolve = [], [], [], len(tkey), int(n_convolve / 2)
    for i in range(n_obs):
        i_max = None if i + n_convolve >= n_obs else i + n_convolve
        i_min = np.max([0, i - n_convolve])
        t_idx = np.argsort(tkey)[i_min:i_max]  # temporal neighbourhood
        # t_idx = np.intersect1d(t_idx, c_idx[i])
        rows.extend(np.ones(len(t_idx), dtype=int) * np.argsort(tkey)[i])
        cols.extend(t_idx)
        vals.extend(np.ones(len(t_idx), dtype=int))

    # c_conn = get_connectivities(adata, recurse_neighbors=True)
    t_conn = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs))
    return normalize(t_conn)  # normalize(t_conn.multiply(c_conn))
