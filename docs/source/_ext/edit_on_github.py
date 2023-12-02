"""Loosely based on gist.github.com/MantasVaitkunas/7c16de233812adcb7028."""

import os
import warnings

__licence__ = "BSD (3 clause)"


# TODO: Add docstrings
def get_github_repo(app, path):
    """TODO."""
    if path.endswith(".ipynb"):
        return app.config.github_nb_repo, "/"
    return app.config.github_repo, "/docs/source/"


# TODO: Add docstrings
def html_page_context(app, pagename, templatename, context, doctree):
    """TODO."""
    if templatename != "page.html":
        return

    if not app.config.github_repo:
        warnings.warn("`github_repo `not specified", stacklevel=2)
        return

    if not app.config.github_nb_repo:
        nb_repo = f"{app.config.github_repo}_notebooks"
        warnings.warn(
            "`github_nb_repo `not specified. Setting to `{nb_repo}`", stacklevel=2
        )
        app.config.github_nb_repo = nb_repo

    path = os.path.relpath(doctree.get("source"), app.builder.srcdir)
    repo, conf_py_path = get_github_repo(app, path)

    # For sphinx_rtd_theme.
    context["display_github"] = True
    context["github_user"] = "theislab"
    context["github_version"] = "main"
    context["github_repo"] = repo
    context["conf_py_path"] = conf_py_path


# TODO: Add docstrings
def setup(app):
    """TODO."""
    app.add_config_value("github_nb_repo", "", True)
    app.add_config_value("github_repo", "", True)
    app.connect("html-page-context", html_page_context)
