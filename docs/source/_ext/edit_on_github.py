"""
Loosely based on gist.github.com/MantasVaitkunas/7c16de233812adcb7028
"""

import os
import warnings

__licence__ = "BSD (3 clause)"


def get_github_repo(app, path):
    if path.endswith(".ipynb"):
        return app.config.github_nb_repo, "/"
    return app.config.github_repo, "/docs/source/"


def html_page_context(app, pagename, templatename, context, doctree):
    if templatename != "page.html":
        return

    if not app.config.github_repo:
        warnings.warn("`github_repo `not specified")
        return

    if not app.config.github_nb_repo:
        nb_repo = f"{app.config.github_repo}_notebooks"
        warnings.warn("`github_nb_repo `not specified. Setting to `{nb_repo}`")
        app.config.github_nb_repo = nb_repo

    path = os.path.relpath(doctree.get("source"), app.builder.srcdir)
    repo, conf_py_path = get_github_repo(app, path)

    # For sphinx_rtd_theme.
    context["display_github"] = True
    context["github_user"] = "theislab"
    context["github_version"] = "master"
    context["github_repo"] = repo
    context["conf_py_path"] = conf_py_path


def setup(app):
    app.add_config_value("github_nb_repo", "", True)
    app.add_config_value("github_repo", "", True)
    app.connect("html-page-context", html_page_context)
