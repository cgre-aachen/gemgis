# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# -- Project information -----------------------------------------------------

project = "GemGIS"
copyright = "2020–2023, GemGIS Developers"
author = 'Alexander Juestel'



# -- GemGIS configuration ---------------------------------------------------
import os
import sys
#sys.path.insert(0, os.path.abspath(".."))
#sys.path.insert(0, os.path.abspath('../..'))
#sys.path.insert(0, os.path.abspath("gemgis/"))
#sys.path.insert(0, os.path.abspath("../gemgis/"))
#sys.path.insert(0, os.path.abspath("../../gemgis/"))

# The full version, including alpha/beta/rc tags

import gemgis as gg

version = gg.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'sphinx_book_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx_markdown_tables',
    'sphinx_copybutton',
    'sphinx.ext.extlinks',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

html_title = "GemGIS - Spatial data processing for geomodeling"
html_logo = "getting_started/images/Modern1.png"
html_favicon = "https://raw.githubusercontent.com/cgre-aachen/gemgis/main/docs/getting_started/images/favicon.ico"

html_theme_options = {
    "navbar_end": ["navbar-icon-links.html", "search-field.html"],

}

nbsphinx_execute = 'never'

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]

copybutton_prompt_text = ">>> "

#autodoc_mock_imports = ["gemgis"]
