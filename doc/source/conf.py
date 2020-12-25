# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import sphinx_rtd_theme
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'ComCH'
copyright = '2020, Anibal M. Medina-Mardones'
author = 'Anibal M. Medina-Mardones'

# The full version, including alpha/beta/rc tags
# release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    #'nbsphinx_link',
    #'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    # 'numpydoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    # 'sphinx.ext.imgconverter',
    # 'sphinx_issues',
    'sphinx_rtd_theme',
    'sphinx.ext.napoleon'
    # 'custom_references_resolver' # custom for sklearn, not sure what it does
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# generate autosummary even if no references
autosummary_generate = True

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'logo_only': True,
}
    
# autodoc options
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': True,
    'undoc-members': True,
    'exclude-members': '__weakref__, __dict__, __module__, __hash__, __str__, __init__'
}
