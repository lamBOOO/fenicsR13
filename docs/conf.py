# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('../tests/heat'))
sys.path.insert(0, os.path.abspath('../tests/stress'))
sys.path.insert(0, os.path.abspath('../tests/r13'))
sys.path.insert(0, os.path.abspath('../examples'))

# -- Project information -----------------------------------------------------

project = 'fenicsR13'
copyright = '2019'
author = 'Lambert Theisen'

# The full version, including alpha/beta/rc tags
release = '0.5'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  'sphinx.ext.autodoc',
  'sphinx.ext.coverage',
  'sphinx.ext.autosummary',
  'sphinx.ext.mathjax',
  'sphinx.ext.viewcode',
  'sphinx.ext.napoleon', # for numpy and Google docstrings
  'fluiddoc.mathmacro',
  'sphinx.ext.imgconverter', # for svg usage in latex, but bad output
  'sphinx.ext.graphviz',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

today_fmt = '%b %d %y at %H:%M' # for |today| directive in index.rst

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = "../media/logo_large.png"
html_favicon = "../media/logo_large.png"


html_theme_options = {
    'sidebar_collapse': False,
    'description': "Release v{}<br>Extended gas dynamics using FEniCS platform.".format(release),
    "touch_icon": "../media/logo_large.png",
    "fixed_sidebar": False, # fails on mobile and with large sidebar
    "note_bg": "#FFF59C",
    "show_relbars": False,
}

# -- LaTeX
latex_logo = "../media/logo_large.png"
latex_elements = {
    'extraclassoptions': 'openany', # skip empty pages
    # 'printindex': '\\footnotesize\\raggedright\\printindex', # small index
    'preamble': r'''
      \usepackage[columns=1]{idxlayout}
    ''' # for onecolumn index layout
}

# -- Autodoc
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'undoc-members': True,
    'private-members': True,
    'special-members': True,
    # 'inherited-members': True,
    'show-inheritance': True,
    # 'ignore-module-all': True,
    # 'imported-members': True,
    'exclude-members': '__dict__,__module__,__weakref__'
}
