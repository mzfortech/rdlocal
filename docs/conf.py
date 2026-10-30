from rdlocal import __version__

project = "rdlocal"
copyright = "2026, Model Zero for Technological Solutions"
author = "Model Zero for Technological Solutions"
release = __version__
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
myst_enable_extensions = ["dollarmath"]
myst_heading_anchors = 3
autodoc_typehints = "description"
autodoc_member_order = "bysource"
html_theme = "alabaster"
exclude_patterns = ["_build"]
