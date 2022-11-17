#!/usr/bin/env python
# https://github.com/kennethreitz/setup.py/blob/master/setup.py

import os

from setuptools import find_packages, setup

# Package meta-data.
NAME = "spgrep"
DESCRIPTION = "On-the-fly generator of space-group irreducible representations"
URL = "https://github.com/spglib/spgrep"
AUTHOR = "Kohei Shinohara"
EMAIL = "kshinohara0508@gmail.com"
REQUIRES_PYTHON = ">=3.8.0"

# What packages are required for this module to be executed?
REQUIRED = [
    "setuptools",
    "setuptools_scm",
    "wheel",
    "typing_extensions",
    "numpy>=1.20.1",
    "spglib>=1.16.5",
]

# What packages are optional?
EXTRAS = {
    "dev": [
        "pytest==7.2.0",
        "pytest-cov==4.0.0",
        "pre-commit",
        "black",
        "mypy",
        "flake8",
        "pyupgrade",
        "pydocstyle",
        "nbqa",
        "phonopy==2.16.3",
        # Jupyter notebook
        "notebook",
        "matplotlib",
        "seaborn",
        "ipython",
        "ipykernel",
    ],
    "docs": [
        "sphinx",
        "sphinx-autobuild",
        "nbsphinx",
        "sphinxcontrib-bibtex",
        "sphinxcontrib-mermaid",
        "myst-parser",
        "sphinx-book-theme",
        "linkify-it-py",
    ],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


# Where the magic happens:
setup(
    name=NAME,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["spgrep"]),
    package_data={},
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    # numpy: https://github.com/numpy/numpy/issues/2434
    setup_requires=["setuptools_scm", "numpy"],
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="BSD",
    test_suite="tests",
    zip_safe=False,
    use_scm_version=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
