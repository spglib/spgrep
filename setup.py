from pathlib import Path

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
    "spglib>=2.0.2",
]

# What packages are optional?
EXTRAS = {
    "dev": [
        "pytest==7.3.1",
        "pytest-cov==4.1.0",
        "pre-commit",
        "black",
        "mypy",
        "flake8",
        "pyupgrade",
        "pydocstyle",
        "nbqa",
        "phonopy==2.19.1",
        # Jupyter notebook
        "notebook",
        "matplotlib",
        "seaborn",
        "ipython",
        "ipykernel",
    ],
    "docs": [
        # We cannot update sphinx to >6 until docutils==0.20 is released: https://github.com/mcmtroffaes/sphinxcontrib-bibtex/issues/322
        # Another issue requires sphinx<4.3 (fixed in pydata-sphinx-theme>=0.8.0): https://github.com/pydata/pydata-sphinx-theme/pull/509
        "docutils<=0.17",
        "sphinx<4.3",
        "sphinx-autobuild",
        "nbsphinx",
        "sphinxcontrib-bibtex",
        "myst-parser",
        "sphinx-book-theme",
        "linkify-it-py",
    ],
}

# Import the README and use it as the long-description.
try:
    with open(Path(__file__).parent.resolve() / "README.md", encoding="utf-8") as f:
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
        "Development Status :: 3 - Alpha",
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
