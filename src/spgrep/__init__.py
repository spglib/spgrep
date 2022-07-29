"""Import top APIs and version."""
from importlib.metadata import PackageNotFoundError, version

from spgrep.core import (  # noqa: F401
    get_crystallographic_pointgroup_irreps_from_symmetry,
    get_spacegroup_irreps,
    get_spacegroup_irreps_from_primitive_symmetry,
)

# https://github.com/pypa/setuptools_scm/#retrieving-package-version-at-runtime
try:
    __version__ = version("spgrep")
except PackageNotFoundError:
    # package is not installed
    pass
