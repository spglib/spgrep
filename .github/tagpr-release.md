Release {{.NextVersion}}

Review the proposed version and update `docs/changelog.md` with a dated release entry when needed. Wait for CI to pass before merging.

Merging this pull request creates the release tag and starts `deploy.yml`, which publishes to PyPI and then creates the GitHub Release. The package version comes from the tag through `setuptools-scm`.
