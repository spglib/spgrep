# Development

## Installation

Install [mise](https://mise.jdx.dev/getting-started.html), [uv](https://docs.astral.sh/uv/getting-started/installation/), and [prek](https://prek.j178.dev/installation/) before running these commands.

```shell
git clone git@github.com:spglib/spgrep.git
cd spgrep
mise trust
mise run install
prek install

# Run prek manually
mise run prek

# Run tests
mise run test

# List available tasks
mise run
```

## Build and preview the documentation

The site uses Zensical and mkdocstrings. Generate the introduction from `README.md`,
the example pages from the notebooks' saved outputs, and the bibliography from
`docs/references.bib` before building:

```shell
mise run docs-build
mise run docs
```

Open the URL printed by `zensical serve`. Rerun `scripts/prepare_docs.py` after
editing the README, notebooks, or bibliography. The generated Markdown and images
are ignored by Git. Notebook cells are not executed during the documentation build;
run and save a notebook to update its displayed results.

Without mise, run `uv sync --extra docs`, `uv run python scripts/prepare_docs.py`,
then `uv run zensical build --clean --strict` or `uv run zensical serve`.

Write docstrings in NumPy style with Markdown links and dollar-delimited math.
The API reference is generated directly from `src/spgrep`.

### How to generate diagram from Mermaid file

https://www.npmjs.com/package/@mermaid-js/mermaid-cli/v/8.9.2

```shell
npm install @mermaid-js/mermaid-cli
./node_modules/.bin/mmdc --input docs/point_group_chain.mmd
```

### How to compile JOSS draft

```shell
# At root directory
docker run --rm --volume $PWD/docs/paper:/data --user $(id -u):$(id -g) --env JOURNAL=joss openjournals/inara
```

## Release

[tagpr](https://github.com/Songmu/tagpr) maintains a release pull request after pushes
to `develop`. Merging that pull request creates a `vX.Y.Z` tag and starts
`.github/workflows/deploy.yml`. The deploy workflow publishes the package to PyPI,
then creates the GitHub Release with generated release notes.

The package version comes from Git tags through `setuptools-scm`, so tagpr does not
update a version file. The curated `docs/changelog.md` remains manually maintained.

### Configure the release app

1. Install a GitHub App on `spglib/spgrep` with **Contents: Read and write**,
   **Pull requests: Read and write**, and **Issues: Read-only** permissions.
   An existing release App can be reused if its installation includes this repository.
2. Set the Actions repository variable `TAGPR_APP_CLIENT_ID` to the App's client ID.
3. Set the Actions repository secret `TAGPR_APP_PRIVATE_KEY` to the App's private key.
4. Run the **Prepare release** workflow on `develop` to prepare the first release PR.

The workflow skips the job until `TAGPR_APP_CLIENT_ID` is set. The App token lets
release PRs trigger CI and tags trigger the separate deploy workflow. Tags created
with the default `GITHUB_TOKEN` do not trigger that workflow.

PyPI publishing continues to use Trusted Publishing through the `release`
environment and `deploy.yml`.

### Publish a release

1. Review the release PR's proposed version. The default is a patch release.
   For a minor or major release, add `tagpr:minor` or `tagpr:major` to the release PR,
   then rerun **Prepare release** on `develop` to refresh the proposal.
   In tag-only mode, these labels determine the final version when the PR is merged.
2. Update `docs/changelog.md` with a dated entry for that version when needed.
   Commit the entry to the release PR branch or merge it into `develop` and let tagpr
   refresh the release PR.
3. Wait for the release PR's CI to pass, then merge it.
4. Check that **Prepare release** creates the tag and **deploy** publishes the package
   and GitHub Release.

For a release candidate, retain the manual `vX.Y.Z-rcN` tag workflow. These tags
publish to TestPyPI and create a prerelease on GitHub.

## Subpages

- [Memo](memo.md)
- [v0.8.0 performance comparison](performance-v0.8.0.md)
