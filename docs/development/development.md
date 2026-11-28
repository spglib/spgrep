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

```shell
# Confirm the version number via `setuptools-scm`
python -m setuptools_scm

# Update changelog here
vim docs/changelog.md

# Push with tag
git tag <next-version>
git push origin <next-version>
```

## Subpages

- [Memo](memo.md)
- [v0.8.0 performance comparison](performance-v0.8.0.md)
