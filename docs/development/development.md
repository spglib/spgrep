# Development

## Installation

```shell
git clone git@github.com:spglib/spgrep.git
cd spgrep
uv sync --all-extras
prek install

# Run prek manually
prek run --all-file
```

## How to compile documents

```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

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

```{toctree}
  :maxdepth: 1
  Memo <memo>
```
