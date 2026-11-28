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

## How to compile documents

```shell
mise run docs
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
  v0.8.0 performance comparison <performance-v0.8.0>
```
