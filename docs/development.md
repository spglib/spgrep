# Development

## Installation

```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
git clone git@github.com:spglib/spgrep.git
cd spgrep
pip install -e ".[dev,docs]"
pre-commit install

# Run pre-commit manually
pre-commit run --all-file 
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

## Release

```shell
# Confirm the version number via `setuptools-scm`
python -m setuptools_scm

# Update changelog here
vim docs/changelog.md

# Push with tag
git tag <next-version>
git push origin main
git push origin --tags
```
