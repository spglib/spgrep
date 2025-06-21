set shell := ["zsh", "-uc"]
set positional-arguments

default:
    just --list

pre-commit:
    pre-commit run --all-files

install:
    uv sync --all-extras

test:
    uv run pytest -v

docs:
    uv run sphinx-autobuild docs docs_build
