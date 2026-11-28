"""Generate the introduction, notebook pages, and bibliography before running Zensical."""

from pathlib import Path

from nbconvert import MarkdownExporter
from pybtex.database import parse_file
from pybtex.plugin import find_plugin

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def main() -> None:
    """Render saved notebook outputs without executing the examples."""
    # README links are relative to the repository root, while index.md is in docs/.
    introduction = (ROOT / "README.md").read_text().replace("](docs/", "](")
    (DOCS / "index.md").write_text(introduction)

    exporter = MarkdownExporter()
    for notebook in sorted((DOCS / "examples").glob("*.ipynb")):
        body, resources = exporter.from_filename(
            str(notebook),
            resources={"output_files_dir": f"{notebook.stem}_files"},
        )
        title, content = body.split("\n", 1)
        download = (
            f"[Download this notebook]({notebook.name}). Outputs below are saved in the notebook."
        )
        notebook.with_suffix(".md").write_text(f"{title}\n\n{download}\n{content}")
        for filename, data in resources["outputs"].items():
            output = notebook.parent / filename
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)

    bibliography = parse_file(str(DOCS / "references.bib"))
    style = find_plugin("pybtex.style.formatting", "plain")()
    backend = find_plugin("pybtex.backends", "html")()
    references = ["# References", "", "[Download the bibliography](references.bib).", ""]
    for entry in style.format_bibliography(bibliography):
        references.extend(
            [f'<p id="{entry.key}">[{entry.label}] {entry.text.render(backend)}</p>', ""]
        )
    (DOCS / "references.md").write_text("\n".join(references))


if __name__ == "__main__":
    main()
