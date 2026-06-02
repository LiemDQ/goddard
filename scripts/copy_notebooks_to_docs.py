#!/usr/bin/env python3
"""Copy Jupyter notebooks from ``examples/`` into ``docs/examples/``.

mkdocs-jupyter renders ``.ipynb`` files that live under the docs source tree.
Rather than symlinking (fragile across platforms) we copy them so the docs build
is self-contained.

Run this as part of the docs-prepare step before ``mkdocs build``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "examples"
DST_DIR = PROJECT_ROOT / "docs" / "examples"


def main() -> None:
    DST_DIR.mkdir(parents=True, exist_ok=True)

    notebooks = sorted(SRC_DIR.glob("*.ipynb"))
    if not notebooks:
        print("[copy_notebooks] no .ipynb files found in examples/")
        return

    for nb in notebooks:
        dst = DST_DIR / nb.name
        shutil.copy2(nb, dst)
        print(f"[copy_notebooks] {nb.relative_to(PROJECT_ROOT)} -> {dst.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
