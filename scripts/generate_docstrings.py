#!/usr/bin/env python3
"""Generate ``goddard_docstrings.h`` from Goddard's public C++ headers.

Thin wrapper around `pybind11_mkdoc <https://github.com/pybind/pybind11_mkdoc>`_.
It parses the requested headers with libclang, extracts Doxygen-style comments,
and emits a header containing ``__doc_*`` string constants plus a ``DOC(...)``
macro that nanobind binding files can use::

    .def("set_state_TP", &Goddard::Gas::set_state_TP,
         "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))

If ``pybind11_mkdoc`` (or its libclang dependency) is not available or fails, the
script keeps an existing header and writes a stub only when there is none; ``--stub``
always writes the stub, which defines ``DOC(...)`` as an empty string. The generated
header is committed, so builds
without the docs toolchain (e.g. plain ``pixi run compile``) use it as is and
only fall back to the stub when it is missing (``--stub --if-missing``).
"""

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path

STUB_CONTENTS = """\
#pragma once
// Stub goddard_docstrings.h emitted by scripts/generate_docstrings.py.
//
// Generated when pybind11_mkdoc is unavailable, or when the committed header is missing
// and -Dgoddard_BUILD_DOCSTRINGS=OFF. DOC(...) expands to an empty string so binding
// sources compile cleanly without embedded docstrings. To regenerate the real docstrings,
// run `pixi run -e docs docs-compile` (it configures with -Dgoddard_BUILD_DOCSTRINGS=ON).

#define DOC(...) ""
"""


def write_stub(out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(STUB_CONTENTS)
    print(f"[generate_docstrings] wrote stub -> {out_path}", file=sys.stderr)


def fall_back(out_path: Path, reason: str) -> None:
    """Keep an existing header when extraction is impossible; write the stub only if none exists.

    The header is committed, so replacing it with the stub would silently strip every docstring
    from the next commit and from packages built from it.
    """
    if out_path.exists():
        print(f"[generate_docstrings] {reason}; keeping the existing {out_path}", file=sys.stderr)
    else:
        print(f"[generate_docstrings] {reason}; falling back to stub", file=sys.stderr)
        write_stub(out_path)


_DOCSTRING = re.compile(r'R"doc\((.*?)\)doc"', re.DOTALL)
_LEADING_MARKER = re.compile(r"^(?://[/!]<|[/!]?<)\s?")


def clean_docstrings(text: str) -> str:
    """Strip comment markers pybind11_mkdoc leaves in: the closing ``*/`` of a one-line
    ``/** ... */`` comment and the ``<`` of a trailing ``///<`` or ``//!<`` comment."""
    def clean(match: re.Match) -> str:
        body = match.group(1)
        body = re.sub(r"\s*\*/\s*$", "", body)
        body = _LEADING_MARKER.sub("", body)
        return f'R"doc({body})doc"'
    return _DOCSTRING.sub(clean, text)


def run_mkdoc(out_path: Path, include_dirs: list[str], headers: list[Path],
              extra_clang_args: list[str]) -> int:
    try:
        import pybind11_mkdoc  # noqa: F401
    except ImportError:
        fall_back(out_path, "pybind11_mkdoc not installed")
        return 0

    # Extract into a temporary file, so that a failed run cannot leave a partial header behind.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".tmp")
    cmd: list[str] = [sys.executable, "-m", "pybind11_mkdoc", "-o", str(tmp_path)]
    for inc in include_dirs:
        cmd.append(f"-I{inc}")
    cmd.extend(extra_clang_args)
    cmd.extend(str(h) for h in headers)

    print("[generate_docstrings] running:",
          " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    try:
        result = subprocess.run(cmd, check=False)
    except FileNotFoundError as e:
        fall_back(out_path, f"failed to launch pybind11_mkdoc: {e}")
        return 0

    if result.returncode != 0 or not tmp_path.exists():
        tmp_path.unlink(missing_ok=True)
        fall_back(out_path, f"pybind11_mkdoc exited with {result.returncode}")
    else:
        tmp_path.write_text(clean_docstrings(tmp_path.read_text()))
        tmp_path.replace(out_path)
        print(f"[generate_docstrings] wrote {out_path}", file=sys.stderr)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", "-o", type=Path, required=True,
                        help="Output header path (typically python/src/goddard_docstrings.h)")
    parser.add_argument("--include", "-I", action="append", default=[], metavar="DIR",
                        help="Include directory passed to libclang (repeatable)")
    parser.add_argument("--header", action="append", default=[], type=Path, metavar="HEADER",
                        help="Public C++ header to scan (repeatable)")
    parser.add_argument("--stub", action="store_true",
                        help="Always emit the stub header without invoking pybind11_mkdoc")
    parser.add_argument("--if-missing", action="store_true",
                        help="Do nothing if the output header already exists")
    parser.add_argument("--clang-arg", action="append", default=[], metavar="ARG",
                        help="Extra argument forwarded to libclang (repeatable)")
    args = parser.parse_args()

    if args.if_missing and args.output.exists():
        return 0

    if args.stub or not args.header:
        write_stub(args.output)
        return 0

    extra = list(args.clang_arg)
    if not any(a.startswith("-std=") for a in extra):
        extra.append("-std=c++20")

    return run_mkdoc(args.output, args.include, args.header, extra)


if __name__ == "__main__":
    sys.exit(main())
