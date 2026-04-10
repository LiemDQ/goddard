#!/usr/bin/env python3
"""Generate ``goddard_docstrings.h`` from Goddard's public C++ headers.

Thin wrapper around `pybind11_mkdoc <https://github.com/pybind/pybind11_mkdoc>`_.
It parses the requested headers with libclang, extracts Doxygen-style comments,
and emits a header containing ``__doc_*`` string constants plus a ``DOC(...)``
macro that nanobind binding files can use::

    .def("set_state_TP", &Goddard::Gas::set_state_TP,
         "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))

If ``pybind11_mkdoc`` (or its libclang dependency) is not available, OR if
``--stub`` is passed, the script writes a stub header that defines
``DOC(...)`` as an empty string. This keeps the bindings compiling cleanly in
environments that don't have the docs toolchain installed (e.g. plain
``pixi run compile``).
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

STUB_CONTENTS = """\
#pragma once
// Stub goddard_docstrings.h emitted by scripts/generate_docstrings.py.
//
// Generated when pybind11_mkdoc is unavailable or when -Dgoddard_BUILD_DOCSTRINGS=OFF.
// DOC(...) expands to an empty string so binding sources compile cleanly without
// embedded docstrings. To populate real docstrings, configure with
// -Dgoddard_BUILD_DOCSTRINGS=ON (the `docs` pixi environment does this automatically).

#define DOC(...) ""
"""


def write_stub(out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(STUB_CONTENTS)
    print(f"[generate_docstrings] wrote stub -> {out_path}", file=sys.stderr)


def run_mkdoc(out_path: Path, include_dirs: list[str], headers: list[Path],
              extra_clang_args: list[str]) -> int:
    try:
        import pybind11_mkdoc  # noqa: F401
    except ImportError:
        print("[generate_docstrings] pybind11_mkdoc not installed; falling back to stub",
              file=sys.stderr)
        write_stub(out_path)
        return 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [sys.executable, "-m", "pybind11_mkdoc", "-o", str(out_path)]
    for inc in include_dirs:
        cmd.append(f"-I{inc}")
    cmd.extend(extra_clang_args)
    cmd.extend(str(h) for h in headers)

    print("[generate_docstrings] running:",
          " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    try:
        result = subprocess.run(cmd, check=False)
    except FileNotFoundError as e:
        print(f"[generate_docstrings] failed to launch pybind11_mkdoc: {e}", file=sys.stderr)
        write_stub(out_path)
        return 0

    if result.returncode != 0:
        print(f"[generate_docstrings] pybind11_mkdoc exited with {result.returncode}; "
              "falling back to stub", file=sys.stderr)
        write_stub(out_path)
    else:
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
    parser.add_argument("--clang-arg", action="append", default=[], metavar="ARG",
                        help="Extra argument forwarded to libclang (repeatable)")
    args = parser.parse_args()

    if args.stub or not args.header:
        write_stub(args.output)
        return 0

    extra = list(args.clang_arg)
    if not any(a.startswith("-std=") for a in extra):
        extra.append("-std=c++20")

    return run_mkdoc(args.output, args.include, args.header, extra)


if __name__ == "__main__":
    sys.exit(main())
