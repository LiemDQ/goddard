---
name: build-and-test
description: How to configure, build and test Goddard (C++ and Python), set up worktrees, and build the docs, including sandbox and resource rules. Use before compiling, running tests, or briefing an agent that will.
---
# Building and testing Goddard

Everything runs through pixi tasks in `pixi.toml`. Use the tasks rather than raw `cmake`/`pytest` invocations: they carry the flags, environment variables and dependencies that make builds reproducible.

## Resource rule: one build or test run at a time

C++ builds (Cantera, Eigen, ASan) and test runs use a lot of memory, and the development machine has 15 GB. **Only one agent may build or run tests at any time.** Parallel agents are fine for reading and editing. When work is split across agents, they edit only (or take turns), and a single integration build and test run happens afterwards.

Build parallelism is already capped by `CMAKE_BUILD_PARALLEL_LEVEL=6` (pixi activation env). Don't pass a higher `-j`.

## Environments

| Env | Use |
|---|---|
| `default` | C++ configure, build, `ctest`, lint |
| `test` | Python tests; adds `cantera` and the `cea` package (the CEA oracle) |
| `docs` | Documentation build (`pybind11_mkdoc`, mkdocs) |

pixi sets `CCACHE_DIR` and `CPM_SOURCE_CACHE` to `.cache/` in the project, so sandboxed builds can write their caches and fresh build directories configure without network access once the CPM cache is filled.

## C++

```bash
pixi run configure-quick          # Debug, no -Werror, no static analysis: for iteration
pixi run configure                # Release with the CI flags (-Werror): run before pushing C++ changes
pixi run compile
pixi run ctest --output-on-failure     # all tests
pixi run ctest-fast               # skips tests labelled "slow" (condensed-phase suites)
./build/test/goddardTests --gtest_filter='Suite.*'   # a subset, fastest while iterating
```

- The ASan options are compiled into the test binary; no `ASAN_OPTIONS` is needed. Other binaries linked against an ASan build (e.g. `tools/moc_sweep`) need `LD_PRELOAD=$CONDA_PREFIX/lib/libasan.so.8`, which the `moc-sweep` task already sets.
- New condensed-phase test suites belong in the `GODDARD_SLOW_TEST_SUITES` list in `test/CMakeLists.txt`.

## Python

```bash
pixi run -e test test-python        # reinstalls the package, then runs all Python tests
pixi run -e test test-python-fast   # reinstalls, then skips tests marked "slow"
```

- Both tasks depend on `pip-install`, so the tests never import a stale `_core`. Run `pip-install` alone after changing bindings if you import `goddard` outside the tests.
- Mark long-running Python tests with `pytest.mark.slow` (see `test_condensed_cea.py`).
- Call the `cea` oracle in a forked child process (see `run_isolated` in `test_condensed_cea.py`): a Fortran `STOP` kills the interpreter.

## Before committing

```bash
pixi run lint-changed <base-ref>   # clang-tidy + cppcheck on changed lines; CI runs it against origin/main
```

It needs a configured `build/`, not a compiled one.

## Worktrees

```bash
pixi run setup-worktree <name> [base-ref]   # branch <name> at base (default HEAD) in .claude/worktrees/<name>
```

The script shares the main checkout's `.pixi` and `.cache` through symlinks and checks the worktree's HEAD. Brief agents with the SHA, and have them confirm `git -C <worktree> rev-parse HEAD` before starting.

Rules inside a worktree, because `.pixi` is shared with the main checkout:

- **Run every pixi command with `--as-is`** (e.g. `pixi run --as-is configure-quick`). Without it, pixi may install or update packages through the worktree path, which writes that path into the shared environment and breaks it once the worktree is deleted.
- **Never run `pip-install`, the Python test tasks, `pixi install` or `pixi reinstall` in a worktree.** `pip-install` refuses to: an editable install would repoint the shared environment at the worktree. Merge first, then run Python tests in the main checkout.
- C++ builds and `ctest` in the worktree's own `build/` are fine, subject to the resource rule.
- Remove finished worktrees with `git worktree remove .claude/worktrees/<name>` and `git branch -d <name>`.

## Docs

```bash
pixi run -e docs docs         # configure build-docs/, compile, generate stubs, mkdocs build
pixi run -e docs docs-serve   # live preview
```

The docs use their own `build-docs/` directory with the docs environment's Python, so they don't touch `build/`. Check that `python/src/goddard_docstrings.h` has real `__doc_` entries, not the stub.

## Sandbox and network

- **Needs the sandbox disabled:** `pixi install`, `pixi reinstall`, and any `pixi run` right after `pixi.toml` dependencies or `pixi.lock` change. They need pixi's global cache.
- **Allowed domains:** the first configure after the CPM cache is cleared downloads `CPM.cmake`; list `release-assets.githubusercontent.com`, `github.com` and `objects.githubusercontent.com`.
- **`$TMPDIR`** is set only inside the sandbox; don't pass files between sandboxed and unsandboxed commands through it.
- **gtest stream capture** (`testing::internal::CaptureStderr`) writes to `/tmp`, which the sandbox blocks. To check Cantera warnings, install a recording `Cantera::Logger` instead (see `recorded_warnings` in `test/test_nozzle.cpp`).

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `ccache: error: Read-only file system` | `CCACHE_DIR` isn't set: the command bypassed pixi. Run it through `pixi run` |
| Build killed, or the machine stalls | Too many parallel compile jobs or concurrent builds; see the resource rule |
| `failed to acquire global cache lock` | A pixi command needs the global cache; rerun with the sandbox disabled |
| CPM download blocked during configure | Add the allowed domains above |
| Python tests don't see a C++ change | The package wasn't reinstalled; use the test tasks, which reinstall first |
| `Error launching '<tool>': No such file or directory` | An env launcher points at a deleted worktree; run `pixi reinstall -e <env>` with the sandbox disabled |
| `could not lock config file .git/config` on branch deletion | Sandbox; harmless, the branch is still deleted |
| Docs build shows stub docstrings | `goddard_docstrings.h` was written by a non-docs build; the `docs` task removes and regenerates it |
