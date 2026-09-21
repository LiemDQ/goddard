#!/usr/bin/env bash
# Editable install of the goddard Python package into the active pixi environment.
#
# Worktrees made by scripts/setup_worktree.sh share the main checkout's .pixi, so an editable install
# from a worktree would repoint the shared environment at the worktree's sources (and write its path
# into the environment). Install and run Python tests from the main checkout instead.
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
if [[ -L "$root/.pixi" ]]; then
    echo "error: $root is a worktree sharing the main checkout's pixi environments;" \
         "run pip-install and the Python tests from the main checkout." >&2
    exit 1
fi

pip install -e "$root" --no-build-isolation
