#!/usr/bin/env bash
# Create a git worktree for an agent or a parallel line of work, ready to build and test.
#
# Usage: scripts/setup_worktree.sh <name> [base-ref]
#
# Creates branch <name> at <base-ref> (default HEAD) in .claude/worktrees/<name>, shares the main
# checkout's pixi environments and build caches with it, and prints the worktree path and HEAD.
# Configure and build inside the worktree as usual; it gets its own build/ directory.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 <name> [base-ref]" >&2
    exit 2
fi

name="$1"
base="${2:-HEAD}"
root="$(git rev-parse --show-toplevel)"
base_sha="$(git -C "$root" rev-parse --verify "${base}^{commit}")"
worktree="$root/.claude/worktrees/$name"

if [[ -e "$worktree" ]]; then
    echo "error: $worktree already exists" >&2
    exit 1
fi

git -C "$root" worktree add -b "$name" "$worktree" "$base_sha" >&2

# pixi environments are large and identical across worktrees; the caches are safe to share.
ln -s "$root/.pixi" "$worktree/.pixi"
mkdir -p "$root/.cache"
ln -s "$root/.cache" "$worktree/.cache"

head_sha="$(git -C "$worktree" rev-parse HEAD)"
if [[ "$head_sha" != "$base_sha" ]]; then
    echo "error: worktree HEAD $head_sha does not match base $base_sha" >&2
    exit 1
fi

echo "worktree: $worktree"
echo "HEAD:     $head_sha"
