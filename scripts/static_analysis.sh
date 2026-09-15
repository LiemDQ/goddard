#!/usr/bin/env bash
# Run clang-tidy and cppcheck on Goddard sources, using the compile database of a configured build.
#
# Usage:
#   scripts/static_analysis.sh [--build-dir DIR]              Analyze all project sources (report only).
#   scripts/static_analysis.sh --changed BASE [--build-dir DIR]
#       Analyze only what changed relative to the git ref BASE, and exit non-zero on any finding.
#       Both tools report only findings on changed lines (cppcheck: .cpp files in src/).
#
# Configuration: .clang-tidy (with python/src/.clang-tidy) and .cppcheck-suppressions.
# Needs only `cmake` configuration (compile_commands.json), not a compiled build.
set -euo pipefail

build_dir="build"
base_ref=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --changed) base_ref="$2"; shift 2 ;;
        --build-dir) build_dir="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [[ ! -f "$build_dir/compile_commands.json" ]]; then
    echo "No $build_dir/compile_commands.json. Configure first, e.g. 'pixi run configure'." >&2
    exit 2
fi

cppcheck_args=(
    --project="$build_dir/compile_commands.json"
    --enable=warning,style,performance,portability
    --inline-suppr
    --suppressions-list=.cppcheck-suppressions
    --std=c++20
    --template=gcc
    --quiet
    -j "$(nproc)"
)

status=0

if [[ -z "$base_ref" ]]; then
    echo "== clang-tidy (all sources)"
    run-clang-tidy -p "$build_dir" -quiet -j "$(nproc)" \
        -source-filter '.*/(src|test|python/src)/[^/]*\.cpp$' || status=1

    echo "== cppcheck (src/)"
    cppcheck "${cppcheck_args[@]}" "--file-filter=$repo_root/src/*" || status=1
    exit "$status"
fi

merge_base="$(git merge-base "$base_ref" HEAD)"
source_paths=(src include test python/src)

echo "== clang-tidy (lines changed since $base_ref)"
# Headers are not in the compile database, so changed header lines are only reported when a
# changed source file includes them.
git diff -U0 "$merge_base" -- "${source_paths[@]}" \
    | python "$CONDA_PREFIX/share/clang/clang-tidy-diff.py" \
        -p1 -path "$build_dir" -quiet -only-check-in-db -warnings-as-errors='*' \
    || status=1

echo "== cppcheck (lines changed since $base_ref, in src/)"
mapfile -t changed_sources < <(git diff --name-only --diff-filter=d "$merge_base" -- 'src/*.cpp')
if [[ ${#changed_sources[@]} -gt 0 ]]; then
    file_filters=()
    for file in "${changed_sources[@]}"; do
        file_filters+=("--file-filter=$repo_root/$file")
    done
    cppcheck_output="$(cppcheck "${cppcheck_args[@]}" "${file_filters[@]}" 2>&1 || true)"
    # Keep only findings on changed lines, so touching a file does not require fixing its old findings.
    git diff -U0 "$merge_base" -- "${changed_sources[@]}" \
        | python scripts/filter_changed_lines.py "$repo_root" <(printf '%s\n' "$cppcheck_output") \
        || status=1
else
    echo "No changed source files."
fi

exit "$status"
