"""Filter compiler-style diagnostics to lines changed in a unified diff.

Usage: git diff -U0 BASE -- FILES | python scripts/filter_changed_lines.py REPO_ROOT DIAGNOSTICS_FILE

Diagnostics are lines of the form `path:line:column: severity: message`. Each diagnostic whose
location is a changed line is printed with the lines that follow it (source excerpt), and the exit
status is 1 if any were printed.
"""
import os
import re
import sys

HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
DIAGNOSTIC = re.compile(r"^(.+?):(\d+):\d+: \w+: ")


def changed_lines(diff_text):
    changed = {}
    current = None
    for line in diff_text.splitlines():
        if line.startswith("+++ "):
            path = line[4:]
            current = path[2:] if path.startswith("b/") else None
            continue
        match = HUNK.match(line)
        if match and current:
            start = int(match.group(1))
            count = int(match.group(2)) if match.group(2) is not None else 1
            changed.setdefault(current, set()).update(range(start, start + count))
    return changed


def main():
    repo_root, diagnostics_path = sys.argv[1], sys.argv[2]
    changed = changed_lines(sys.stdin.read())
    with open(diagnostics_path) as handle:
        diagnostics = handle.read().splitlines()

    found = False
    keep = False
    for line in diagnostics:
        match = DIAGNOSTIC.match(line)
        if match:
            path = os.path.relpath(match.group(1), repo_root)
            keep = int(match.group(2)) in changed.get(path, set())
            found = found or keep
        if keep:
            print(line)
    sys.exit(1 if found else 0)


if __name__ == "__main__":
    main()
