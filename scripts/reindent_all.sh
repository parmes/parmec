#!/usr/bin/env bash
# Reindent all source files recursively using Vim's smartindent rules.
# Affects: .cpp, .ispc, .h, .py

set -euo pipefail

# Optional: preview only (dry-run)
# Set to 1 to see which files will be processed without modifying them
DRY_RUN=0

# Root directory (default: current)
ROOT_DIR="${1:-.}"

echo "Reindenting files under: $ROOT_DIR"
echo

find "$ROOT_DIR" -type f \( \
  -name "*.cpp" -o -name "*.ispc" -o -name "*.h" \
\) | while IFS= read -r file; do
  echo "→ $file"
  if [ "$DRY_RUN" -eq 0 ]; then
    # Run Vim silently with user's .vimrc and auto-indent the file
    vim -u ~/.vimrc -E -s "$file" \
      -c 'normal! gg=G' \
      -c 'wq'
  fi
done

echo
echo "✅ Reindent complete."

