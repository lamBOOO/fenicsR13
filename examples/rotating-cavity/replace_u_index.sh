#!/usr/bin/env bash
# Replace occurrences of a token (default: u_0) with another (default: u_1)
# in-place for a given file, with macOS/BSD and GNU sed compatibility.
#
# Usage:
#   bash examples/rotating-cavity/replace_u_index.sh \
#     /path/to/file.xdmf [FROM] [TO] [NEW_TIME]
# Example:
#   bash examples/rotating-cavity/replace_u_index.sh \
#     examples/rotating-cavity/collected_files/0/u_0.xdmf u_0 u_1 3.5

set -euo pipefail

if [ "${1:-}" = "" ]; then
  echo "Usage: $0 /path/to/file [FROM] [TO] [NEW_TIME]" >&2
  exit 1
fi

FILE="$1"
FROM_TOKEN="${2:-u_0}"
TO_TOKEN="${3:-u_1}"
NEW_TIME="${4:-}"

if [ ! -f "$FILE" ]; then
  echo "ERROR: File not found: $FILE" >&2
  exit 2
fi

replace_inplace() {
  local pattern="$1"
  local replacement="$2"
  # Detect GNU sed vs BSD sed (macOS)
  if sed --version >/dev/null 2>&1; then
    sed -i -e "s|${pattern}|${replacement}|g" "$FILE"
  else
    sed -i '' -e "s|${pattern}|${replacement}|g" "$FILE"
  fi
}

# Replace token FROM_TOKEN -> TO_TOKEN
# Escape '|' in tokens to avoid delimiter collision
SAFE_FROM=${FROM_TOKEN//|/\|}
SAFE_TO=${TO_TOKEN//|/\|}
replace_inplace "$SAFE_FROM" "$SAFE_TO"

# Optionally replace <Time Value="0" /> with provided NEW_TIME
if [ -n "$NEW_TIME" ]; then
  # Escape '|' in NEW_TIME
  SAFE_TIME=${NEW_TIME//|/\|}
  replace_inplace '<Time Value="0" />' "<Time Value=\"${SAFE_TIME}\" />"
fi

if [ -n "$NEW_TIME" ]; then
  echo "Replaced '${FROM_TOKEN}' -> '${TO_TOKEN}' and Time=0 -> ${NEW_TIME} in: $FILE"
else
  echo "Replaced '${FROM_TOKEN}' -> '${TO_TOKEN}' in: $FILE"
fi
