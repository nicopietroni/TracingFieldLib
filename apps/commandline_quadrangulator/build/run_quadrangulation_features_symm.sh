#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <subfolder_path>"
  exit 1
fi

SUBFOLDER="$1"
DONEFOLDER="$SUBFOLDER/done"
FAILFILE="fail.txt"

mkdir -p "$DONEFOLDER"

find "$SUBFOLDER" -maxdepth 1 -type f -name "*.obj" | while read -r file; do
  echo "Processing: $file"
  ./quad_commandline "$file" -symmetry

  STATUS=$?
  BASENAME=$(basename "$file")

  if [ $STATUS -ne 0 ]; then
    echo "❌ Failed: $BASENAME (exit code $STATUS)"
    echo "$BASENAME" >> "$FAILFILE"
  else
    mv "$file" "$DONEFOLDER/$BASENAME"
    echo "✅ Moved $BASENAME to $DONEFOLDER"
  fi
done
