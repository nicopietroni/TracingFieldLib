#!/bin/bash

# usage: $0 <subfolder_path> [angle]
# angle is optional: if given it is passed to quad_commandline with ANGLE_FLAG

# command line flag of quad_commandline for the angle (change it if different)
ANGLE_FLAG="-angle"

if [ -z "$1" ]; then
  echo "Usage: $0 <subfolder_path> [angle]"
  exit 1
fi

SUBFOLDER="$1"
ANGLE="$2"
DONEFOLDER="$SUBFOLDER/done"
FAILFILE="fail.txt"

# extra arguments for quad_commandline, empty if no angle is given
EXTRA_ARGS=()
if [ -n "$ANGLE" ]; then
  EXTRA_ARGS=("$ANGLE_FLAG" "$ANGLE")
  echo "Using angle: $ANGLE"
fi

mkdir -p "$DONEFOLDER"

find "$SUBFOLDER" -maxdepth 1 -type f -name "*.obj" | while read -r file; do
  echo "Processing: $file"
  ./quad_commandline "$file" "${EXTRA_ARGS[@]}"

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
