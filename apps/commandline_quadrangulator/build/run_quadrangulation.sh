#!/bin/bash

# Check if exactly one argument is passed
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folder>"
    exit 1
fi

INPUT_FOLDER="$1"
EXECUTABLE="./quad_commandline"  # Modify this path if needed

# Loop through all .obj files in the folder
for OBJ_FILE in "$INPUT_FOLDER"/*.obj; do
    # Skip if no .obj files are found
    [ -e "$OBJ_FILE" ] || continue

    BASE_NAME="${OBJ_FILE%.obj}"  # Remove .obj extension

    # Start building the command
    CMD="$EXECUTABLE \"$OBJ_FILE\""

    # Add .field if it exists
    if [ -f "$BASE_NAME.field" ]; then
        CMD="$CMD \"$BASE_NAME.field\""
    fi

    # Add .ffield if it exists
    if [ -f "$BASE_NAME.ffield" ]; then
        CMD="$CMD \"$BASE_NAME.ffield\""
    fi

    # Add .features if it exists
    if [ -f "$BASE_NAME.features" ]; then
        CMD="$CMD \"$BASE_NAME.features\""
    fi

    # Add .feature if it exists
    if [ -f "$BASE_NAME.feature" ]; then
        CMD="$CMD \"$BASE_NAME.feature\""
    fi

    echo "Running: $CMD"
    eval $CMD
done
