#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folder_name>"
    exit 1
fi

# Get folder name
folder_name="$1"

# Check if the folder exists
if [ ! -d "$folder_name" ]; then
    echo "Folder not found: $folder_name"
    exit 1
fi

# Find and process pairs of .obj and .ffield files in the folder
for obj_file in "$folder_name"/*.obj; do
    if [ -f "$obj_file" ]; then
        # Get the corresponding .ffield file
            ./quad_commandline "$obj_file" dual smoothtangent no_match_sing_val
    fi
done
