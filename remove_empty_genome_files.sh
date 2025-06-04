#!/bin/bash

# Usage: ./clean_empty_files.sh /path/to/folder list.txt

FOLDER="${1:-.}"
LIST_FILE="$2"

# Sanity check
if [ ! -d "$FOLDER" ]; then
    echo "Error: $FOLDER is not a valid directory"
    exit 1
fi

if [ ! -f "$LIST_FILE" ]; then
    echo "Error: $LIST_FILE is not a valid file"
    exit 1
fi

# Find empty files and store their names
empty_files=$(find "$FOLDER" -maxdepth 1 -type f -empty)

# Delete each empty file and remove matching lines from the list
while IFS= read -r file; do
    echo "Removing: $file"
    rm -f "$file"
    
    # Escape file path for sed
    escaped_file=$(printf '%s\n' "$file" | sed 's/[.[\*^$/]/\\&/g')
    
    # Remove matching lines from the list file
    sed -i "\|$escaped_file|d" "$LIST_FILE"
done <<< "$empty_files"