#!/bin/bash

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 environment_file"
    exit 1
fi

ENVIRONMENT_FILE=$1

# Function to process each directory
process_directory() {
    local DIR=$1

    if [ -f "$DIR/Makefile" ]; then
        echo "Processing $DIR..."
        (cd "$DIR" && ./mk $2)
    fi
}

export -f process_directory

# Find all directories and process them
find . -type d -print0 | while IFS= read -r -d '' DIR; do
    process_directory "$DIR" "$ENVIRONMENT_FILE"
done

echo "Done."

exit 0
