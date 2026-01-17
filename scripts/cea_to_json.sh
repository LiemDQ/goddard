#!/bin/bash

# Shell script to run cea_parser.py on CEA output files
# Usage: cea_to_json.sh [--exclude|-E FILE...] [--list|-L FILE...]

set -e

# Default values
EXCLUDE_FILES=()
LIST_FILES=()
USE_LIST=false
USE_PRETTY=false
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../data/cea_results"

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "OPTIONS:"
    echo "  --exclude, -E FILE...    Skip specific .output files"
    echo "  --list, -L FILE...       Run only on listed .output files"
    echo "  --pretty, -P             Pass --pretty flag to cea_parser.py"
    echo "  --help, -h               Show this help message"
    echo ""
    echo "If both --exclude and --list are used, --list takes precedence."
    echo "By default, runs on all .output files in data/cea_results/"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --exclude|-E)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                EXCLUDE_FILES+=("$1")
                shift
            done
            ;;
        --list|-L)
            shift
            USE_LIST=true
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                LIST_FILES+=("$1")
                shift
            done
            ;;
        --pretty|-P)
            USE_PRETTY=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Check if cea_parser.py exists
if [[ ! -f "$SCRIPT_DIR/cea_parser.py" ]]; then
    echo "Error: cea_parser.py not found in $SCRIPT_DIR"
    exit 1
fi

# Check if data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: CEA results directory not found: $DATA_DIR"
    exit 1
fi

# Function to check if file should be excluded
should_exclude() {
    local file="$1"
    for exclude_file in "${EXCLUDE_FILES[@]}"; do
        if [[ "$(basename "$file")" == "$exclude_file" ]]; then
            return 0
        fi
    done
    return 1
}

# Get list of files to process
if [[ "$USE_LIST" == true ]]; then
    # Use only listed files
    FILES_TO_PROCESS=()
    for file in "${LIST_FILES[@]}"; do
        # Add .output extension if not present
        if [[ ! "$file" =~ \.output$ ]]; then
            file="${file}.output"
        fi
        
        full_path="$DATA_DIR/$file"
        if [[ -f "$full_path" ]]; then
            FILES_TO_PROCESS+=("$full_path")
        else
            echo "Warning: File not found: $full_path"
        fi
    done
else
    # Use all .output files, excluding specified ones
    FILES_TO_PROCESS=()
    for file in "$DATA_DIR"/*.output; do
        if [[ -f "$file" ]]; then
            if ! should_exclude "$file"; then
                FILES_TO_PROCESS+=("$file")
            else
                echo "Excluding: $(basename "$file")"
            fi
        fi
    done
fi

# Check if we have any files to process
if [[ ${#FILES_TO_PROCESS[@]} -eq 0 ]]; then
    echo "No .output files found to process"
    exit 0
fi

# Process each file
echo "Processing ${#FILES_TO_PROCESS[@]} file(s)..."
for file in "${FILES_TO_PROCESS[@]}"; do
    # echo "Processing: $(basename "$file")"
    if [[ "$USE_PRETTY" == true ]]; then
        python "$SCRIPT_DIR/cea_parser.py" --pretty "$file"
    else
        python "$SCRIPT_DIR/cea_parser.py" "$file"
    fi
done

echo "Done!"