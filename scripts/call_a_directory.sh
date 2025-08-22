#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# call_a_directory.sh – batch process multiple samples with testCallPerturb.sh
#
# Usage:
#   ./call_a_directory.sh BASE_DIR METHOD
#
# Arguments:
#   BASE_DIR  - Directory containing subdirectories (one per sample)
#   METHOD    - Method to pass to testCallPerturb.sh
#               (flex|em|em-fixed|em-noncells|simple|otsu|quantile|model3)
#
# Example:
#   ./call_a_directory.sh /storage/gene_features em-fixed
#
# This will:
#   1. Find all subdirectories in /storage/gene_features
#   2. For each subdirectory (e.g., "sample1"), call:
#      testCallPerturb.sh em-fixed /storage/gene_features/sample1 sample1
# ---------------------------------------------------------------------------

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 BASE_DIR METHOD"
    echo ""
    echo "Arguments:"
    echo "  BASE_DIR  - Directory containing sample subdirectories"
    echo "  METHOD    - Demultiplexing method to use"
    echo ""
    echo "Valid methods: flex, em, em-fixed, em-noncells, simple, otsu, quantile, model3"
    echo ""
    echo "Example:"
    echo "  $0 /storage/gene_features em-fixed"
    echo ""
    exit 1
fi

BASE_DIR="$1"
METHOD="$2"

# Validate that base directory exists
if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: Base directory '$BASE_DIR' does not exist"
    exit 1
fi

# Validate method (basic check - testCallPerturb.sh will do full validation)
if [[ ! "$METHOD" =~ ^(flex|em|em-fixed|em-noncells|simple|otsu|quantile|model3)$ ]]; then
    echo "ERROR: Invalid method '$METHOD'"
    echo "Valid methods: flex, em, em-fixed, em-noncells, simple, otsu, quantile, model3"
    exit 1
fi

# Get script directory to find testCallPerturb.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_SCRIPT="$SCRIPT_DIR/testCallPerturb.sh"

if [[ ! -f "$TEST_SCRIPT" ]]; then
    echo "ERROR: Cannot find testCallPerturb.sh at '$TEST_SCRIPT'"
    exit 1
fi

# Find all subdirectories in BASE_DIR
echo "Scanning for subdirectories in: $BASE_DIR"
echo "Method: $METHOD"
echo "----------------------------------------"

# Counter for tracking progress
total_dirs=0
processed_dirs=0
failed_dirs=0

# First, count total directories
for subdir in "$BASE_DIR"/*; do
    if [[ -d "$subdir" ]]; then
        ((total_dirs++))
    fi
done

if [[ $total_dirs -eq 0 ]]; then
    echo "No subdirectories found in '$BASE_DIR'"
    exit 0
fi

echo "Found $total_dirs subdirectories to process"
echo ""

# Process each subdirectory
for subdir in "$BASE_DIR"/*; do
    if [[ -d "$subdir" ]]; then
        # Extract just the subdirectory name (not full path)
        subdir_name=$(basename "$subdir")
        
        echo "[$((processed_dirs + 1))/$total_dirs] Processing: $subdir_name"
        echo "  MTX_DIR: $subdir"
        echo "  SC_ID: $subdir_name"
        
        # Call testCallPerturb.sh with the three arguments
        if "$TEST_SCRIPT" "$METHOD" "$BASE_DIR" "$subdir_name"; then
            echo "  ✓ Success"
            ((processed_dirs++))
        else
            echo "  ✗ Failed (exit code: $?)"
            ((failed_dirs++))
        fi
        
        echo ""
    fi
done

# Summary
echo "========================================="
echo "Batch processing complete"
echo "Total directories: $total_dirs"
echo "Successfully processed: $processed_dirs"
echo "Failed: $failed_dirs"

if [[ $failed_dirs -gt 0 ]]; then
    echo ""
    echo "WARNING: $failed_dirs directories failed processing"
    exit 1
else
    echo ""
    echo "All directories processed successfully!"
    exit 0
fi
