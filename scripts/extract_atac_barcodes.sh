#!/bin/bash

# Extract barcodes from atac_fragments.tsv.gz files in subdirectories
# and save them to a single CSV file with directory-based prefixes

# Set the base directory
BASE_DIR="/cluster/projects/wouterslab/ArchR103_4/data/human_multiome"
OUTPUT_FILE="combined_barcodes.csv"
TEMP_DIR=$(mktemp -d)

# Function to clean up temporary files
cleanup() {
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

# Check if base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory $BASE_DIR does not exist."
    echo "Please update the BASE_DIR variable with the correct path."
    exit 1
fi

echo "Extracting barcodes from: $BASE_DIR"
echo "========================================================"

# Create CSV header
echo "original_directory,original_barcode,prefixed_barcode" > "$OUTPUT_FILE"

# Counter for total barcodes
total_barcodes=0

# Process each subdirectory
for subdir in "$BASE_DIR"/*; do
    if [ -d "$subdir" ]; then
        fragments_file="$subdir/atac_fragments.tsv.gz"
        
        if [ -f "$fragments_file" ]; then
            # Get directory name for prefix
            prefix=$(basename "$subdir")
            echo "Processing $fragments_file"
            
            # Extract unique barcodes from the 4th column
            # Skip header lines starting with #, extract 4th column, sort and get unique values
            temp_barcodes="$TEMP_DIR/${prefix}_barcodes.txt"
            
            echo "  Extracting barcodes..."
            zcat "$fragments_file" | \
            awk '!/^#/ && NF >= 4 {print $4}' | \
            sort -u > "$temp_barcodes"
            
            # Count barcodes for this directory
            barcode_count=$(wc -l < "$temp_barcodes")
            echo "  Found $barcode_count unique barcodes in $prefix"
            
            # Add prefixed barcodes to output file
            echo "  Adding prefixed barcodes to output..."
            while IFS= read -r barcode; do
                echo "$prefix,$barcode,$prefix#$barcode"
            done < "$temp_barcodes" >> "$OUTPUT_FILE"
            
            # Update total count
            total_barcodes=$((total_barcodes + barcode_count))
            
            # Clean up temp file
            rm "$temp_barcodes"
            
        else
            echo "No atac_fragments.tsv.gz found in $subdir"
        fi
    fi
done

echo "========================================================"
echo "Process completed successfully!"
echo "Total barcodes processed: $total_barcodes"
echo "Output file: $OUTPUT_FILE"
echo ""

# Display summary statistics
echo "Summary by directory:"
if [ -f "$OUTPUT_FILE" ]; then
    # Skip header line and count occurrences of each directory
    tail -n +2 "$OUTPUT_FILE" | cut -d',' -f1 | sort | uniq -c | \
    while read count dir; do
        printf "  %-25s: %'d barcodes\n" "$dir" "$count"
    done
fi

echo ""
echo "First few lines of output:"
head -5 "$OUTPUT_FILE"
echo ""
echo "Output file saved as: $OUTPUT_FILE"