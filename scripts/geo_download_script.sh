#!/bin/bash
# ==============================================================================
# Purpose: Downloads 10x Multiome supplementary files from NCBI GEO FTP for 26
#          glioblastoma samples (GSM8820702 - GSM8820727; Wu et al.).
# Downloads: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz, atac_fragments.tsv.gz
# Output:  Saved into sample-specific folders under 'geo_downloads/'
# Next steps: Processed by scripts/convert_10x_geo_files_to_h5.sh for ArchR/Seurat
# ==============================================================================

# List of samples (extracted from atac_fragments filenames)
samples=(
    "GSM8820702_6419_A"
    "GSM8820703_6419_B"
    "GSM8820704_6419_C"
    "GSM8820705_6419_D"
    "GSM8820706_6425_A"
    "GSM8820707_6425_C"
    "GSM8820708_6425_D"
    "GSM8820709_6425_E"
    "GSM8820710_6425_F"
    "GSM8820711_6425_G"
    "GSM8820712_6434_A"
    "GSM8820713_6467_A"
    "GSM8820714_6467_B"
    "GSM8820715_6467_C"
    "GSM8820716_6467_D"
    "GSM8820717_6467_E"
    "GSM8820718_6467_F"
    "GSM8820719_6509_A"
    "GSM8820720_6509_B"
    "GSM8820721_6509_C"
    "GSM8820722_6509_D"
    "GSM8820723_6509_E"
    "GSM8820724_6514_A"
    "GSM8820725_6514_B"
    "GSM8820726_6514_C"
    "GSM8820727_6514_D"
)

# Base URL pattern
base_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8820nnn"

# Files to download for each sample
file_types=("barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz" "atac_fragments.tsv.gz")

# Create output directory if it doesn't exist
output_dir="geo_downloads"
mkdir -p "$output_dir"

# Counter for progress
total_samples=${#samples[@]}
current=0

echo "Starting download of files for $total_samples samples..."
echo "Output directory: $output_dir"
echo ""

# Loop through each sample
for sample in "${samples[@]}"; do
    current=$((current + 1))
    
    # Extract GSM ID (e.g., GSM8820702)
    gsm_id=$(echo "$sample" | cut -d'_' -f1)
    
    echo "[$current/$total_samples] Processing $sample (${gsm_id})..."
    
    # Create sample-specific directory
    sample_dir="$output_dir/$sample"
    mkdir -p "$sample_dir"
    
    # Download each file type
    for file_type in "${file_types[@]}"; do
        # Output filename (without URL encoding)
        output_file="$sample_dir/${sample}_${file_type}"
        
        # Check if the file already exists
        if [ -f "$output_file" ]; then
            echo "  - Already exists: $file_type"
            continue
        fi
        
        # Construct the filename with URL encoding (%5F for underscore)
        encoded_filename="${sample//_/%5F}_${file_type}"
        
        # Construct the full URL
        url="$base_url/$gsm_id/suppl/$encoded_filename"
        
        # Download the file
        echo "  Downloading $file_type..."
        if command -v wget >/dev/null 2>&1; then
            wget -q -O "$output_file" "$url"
        else
            curl -s -L -o "$output_file" "$url"
        fi
        
        if [ $? -eq 0 ]; then
            echo "    ✓ Success: $file_type"
        else
            echo "    ✗ Failed: $file_type"
        fi
    done
    
    echo ""
done

echo "Download complete!"
echo "Files saved to: $output_dir"