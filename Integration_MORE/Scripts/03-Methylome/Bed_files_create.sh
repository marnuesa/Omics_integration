#!/bin/bash
################################################################################
#                                                                              #
#                              Bed_files_create.sh                             #
#                                                                              #
#   Description:                                                               #
#   ------------                                                               #
#   This script generates BED files from a GFF3 annotation file, specifically  #
#   extracting upstream regulatory regions and gene body regions for a given   #
#   set of genes. It processes gene coordinates based on strand orientation to #
#   accurately define upstream regions and gene bodies.                        #
#                                                                              #
#   Key Features:                                                              #
#   -------------                                                              #
#   1. **Input Files**:                                                        #
#      - A gene list file (`tx2gene.txt`) containing gene identifiers.         #
#      - A GFF3 annotation file with genomic features.                         #
#                                                                              #
#   2. **Region Extraction**:                                                  #
#      - For genes on the '+' strand: Extracts 2kb upstream of the start codon.#
#      - For genes on the '-' strand: Extracts 2kb downstream of the end codon.#
#      - Correctly adjusts gene boundaries based on 5' and 3' UTR annotations. #
#                                                                              #
#   3. **Output**:                                                             #
#      - BED file for upstream regions (`CMelon_DHL92_v4_upstream.bed`).       #
#      - BED file for gene body regions (`CMelon_DHL92_v4_genes.bed`).         #
#                                                                              #
#   4. **Usage**:                                                              #
#      The script allows customization of input and output paths through       #
#      command-line arguments. It ensures proper directory creation and        #
#      cleans previous outputs before generating new ones.                     #
#                                                                              #
#   Arguments:                                                                 #
#   -----------                                                                #
#   -g  Path to the gene list file (`tx2gene.txt`).                            #
#   -f  Path to the GFF3 annotation file.                                      #
#   -o  Output directory for the generated BED files.                          #
#   -h  Display help and usage information.                                    #
#                                                                              #
#   Example Command:                                                           #
#   ----------------                                                           #
#   bash Bed_files_create.sh -g data/tx2gene.txt -f annotations/genes.gff3 \   #
#                             -o output/BED_files                              #
#                                                                              #
#   Author: Marta Núñez Salvador                                              #
#   Version: 1.0                                                              #
#   Date: 07/01/2025                                                          #
#                                                                              #
################################################################################

# Function to display help message
usage() {
  echo "Usage: $0 -g <gene_file> -f <gff_file> -o <output_directory>"
  echo
  echo "  -g  Path to the tx2gene.txt file"
  echo "  -f  Path to the GFF3 file"
  echo "  -o  Output directory for the BED files"
  echo "  -h  Show this help message"
  exit 1
}

# Process command-line arguments
while getopts ":g:f:o:h" opt; do
  case ${opt} in
    g ) gene_file=$OPTARG ;;   # Path to the tx2gene.txt file
    f ) gff_file=$OPTARG ;;    # Path to the GFF3 file
    o ) bed_path=$OPTARG ;;    # Output directory path
    h ) usage ;;               # Show help message
    \? ) echo "Invalid option: -$OPTARG" >&2; usage ;;  # Invalid option handling
    : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;  # Missing argument handling
  esac
done

# Check if the required arguments were provided
if [[ -z "$gene_file" || -z "$gff_file" || -z "$bed_path" ]]; then
  echo "Error: Missing required arguments." >&2
  usage
fi

# Define output file paths for upstream and gene regions
output_up="$bed_path/CMelon_DHL92_v4_upstream.bed"
output_genes="$bed_path/CMelon_DHL92_v4_genes.bed"

# Create output directory if it doesn't exist
mkdir -p "$bed_path"

# Clear the output files if they already exist
> "$output_up"
> "$output_genes"

# Process the GFF3 file
echo "...Processing GFF3 file..."

# Iterate through each gene ID in the gene file
cut -f 2 "$gene_file" | while read -r gene_id; do
  # Extract lines from the GFF3 file matching the gene ID
  grep "$gene_id" "$gff_file" > "temp_gene.gff3"
  lines=$(wc -l < "temp_gene.gff3")  # Count the number of matching lines

  # Parse the GFF3 file for relevant columns
  awk -F '\t' '{print NR, $1, $3, $4, $5, $7, $9}' "temp_gene.gff3" | while read line_num chr types start end strand gene_name; do
    if [ "$strand" == "+" ]; then  # If the gene is on the forward strand
      three_prime_UTR_processed=false

      if [ "$types" == "gene" ]; then
        # Define gene boundaries
        start_gene=$start
        end_gene=$end
        # Define upstream region (2kb before the start codon)
        start_upstream=$(($start > 2000 ? $start - 2000 : 1))
        end_upstream=$(($start - 1))
      elif [ "$types" == "five_prime_UTR" ]; then
        # Adjust upstream region if a 5' UTR exists
        end_upstream=$end
        start_gene=$(($end_upstream + 1))
      elif [ "$types" == "three_prime_UTR" ] && [ "$three_prime_UTR_processed" == false ]; then
        # Adjust gene end if a 3' UTR exists
        three_prime_UTR_processed=true
        start_downstream=$start
        end_gene=$(($start_downstream - 1))
      fi
    else  # If the gene is on the reverse strand
      five_prime_UTR_processed=false

      if [ "$types" == "gene" ]; then
        # Define gene boundaries
        start_gene=$start
        end_gene=$end
        # Define upstream region (2kb after the end codon for reverse strand)
        start_upstream=$(($end + 1))
        end_upstream=$(($end + 2000))
      elif [ "$types" == "five_prime_UTR" ] && [ "$five_prime_UTR_processed" == false ]; then
        # Adjust upstream region if a 5' UTR exists
        five_prime_UTR_processed=true
        start_upstream=$start
        end_gene=$(($start_upstream - 1))
      elif [ "$types" == "three_prime_UTR" ]; then
        # Adjust gene start if a 3' UTR exists
        end_downstream=$end
        start_gene=$(($end_downstream + 1))
      fi
    fi

    # Write the regions to BED files when processing the last line
    if [ "$line_num" -eq "$lines" ]; then
      # Write upstream region to the upstream BED file (0-based indexing)
      echo -e "$chr\t$(($start_upstream - 1))\t$(($end_upstream - 1))\t$strand\t$gene_id" >> "$output_up"
      # Write gene region to the gene BED file (0-based indexing)
      echo -e "$chr\t$(($start_gene - 1))\t$(($end_gene - 1))\t$strand\t$gene_id" >> "$output_genes"
    fi
  done
done

echo "...Feature extraction completed..."

# Remove temporary file
rm "temp_gene.gff3"

cat /home/nuezsal/Omics_integration/Integration_methylome/Results/02-Total_features_DMRs/*.bed > "$bed_path/merged_files.bed"
bedtools intersect -a $output_up -b "$bed_path/merged_files.bed" -wa > "$bed_path/CMelon_DHL92_v4_upstream_filter.bed"
bedtools intersect -a $output_genes -b "$bed_path/merged_files.bed" -wa > "$bed_path/CMelon_DHL92_v4_genes_filter.bed"

sort "$bed_path/CMelon_DHL92_v4_upstream_filter.bed" | uniq > "$bed_path/CMelon_DHL92_v4_upstream_filter_uniq.bed"
sort "$bed_path/CMelon_DHL92_v4_genes_filter.bed" | uniq > "$bed_path/CMelon_DHL92_v4_genes_filter_uniq.bed"
