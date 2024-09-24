#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=1               # Number of cpus
#SBATCH -t 0-03:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=30G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o intersect_DMRs_%j.log # Standard output goes to this file

################################################################################
#
#               03-Intersect_DMRs_features.sh
#
#    Script for Processing and Annotating DMRs
# 
#    This script performs the following tasks:
#    1. Converts TSV files from the specified input path to BED format.
#    2. Iterates over BED files to annotate regions based on intersections with other BED files.
#    3. Handles unmatched regions and updates output files accordingly.
#    4. Creates a list of unique features from the processed BED files.
#    5. Merges and formats the final output into a comprehensive gene features file.
#
################################################################################


path_bed='/home/nuezsal/Omics_integration/Integration_methylome/Results/Bed_files'
path_tsv='/home/nuezsal/Omics_integration/Methylome/Results/03-DMRcaller'
path_in='/home/nuezsal/Omics_integration/Integration_methylome/Libraries/DMRs_bed'
path_out='/home/nuezsal/Omics_integration/Integration_methylome/Results/01-Gene_features_DMRs'

# Ensure path_out and in exists
mkdir -p "$path_in"
mkdir -p "$path_out"

# Transform TSV files in bed files

for file in "$path_tsv"/*duplicates_removed.tsv
do
    name_file=$(basename "$file")
    awk -F'\t' 'NR==1 {OFS="\t"; print $1, $2, $3, $16, "difference"} NR>1 {OFS="\t"; print $1, $2, $3, $16, $11 - $8}' $file | tail -n +2 > "${path_in}/${name_file%.tsv}.bed"
done

# Iterate over .bed files in path_in
for file in "$path_in"/*.bed; do
    if [ -f "$file" ]; then  # Check if $file is a regular file
        name_file=$(basename "$file")

        # Temporary file to hold unmatched regions of A
	unmatched_file="$path_out/${name_file%.bed}_undetermined.bed"
	cp $file $unmatched_file
        
        > "$path_out/${name_file%.bed}_annotation.bed" # Create or clear the output file
        
        # Loop over each B file in order of priority
	for b_file in "$path_bed"/*; do
		echo "$b_file"
		# Remove the .bed extension
		file_name_without_ext="${name_file%.bed}"

		# Split by _ and keep the last value
		name_out=$(echo "$file_name_without_ext" | rev | cut -d'_' -f1 | rev)

		# Find intersections with current B file and write to the output
		bedtools intersect -wa -wb -a $unmatched_file -b $b_file >> "$path_out/${name_file%.bed}_annotation.bed"
		    
		# Update unmatched file to exclude already matched regions
		unmatched_tmp="$path_out/unmatched_tmp.bed"
		bedtools intersect -a $unmatched_file -b $b_file -v > $unmatched_tmp
		mv $unmatched_tmp $unmatched_file
	done

	            
    else
        echo "Warning: $file is not a regular file or does not exist."
    fi
done

# UNION OF UNMATCHED AND MATCHED REGIONS
path_out_2='/home/nuezsal/Omics_integration/Integration_methylome/Results/02-Total_features_DMRs'

# Ensure path_out and in exists
mkdir -p "$path_out_2"

# Create a list of names
list=()
for file in "$path_out"/*.bed; 
do
    	# Get the file name without the path
    	file_name=$(basename "$file")
    
    	# Get the first, second, third, and fourth values separated by '_'
    	first_value=$(echo "$file_name" | cut -d'_' -f1)
    	second_value=$(echo "$file_name" | cut -d'_' -f2)
    	third_value=$(echo "$file_name" | cut -d'_' -f3)
    	fourth_value=$(echo "$file_name" | cut -d'_' -f4)
    
    	# Combine the first, second, third, and fourth values
    	combined_value="${first_value}_${second_value}_${third_value}_${fourth_value}"
    
    	# Add to the list
    	list+=("$combined_value")
done

# Remove duplicates using 'sort' and 'uniq'
unique_list=($(echo "${list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

# Iterate over the unique list
for element in "${unique_list[@]}"; 
do	
	cut -f 1,2,3,4,5,9,10 "$path_out/${element}_duplicates_removed_annotation.bed" > temp_file.bed
	awk '{print $0 "\tUR\tunknown_region"}' "$path_out/${element}_duplicates_removed_undetermined.bed" > temp_IR.bed
	cat temp_file.bed temp_IR.bed > "$path_out_2/${element}_genes.bed"	
   	rm temp_file.bed
   	rm temp_IR.bed
done

