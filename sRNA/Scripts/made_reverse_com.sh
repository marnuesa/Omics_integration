#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=6               # Number of cpus
#SBATCH -t 0-00:10:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=10G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o pre-fastqc_%j.log # Standard output goes to this file

################################################################################
#                       made_reverse_com.sh
#
#   This script make the reverse complementary sequence of T3 sRNA because
#   this fastq are the R2 libraries
#
################################################################################

#MODULE
module load anaconda

#PATHs
path='/home/nuezsal/omics_integration_G3/sRNA_reverse/Libraries/00-Raw_data'

# Execute
# Iterate over the files in the directory
for file in "$path"/*; do
    if [ -f "$file" ]; then
        # Get the base name of the file without the path
        file_name=$(basename "$file")

        # Split the file name by '_'
        IFS='_' read -ra parts <<< "$file_name"

        # Check if the second value is 'T3'
        if [ "${parts[1]}" = "T3" ]; then
            # Execute seqkit on the file
            seqkit seq "$file" -r -p | gzip -c > "${file%_old.fastq.gz}.fastq.gz"
            echo "Processed file: $file"
        fi
    fi
done
