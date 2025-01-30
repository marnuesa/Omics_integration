#!/bin/bash
#SBATCH --ntasks=1              # Number of task
#SBATCH --cpus-per-task=4               # Number of cpus
#SBATCH -t 0-00:20:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=5G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o fastqtofasta_%j.log # Standard output goes to this file

################################################################################
#
#               03-fastqtofasta.sh
#   This program create fasta from fastq.gz because the next 
#   program needs fasta files
#
#  Author: Marta Núñez Salvador
#   Date: 20/09/2024
#   Version: 1.1 
################################################################################
# PATHs

path_in='/home/nuezsal/Omics_integration/sRNA/Libraries/01-Clean_data'
path_out='/home/nuezsal/Omics_integration/sRNA/Libraries/02-Clean_data_fasta'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

# List species directories
sample_list=$( ls ${path_in}/ )

# Iterate files list
for sample in $sample_list
do
    # fastqtofasta
    first_part=$(echo $sample | cut -d'.' -f1)
    zcat "${path_in}/$sample" | sed -n '1~4s/^@/>/p;2~4p' > ${path_out}/${first_part}.fasta
done
wait
exit 0
