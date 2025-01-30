#!/bin/bash
#SBATCH --ntasks=23               # Number of task
#SBATCH --cpus-per-task=2               # Number of cpus
#SBATCH -t 0-07:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=5G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o Trimming_%j.log # Standard output goes to this file

################################################################################
#
#               01-Trimming.sh
#   This program trims the sequences from a
#   specific library using the fastp program.
#   
#   It is important to remember that T3 samples
#   was transformed from reverse to forward with seqkit
#
#   Arguments:
#       Fastq.gz file path
#       Output directory path
#       Path to the file with Illumina adapters.
#
#   Author: Marta Núñez Salvador
#   Date: 19/09/2024
#   Version: 1.1 
################################################################################

# MODULE
module load biotools

# PATHs
path_in='/home/nuezsal/Omics_integration/sRNA/Libraries/00-Raw_data'
path_out='/home/nuezsal/Omics_integration/sRNA/Libraries/01-Clean_data'
path_adapters='/home/nuezsal/Omics_integration/sRNA/Additional_info/all_adapters.fa'

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
    # Trimming
    srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive fastp --adapter_fasta $path_adapters \
        -i ${path_in}/$sample \
        -o ${path_out}/$sample \
        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
        --length_required 20 --trim_poly_x --poly_x_min_len 10 \
        --length_limit 25 --n_base_limit 0 & # --n_base_limit is for avoiding Ns

done
wait
exit 0
