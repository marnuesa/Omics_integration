#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=6               # Number of cpus
#SBATCH -t 0-00:10:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=10G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o pre-fastqc_%j.log # Standard output goes to this file

################################################################################
#                       00-PreFastQC.sh
#
#   This script execute fastqc to do a preliminar quality control with
#   Raw data
#
#   Author: Marta Núñez Salvador
#   Date: 19/09/2024
#   Version: 1.1 
################################################################################

#MODULE
module load biotools

#PATHs
path_in='/home/nuezsal/Omics_integration/sRNA/Libraries/00-Raw_data'
path_out='/home/nuezsal/Omics_integration/sRNA/Results/00-PreFastQC'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

# Execute
srun -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastqc ${path_in}/*fastq.gz -t 6 -o ${path_out}
