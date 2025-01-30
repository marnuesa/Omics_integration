#!/bin/bash
#SBATCH --ntasks=10               # Number of task
#SBATCH --cpus-per-task=6              # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=2G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o post_fastqc_%j.log # Standard output goes to this file


################################################################################
#                       00-QC_posttrim.sh
#
#   This script execute fastqc to+ do a preliminar quality control with
#   Clean data
#
#   Author: Marta Núñez Salvador
#   Date: 01/10/2024
#   Version: 1.1 
#
################################################################################

#MODULES
module load biotools

#PATHs
path_in='/home/nuezsal/Omics_integration/RNA_seq/Libraries/01-Clean_data'
path_out='/home/nuezsal/Omics_integration/RNA_seq/Results/00-QC/QC_posttrim'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir -p "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

cd $path_in

sample_list=$(ls -d *)
for sample in $sample_list;
do

	srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastqc $sample/*.fq.gz -t 2 -o "$path_out" &
done
wait

exit 0
