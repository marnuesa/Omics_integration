#!/bin/bash

#SBATCH --job-name=transcripts           # Job name to show with squeue
#SBATCH --output=transcripts_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-03:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=4           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   transcripts_submit.sh
#
#   This program generate normalize count tables of genes
#
#   Author: Marta Nuñez Salvador
#   Date: 18/12/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

##################### TRANSCRIPTS ###########################            
# Paths
path_table_trans=/home/nuezsal/Omics_integration/RNA_seq/Results/02-Salmon
path_metadata_trans=/home/nuezsal/Omics_integration/RNA_seq/Additional_info
path_out_trans=/home/nuezsal/Omics_integration/Integration_MORE/Results/02-transcripts_normalize
path_annot_trans=/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript normalize_transcripts.R \
            --input $path_table_trans \
            --metadata $path_metadata_trans \
            --output $path_out_trans \
            --annotation $path_annot_trans 
            
echo -e "Gene normalization finish..."

