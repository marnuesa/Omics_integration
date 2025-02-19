#!/bin/bash

#SBATCH --job-name=diffea           # Job name to show with squeue
#SBATCH --output=diffea_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-03:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=2           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   DEA_submit.sh
#
#   This program executes the DEA.r program to perform differential 
#   expression analysis using the absolute count tables.
#
#   Author: Marta Nuñez Salvador
#   Date: 18/07/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_table=/home/nuezsal/Omics_integration/sRNA/Results/03-Fusion_count_tables_RF
path_metadata=/home/nuezsal/Omics_integration/sRNA/Additional_info/metadata_sRNA.tsv
path_dea=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA
path_graph=/home/nuezsal/Omics_integration/sRNA/Results/DESeq_graphs
alpha=0.05


# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript DEA.R \
            --input $path_table \
            --metadata $path_metadata \
            --output $path_dea \
            --alpha $alpha \
            --specie "cume" \
            --project "Omics_project" \
            --graphs $path_graph
