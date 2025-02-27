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

################################# With TH ######################################
# Paths
path_in=/home/nuezsal/Omics_integration/RNA_seq/Results/02-Salmon
path_metadata=/home/nuezsal/Omics_integration/RNA_seq/Additional_info
path_dea=/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH
path_graph=/home/nuezsal/Omics_integration/RNA_seq/Results/DESeq_graphs_TH
alpha=0.05
TH=0.585


# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript DEA_transcripts.R \
            --input $path_in \
            --metadata $path_metadata \
            --output $path_dea \
            --alpha $alpha \
            --specie "cume" \
            --project "Omics_project" \
            --graphs $path_graph \
            --threshold $TH

################################# Without TH ######################################
# Paths
path_dea=/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA
path_graph=/home/nuezsal/Omics_integration/RNA_seq/Results/DESeq_graphs
TH=0


# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript DEA_transcripts.R \
            --input $path_in \
            --metadata $path_metadata \
            --output $path_dea \
            --alpha $alpha \
            --specie "cume" \
            --project "Omics_project" \
            --graphs $path_graph \
            --threshold $TH
