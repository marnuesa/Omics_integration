#!/bin/bash

#SBATCH --job-name=counts_d         # Job name to show with squeue
#SBATCH --output=counts_d_%j.out    # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=04-00:00:00          # Time limit to execute the job
#SBATCH --mem=2G            # Required Memory per core
#SBATCH --cpus-per-task=10           # CPUs assigned per task.
#SBATCH --qos=medium                # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   04-Counts_and_sRNADatabase_submit.sh
#
#   This program executes the 04-Counts_and_sRNADatabase.py program in
#   parallel to generate the absolute counts and Reads Per Million
#   (RPM) tables for a specific project using the trimmed and filtered
#   libraries, also calculating the averages of both types of counts in
#   the different replicates of each condition for the sequences analysed.
#
#   Author: Marta Núñez Salvador
#   Date: 20/09/2024
#   Version: 1.1 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_lib='/home/nuezsal/Omics_integration/sRNA/Libraries/02-Clean_data_fasta'
path_results='/home/nuezsal/Omics_integration/sRNA/Results'
path_rnacentral='/home/nuezsal/Omics_integration/sRNA/Additional_info/01-RNAcentral/RNAcenral_viridiplantae_snrna_snorna_rrna_trna.fasta'

# Execute
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive python3 sRNA_counts.py \
            --path-project $path_lib \
            --path-results $path_results \
            --threads $SLURM_CPUS_PER_TASK \
            --path-rnacentral $path_rnacentral

