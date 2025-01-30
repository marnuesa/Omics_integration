#!/bin/bash

#SBATCH --job-name=ident_miRNAs         # Job name to show with squeue
#SBATCH --output=ident_miRNAs_%j.out    # Output file
#SBATCH --ntasks=1                      # Maximum number of cores to use
#SBATCH --time=00-01:00:00              # Time limit to execute the job
#SBATCH --mem=2G                        # Required Memory
#SBATCH --cpus-per-task=12              # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   miRNAs_annotation_submit.sh
#
#   The program is designed to identify miRNAs from tables of differentially
#   expressed sequences of Cucumis melo using the miRNAS_annotation.sh
#   programm.
#
#   Author: Marta Núñez Salvador
#   Date: 22/09/2024
#   Version: 2.0 
#
#******************************************************************************
############################## With 1 mismatch #################################
# Input paths
# All sequences annotation
path_in=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA/01-DEA_raw
path_mirbase=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/miRBase
path_PmiREN=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/PmiREN
path_sRNAanno=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/sRNAanno
path_ids_table=/home/nuezsal/Omics_integration/sRNA/Additional_info/species_id.csv
mismatches=1

# Ouput paths
path_out=/home/nuezsal/Omics_integration/sRNA/Results/05-Identification_miRNAs

# Threads
num_threads=12

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mismatches $mismatches \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads


# Significative DE sequences annotation
# Input paths
path_in=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA/02-DEA_sig

# Threads
num_threads=12

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mismatches $mismatches \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads 

############################## With 0 mismatch #################################
# Input paths
# All sequences annotation
path_in=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA/01-DEA_raw
path_mirbase=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/miRBase
path_PmiREN=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/PmiREN
path_sRNAanno=/home/nuezsal/Omics_integration/sRNA/Additional_info/02-Mod_databases/sRNAanno
path_ids_table=/home/nuezsal/Omics_integration/sRNA/Additional_info/species_id.csv
mismatches=0

# Ouput paths
path_out=/home/nuezsal/Omics_integration/sRNA/Results/05-Identification_miRNAs_simm

# Threads
num_threads=12

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mismatches $mismatches \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads


# Significative DE sequences annotation
# Input paths
path_in=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA/02-DEA_sig

# Threads
num_threads=12

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mismatches $mismatches \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads 
exit 0

