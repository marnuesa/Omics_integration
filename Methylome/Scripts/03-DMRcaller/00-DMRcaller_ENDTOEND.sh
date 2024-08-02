#!/bin/bash

#SBATCH --job-name=dmrE								# Job name.
#SBATCH --output=DMRs_prediction_ENDTOEND.log						# Standard output and error log.
#SBATCH --partition=long								# Partition (queue)
#SBATCH --ntasks=1									# Run on one mode. 
#SBATCH --cpus-per-task=40								# Number of tasks = cpus. 
#SBATCH --time=4-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=30gb								# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/Integracion_omicas_melon/Metiloma"

####### DIRECTORY
mkdir -p $WD/Results/02-DMRcaller
mkdir -p $WD/Results/02-DMRcaller/ENDTOEND

####### PIPELINE
Rscript 00-DMRcaller.R $WD/Results/02-DMRcaller/ENDTOEND $WD/Results/01-Bismark/ENDTOEND/04-Methylation_extractor $WD/Additional_info/Summary_samples/summary_samples.tsv $SLURM_CPUS_PER_TASK

