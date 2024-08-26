#!/bin/bash
#SBATCH --output=DMRs_%j.log                                                    # Standard output and error log.
#SBATCH --qos long                                                              # Partition (queue)
#SBATCH --ntasks=1                                                                      # Run on one mode. 
#SBATCH --cpus-per-task=40                                                              # Number of tasks = cpus. 
#SBATCH --time=4-00:00:00                                                              # Time limit days-hrs:min:sec.
#SBATCH --mem=300gb                                                               # Job memory request.


####### MODULES
module load anaconda
source activate group_sRNA

####### VARIABLES
path="/home/nuezsal/Omics_integration/Methylome"


####### PIPELINE
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript 01-DMRcaller.R $path/Results/03-DMRcaller $path/Results/02-Bismark/04-Methylation_extractor $path/Additional_info/Summary_samples/summary_samples.tsv $SLURM_CPUS_PER_TASK
