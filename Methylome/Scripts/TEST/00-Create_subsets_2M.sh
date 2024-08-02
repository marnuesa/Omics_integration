#!/bin/bash

#SBATCH --job-name=subset								# Job name.
#SBATCH --output=Subset.log								# Standard output and error log.
#SBATCH --partition=short								# Partition (queue)
#SBATCH --ntasks=20									# Run on one mode.
#SBATCH --cpus-per-task=1								# Cpus per task.
#SBATCH --time=1-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=8gb								# Job memory request.


####### VARIABLES
WD="/storage/ncRNA/Projects/Integracion_omicas_melon/Metiloma"
F="$WD/Scripts/TEST/Functions.sh"
summary="$WD/Additional_info/Summary_samples/summary_samples_PRUEBA.tsv"

####### DIRECTORY
mkdir -p $WD/Libraries/SUBSET

####### PIPELINE
echo -e "\n\nCreate subsets:\n"

n_samples=$(wc -l $summary | cut -d" " -f1)
for j in `seq 1 $n_samples`; do
	sample=$(cat $summary | head -n $j | tail -n 1)
	sample_name=$(echo $sample | cut -d" " -f2- | sed 's/ /-/g')
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_subset $sample_name $WD/Libraries &
done
wait


