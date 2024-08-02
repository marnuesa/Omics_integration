#!/bin/bash

#SBATCH --job-name=bisED								# Job name.
#SBATCH --output=Bismark_End_Default.log								# Standard output and error log.
#SBATCH --partition=long								# Partition (queue)
#SBATCH --ntasks=4									# Run on one mode. 
#SBATCH --cpus-per-task=24								# Number of tasks = cpus. 
#SBATCH --time=6-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=4gb								# Job memory request.


####### MODULES
module load biotools

####### VARIABLES
WD="/storage/ncRNA/Projects/Integracion_omicas_melon/Metiloma"
F="$WD/Scripts/TEST/TEST_END_DEFAULT/Functions.sh"
summary="$WD/Additional_info/Summary_samples/summary_samples_PRUEBA.tsv"

####### DIRECTORY
mkdir -p $WD/Results/01-Bismark
mkdir -p $WD/Results/01-Bismark/TEST
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/01-Genome_preparation
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/02-Alignment
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/03-Deduplication
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/04-Methylation_extractor
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/05-Global_reports
mkdir -p $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/Logs

####### PIPELINE
#### BISMARK: GENOME PREPARATION (INDEXING)
echo -e "\n\n--------------------------------------------------"
echo -e "--------------- GENOME PREPARATION ---------------"
echo -e "--------------------------------------------------\n"

srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Genome_preparation $WD/Additional_info $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT

#### BISMARK: ALIGNMENT
echo -e "\n\n--------------------------------------------------"
echo -e "------------------- ALIGNMENT --------------------"
echo -e "--------------------------------------------------\n"

n_samples=$(wc -l $summary | cut -d" " -f1)
for j in `seq 1 $n_samples`; do
	sample=$(cat $summary | head -n $j | tail -n 1)
	sample_name=$(echo $sample | cut -d" " -f2- | sed 's/ /-/g')
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Alignment_paired_end $sample_name $WD/Libraries $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT &
done
wait

#### BISMARK: DEDUPLICATE
echo -e "\n\n--------------------------------------------------"
echo -e "----------------- DEDUPLICATION ------------------"
echo -e "--------------------------------------------------\n"

n_samples=$(wc -l $summary | cut -d" " -f1)
for j in `seq 1 $n_samples`; do
	sample=$(cat $summary | head -n $j | tail -n 1)
	sample_name=$(echo $sample | cut -d" " -f2- | sed 's/ /-/g')
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Deduplication_paired_end $sample_name $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT &
done
wait

#### BISMARK: METHYLATION EXTRACTOR
echo -e "\n\n--------------------------------------------------"
echo -e "------------- METHYLATION EXTRACTOR --------------"
echo -e "--------------------------------------------------\n"

n_samples=$(wc -l $summary | cut -d" " -f1)
for j in `seq 1 $n_samples`; do
	sample=$(cat $summary | head -n $j | tail -n 1)
	sample_name=$(echo $sample | cut -d" " -f2- | sed 's/ /-/g')
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Methylation_extractor_paired_end $sample_name $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT &
done
wait

#### BISMARK: GLOBAL REPORT
echo -e "\n\n--------------------------------------------------"
echo -e "----------------- GLOBAL REPORT ------------------"
echo -e "--------------------------------------------------\n"

n_samples=$(wc -l $summary | cut -d" " -f1)
for j in `seq 1 $n_samples`; do
	sample=$(cat $summary | head -n $j | tail -n 1)
	sample_name=$(echo $sample | cut -d" " -f2- | sed 's/ /-/g')
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_QC $sample_name $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT &
done
wait

multiqc $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/05-Global_reports/ -o $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/05-Global_reports/ >> $WD/Results/01-Bismark/TEST/TEST_END_DEFAULT/Logs/MultiQC_stdout.log 2>&1

echo -e "\n\n*************** END BISMARK SCRIPT ***************\n\n"


