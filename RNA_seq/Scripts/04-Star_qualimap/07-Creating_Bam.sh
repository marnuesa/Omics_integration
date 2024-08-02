#!/bin/bash
#SBATCH --ntasks=25               # Number of task
#SBATCH --cpus-per-task=4              # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=3G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o 07-Creating_Bam_%j.log # Standard output goes to this file

#VARIABLES
clean_dir="/home/nuezsal/omics_integration_G3/RNAseq/Data/clean_data"
index_dir="/home/nuezsal/omics_integration_G3/RNAseq/Results/02-Salmon/index_STAR"
result_dir="/home/nuezsal/omics_integration_G3/RNAseq/Results/02-Salmon/STAR"

cd $clean_dir
sample_list=$(ls -d RNA*)

cd $result_dir
for dir in $sample_list;
do
	mkdir ${dir}
done


for dir in $sample_list;
do
  	echo -e "\nSample: $dir\n"
	srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive STAR --genomeDir ${index_dir} \
        --runThreadN 4 \
        --readFilesCommand gunzip -c --readFilesIn ${clean_dir}/${dir}/${dir}_1.fq.gz ${clean_dir}/${dir}/${dir}_2.fq.gz \
        --outFileNamePrefix ${result_dir}/${dir}/${dir}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard & 
done
wait

exit 0
