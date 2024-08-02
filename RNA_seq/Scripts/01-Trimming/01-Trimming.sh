#!/bin/bash
#SBATCH --ntasks=3               # Number of task
#SBATCH --cpus-per-task=2              # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=3G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o 01-trimming_fastp_%j.log    # Standard output goes to this file

################################################################################
#
#               01-Trimming.sh
#   This program trims the sequences from a
#   specific library using the fastp program.
#   
#
#   Arguments:
#       Fastq.gz file path
#       Output directory path
#
################################################################################

# MODULE
module load biotools

#PATHs
path_in="/home/nuezsal/Omics_integration/RNA_seq/Libraries/00-Raw_data"
path_out="/home/nuezsal/Omics_integration/RNA_seq/Libraries/01-Clean_data"

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

cd $path_in

samples_list=$(ls)

# Iterate files list
for sample in $samples_list;
do
	
	# Trimming
	echo -e "\nSample: $sample\n"
	srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastp -i ${path_in}/${sample}/${sample}_1.fq.gz -I ${path_in}/${sample}/${sample}_2.fq.gz \
	-o ${path_out}/${sample}_1.fq.gz -O ${path_out}/${sample}_2.fq.gz \
	--correction --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
	--n_base_limit 0 &

done
wait

# Create directories to each pair of samples
cd $path_out

for sample in $samples_list;
do
	mkdir ${sample}
	mv ${sample}* ${sample}/

done
exit 0
