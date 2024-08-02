#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3               # Number of task
#SBATCH --cpus-per-task=2              # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=3G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o 05-salmon_quant_%j.log # Standard output goes to this file


################################################################################
#
#               02-Salmon.sh
#   This program execute Salmon pseudo-aligments to quantify each transcript
#   using a reference genome
#
################################################################################

#PATHs
path_in="/home/nuezsal/Omics_integration/RNA_seq/Libraries/01-Clean_data"
path_out="/home/nuezsal/Omics_integration/RNA_seq/Results/02-Salmon"
path_index="/home/nuezsal/Omics_integration/RNA_seq/Additional_info/Melon_salmon_ind/"

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

cd ${path_in}
samples_list=$(ls)

for sample in $samples_list; do
	echo -e "Processing sample ${sample}\n"
	srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK -Q --exclusive salmon quant -i ${path_index} -l A \
         -1 ${sample}/${sample}_1.fq.gz \
         -2 ${sample}/${sample}_2.fq.gz \
         -p 8 --validateMappings -o ${path_out}/${sample} &
done
wait

exit 0
