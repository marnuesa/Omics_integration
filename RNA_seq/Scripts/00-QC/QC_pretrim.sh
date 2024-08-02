#!/bin/bash
#SBATCH --ntasks=10               # Number of task
#SBATCH --cpus-per-task=2                # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=5G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o pretrim_fastq_%j.log  # Standard output goes to this file

################################################################################
#                       QC_pretrim.sh
#
#   This script execute fastqc to do a preliminar quality control with
#   Raw data
#
################################################################################

#MODULE
module load biotools

#PATHs
path_in='/home/nuezsal/Omics_integration/RNA_seq/Libraries/00-Raw_data'
path_out='/home/nuezsal/Omics_integration/RNA_seq/Results/00-QC/QC_pretrim'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "The output field does not exist. Creating..."
    mkdir "$path_out"
    echo "Successfully created field."
else
    echo "Output field already exist."
fi

#EXECUTE FASTQC
cd $path_in

sample_list=$(ls)

for sample in $sample_list;
do
	echo -e "\nSample: $sample\n"
	# Library's quality control before trimming
	srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastqc ${path_in}/${sample}/*.fq.gz -t 2 -o "$path_out" &
done
wait

exit 0
