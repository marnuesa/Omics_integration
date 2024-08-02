#!/bin/bash
#SBATCH --ntasks=10               # Number of task
#SBATCH --cpus-per-task=6               # Number of cpus
#SBATCH -t 0-03:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=10G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o post-fastqc_%j.log # Standard output goes to this file

################################################################################
#                       00-PreFastQC.sh
#
#   This script execute fastqc to do a preliminar quality control with
#   Raw data
#
################################################################################

#MODULE
module load biotools

#PATHs
path_in='/home/nuezsal/Omics_integration/Methylome/Libraries/00-Raw_data'
path_out='/home/nuezsal/Omics_integration/Methylome/Results/00-PreFastQC'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "El directorio de salida no existe. Creando directorio..."
    mkdir "$path_out"
    echo "Directorio creado correctamente."
else
    echo "El directorio de salida ya existe."
fi


#EXECUTE FASTQC
cd $path_in

sample_list=$(ls)

for sample in $sample_list;
do
	# Execute
	srun -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastqc ${path_in}/$sample/*fq.gz -t 6 -o ${path_out} &
done
wait

exit 0
