#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 36 processor core(s) per node
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="<JOB NAME>"
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

module load vital-it
module add UHTS/SingleCell/cellranger/6.0.1

localcores=$SLURM_CPUS_PER_TASK

cellranger mkref --genome=<NAME OUTPUT FOLDER> \
--fasta=<PATH TO GENOME REFERENCE FASTA FILE> \
--genes=<PATH TO GTF FILE> 


