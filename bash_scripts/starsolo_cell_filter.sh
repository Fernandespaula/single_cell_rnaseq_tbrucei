#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --job-name="<JOB NAME>"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=00:20:00

#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

localcores=$SLURM_CPUS_PER_TASK
localmem=8

module load vital-it
module load UHTS/Aligner/STAR/2.7.10a_alpha_220517
STAR --runThreadN 64 \
--runMode soloCellFiltering <PATH TO RAW MULTIMAPPING DIRECTORY> <PATH TO FILTERED PREFIX> \  
--soloCellFilter EmptyDrops_CR


