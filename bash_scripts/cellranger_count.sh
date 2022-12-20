#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --job-name="<JOB NAME>"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

module load vital-it
module add UHTS/SingleCell/cellranger/6.0.1

localcores=$SLURM_CPUS_PER_TASK
localmem=64

cd <PATH TO FINAL/RESULT DIRECTORY>

cellranger count --id=<NAME> \
--fastqs=<PATH TO RAW DATA> \
--sample=<SAMPLE DIRECTORY> \
--transcriptome=<PATH TO REFERENCE GENOME DIRECTORY> \
--no-bam \
--r1-length 26 \
--nosecondary \
--localcores="${localcores}" \
--localmem="${localmem}"
