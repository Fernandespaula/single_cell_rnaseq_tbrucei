#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --job-name="<JOB NAME>"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

localcores=$SLURM_CPUS_PER_TASK
localmem=128

module load vital-it
module add UHTS/Aligner/STAR/2.7.10a_alpha_220517


STAR --soloType CB_UMI_Complex \
--soloCBmatchWLtype 1MM \
--soloCBwhitelist <PATH TO BARCODE WHITELIST> \
--soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
--soloAdapterMismatchesNmax 2 \
--soloCBposition 0_0_2_-1 3_1_3_8 \
--soloUMIposition 3_9_3_14 \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--readFilesIn <PATH TO R2> <PATH TO R1> \
--genomeDir <PATH TO REFERENCE GENOME> \
--runThreadN 64 \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outFileNamePrefix <PATH TO RESULT DIRECTORY> \
--soloMultiMappers Uniform EM PropUnique Rescue \
--outSAMattributes GX GN gx gn \
--genomeSAindexNbases 11












