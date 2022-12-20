#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --job-name="<JOB NAME>"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4:00:00

#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)


localcores=$SLURM_CPUS_PER_TASK
localmem=128

module load vital-it
module add UHTS/Aligner/STAR/2.7.10a_alpha_220517

STAR --soloType CB_UMI_Simple \
--soloUMIfiltering MultiGeneUMI_CR \
--genomeSAsparseD 3 \
--soloUMIlen 12 \
--soloUMIdedup 1MM_CR \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloCBwhitelist <PATH TO WHITELIST> \
--readFilesIn <PATH TO R2> <PATH TO R1> \
--genomeDir <PATH TO STAR GENOME INDEX> \
--runThreadN 128 \
--readFilesCommand zcat \
--clipAdapterType CellRanger4 \
--outFilterScoreMin 30 \
--soloCellFilter EmptyDrops_CR \
--soloCellReadStats Standard \
--quantMode GeneCounts \
--soloBarcodeReadLength 0 \
--outSAMtype BAM Unsorted \
--soloMultiMappers Uniform EM PropUnique Rescue \
--outSAMattributes GX GN gx gn \
--outFileNamePrefix <PATH TO OUTPUT DIRECTORY> \
--genomeSAindexNbases 11 






