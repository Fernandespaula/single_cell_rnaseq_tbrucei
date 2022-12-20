#!/bin/bash
#SBATCH --mail-user=<USER EMAIL>
#SBATCH --mail-type=end,fail
#SBATCH --job-name="<JOB NAME>"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00

#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

module load vital-it
module add UHTS/Aligner/STAR/2.7.10a_alpha_220517

STAR --runMode genomeGenerate --runThreadN 48 \
--genomeDir <PATH TO INDEX DIRECTORY OUTPUT> \
--genomeFastaFiles <PATH TO REFERENCE GENOME FASTA FILE> \
--sjdbGTFfile <PATH TO GTF FILE> \
--genomeSAindexNbases 11





