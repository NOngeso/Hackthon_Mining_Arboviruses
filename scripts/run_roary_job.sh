#!/usr/bin/bash

#SBATCH --job-name='roary pangenomics'
#SBATCH --cpus-per-task=32
#SBATCH --mem=132GB
#SBATCH --output=testjob-%j-stdout.log
#SBATCH --error=testjob-%j-stderr.log
#SBATCH --time=224:00:00
#SBATCH --partition=qbatch

echo "Submitting Slurm job"
conda activate pangenomics
roary -e --mafft -p 32 -f /home/kinyugo/pangenomics/output/roary /home/kinyugo/pangenomics/output/prokka/*/*.gff

