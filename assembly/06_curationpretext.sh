#!/bin/bash --login

#Requested resources:
#SBATCH --job-name=nextflow-master
#SBATCH --time=24:00:00
#SBATCH --partition=work
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=pawsey1088

module load nextflow/25.04.6
module load singularity/4.1.0-slurm

nextflow run sanger-tol/curationpretext -r 1.5.1 \
  --input $(pwd)/04.build/scaffolds.fa \
  --cram $(pwd)/cram \
  --reads $(pwd)/reads \
  --read_type hifi \
  --sample santarcticus \
  --outdir $(pwd)/santarcticus_cpxt/ \
  --aligner minimap2 \
  -profile singularity,pawsey_setonix \
  -resume \
  -c ira.config
