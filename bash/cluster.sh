#!/usr/bin/env bash

#SBATCH --job-name=innateness
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100g
#SBATCH --time=80:00:00
#SBATCH --output=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL

bash bash/lmm_job_newcluster.sh -s mouse -p 0.05 -b 10 -c y

