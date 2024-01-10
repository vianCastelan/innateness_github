#!/usr/bin/env bash

#SBATCH --job-name=gviz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60g
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/BioAdHoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/BioAdHoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL

## load conda environment
module load R/4.2.2

## run Rscript
cd /mnt/BioAdHoc/Groups/KronenbergLab/gascui/innateness_github

Rscript R/gviz_tracks.R

echo "chech results"