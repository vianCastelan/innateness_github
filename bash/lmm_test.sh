#!/usr/bin/env bash

#SBATCH --job-name=lmm_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80g
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/bioadhoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/bioadhoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL


####################################################
# Testing different conditions for PCA coordinates #
####################################################

# module
module load R/4.2.2 

# variables
proj_dir='/mnt/BioAdHoc/Groups/KronenbergLab/gascui/innateness/innateness_github'
pval=0.05
beta=10

## change directory
cd $proj_dir
echo "`pwd -P`"
#======================# ITERATE

# mvgs
yq '.test.mvgs[]' test.yaml | while IFS= read -r mvg; do
  #echo " - $mvg"
  # cells
  yq '.test.cells[]' test.yaml | while IFS= read -r cell; do
    #echo "  - $cell"
    echo "Running lmm test R script with ${mvg} on ${cell}"
    ## R script 
    Rscript R/rna_m_lmm_tests.R $pval $beta $mvg $cell
  done
done


# 
echo "test script has run!"