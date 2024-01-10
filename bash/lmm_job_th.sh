#!/bin/bash

###############################
# beta_lmm_th #
###############################

#PBS -q default
#PBS -N beta_lmms
#PBS -l nodes=1:ppn=15
#PBS -l mem=80GB
#PBS -l walltime=80:00:00
#PBS -o lmm_th_out.txt
#PBS -e lmm_th_err.txt
#PBS -V
#PBS -M gascui@lji.org
#PBS -m abe

## activate R4 environment from conda
source activate r4
cd /mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github
echo "Generating beta tables"
Rscript R/rna_m_thelper_lmm.R





