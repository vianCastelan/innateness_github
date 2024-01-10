#!/bin/bash

###############################
# atac_lmm_job3 #
###############################

#PBS -q default
#PBS -N atac_lmm_job3
#PBS -l nodes=1:ppn=10
#PBS -l mem=100GB
#PBS -l walltime=80:00:00
#PBS -o /mnt/BioAdHoc/Groups/KronenbergLab/vcastelan/innateness/mouse/atac-seq/counts/atac_lmm_job3_out.txt
#PBS -e /mnt/BioAdHoc/Groups/KronenbergLab/vcastelan/innateness/mouse/atac-seq/counts/atac_lmm_job3_err.txt
#PBS -V
#PBS -M vcastelan@lji.org
#PBS -m abe

cd /mnt/BioAdHoc/Groups/KronenbergLab/vcastelan/innateness/mouse/atac-seq/counts/

Rscript atac_lmm3.R

