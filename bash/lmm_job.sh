#!/bin/bash

###############################
# beta_lmm #
###############################

#PBS -q default
#PBS -N beta_lmms
#PBS -l nodes=1:ppn=15
#PBS -l mem=80GB
#PBS -l walltime=80:00:00
#PBS -o lmm_out.txt
#PBS -e lmm_err.txt
#PBS -V
#PBS -M vcastelan@lji.org
#PBS -m abe

## activate R4 environment from conda
source activate r4
cd /mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github
## get flag variables
while getopts s:p:b: flag
do
    case ${flag} in
        s) species=${OPTARG};;
        p) pval=${OPTARG};;
        b) beta=${OPTARG};;
    esac
done
echo "----------Running linear mixed models of innateness----------"
echo "Settings:";
echo "species: $species";
echo "sig. p-value: $pval";
echo "sig. beta level: $beta";
## run Rscript 
if [ $species == 'mouse' ]
then
    echo "initiating Rscript"
    Rscript R/rna_m_lmm.R
elif [ $species == 'human' ]
then
    echo "initiating Rscript"
    Rscript R/rna_h_lmm.R $pval $beta
else 
    echo "not correct species, input either mouse or human"
fi







### options for scripts 
#
# bash bash/lmm_job.sh -s mouse -p 0.05 -b 10
#
#
