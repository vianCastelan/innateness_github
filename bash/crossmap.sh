#!/usr/bin/env bash

#SBATCH --job-name=crossmap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30g
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL

## parameters
READDIR="/mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github"
DATADIR="$READDIR/data/bigwig"
OUTDIR="$READDIR/data/mm39bw"
if [ ! -e ${OUTDIR} ]; then
    mkdir $OUTDIR
fi

## access directory
cd $READDIR

## datasets with .bed.gz files 
readarray -t FILES <<< "$(find $DATADIR -type f -name "*.bw")"
for file in "${FILES[@]}"; do
    echo "${file}"
    outfile="${file%.bw}_mm39"
    ## crossmap 
    CrossMap.py bigwig "${READDIR}/data/annotations/mm10ToMm39.over.chain.gz" ${file} ${outfile}
    ## echo 
    echo "${file} done"
done

## move files to correct directory 
readarray -t FILES <<< "$(find $DATADIR -type f -name "*_mm39.bw")"
for file in "${FILES[@]}"; do
    mv $file $OUTDIR
done