#!/usr/bin/env bash

#SBATCH --job-name=beta_exhausted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100g
#SBATCH --time=05:00:00
#SBATCH --output=/mnt/bioadhoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/bioadhoc/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL

## load python module
module load python/3.10.10 R/4.2.2
source $HOME/virtualenvs/python3/bin/activate

# parameters
COUNTMATRIX='data/GSE131847/exhaustion.csv'
BETA='output/b_levels/results_mouse_beta_table.csv'
OUTPUT='output/b_scores/'
NAME='GSE141847_scores.tsv'

# Rscript run 
Rscript R/tissue_residency.R

## fix count matrix
# cut -d',' -f1-5 $COUNTMATRIX | head -n 5
# cut -f1-5 "data/single_cell/salmonella_counts.tsv" | head -n 5
# sed -i '1s/^[^,]*/gene/' $COUNTMATRIX ## no quotes
# sed -i '1s/^[^,]*/"gene"/' $COUNTMATRIX  ## quotes
# cut -d',' -f1-5 $COUNTMATRIX | head -n 5

# run
python3 python/innate_score.py -c $COUNTMATRIX \
      -b $BETA \
      -o $OUTPUT \
      -n $NAME