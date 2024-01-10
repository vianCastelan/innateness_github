#!/usr/bin/env bash

#SBATCH --job-name=beta_lmm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100g
#SBATCH --time=80:00:00
#SBATCH --output=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.out
#SBATCH --error=/mnt/bioadhoc-temp/Groups/KronenbergLab/%u/slurm_logs/job-%j.err
#SBATCH --mail-type=ALL

###############################
# beta_lmm #
###############################

## activate R4 and python3 environment from modules
module load R/4.2.2 R-libs/4.2.2/seurat-v5 python/3.10.10
source virtualenvs/python3/bin/activate

## access directory 
cd /mnt/bioadhoc-temp/Groups/KronenbergLab/gascui/innateness_github
## get flag variables
while getopts s:p:b:c: flag
do
    case ${flag} in
        s) species=${OPTARG};;
        p) pval=${OPTARG};;
        b) beta=${OPTARG};;
        c) atac=${OPTARG};;
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
    if [ $atac == 'y' ]
    then
        echo "running ATAC-seq analysis as well"
        Rscript R/atac_m_lmm.R
    else 
        echo "ATAC-seq analysis will not be ran"
    fi
elif [ $species == 'human' ]
then
    echo "initiating Rscript for human innateness"
    Rscript R/rna_h_lmm.R $pval $beta
else 
    echo "not correct species, input either mouse or human"
fi

#=============#  R scripts #=============# 
echo "------------------Running downstream analysis in R scripts------------------"
Rscript R/PCA_plots.R  ## plot PCA coordinates
Rscript R/b_levels.R  ## plot b_levels for Fig. 1
Rscript R/volcano_plot.R ## Fig. S1 
Rscript R/venn_plots.R ## generate venn plots for Fig. 1 
Rscript R/ipa_plots_innateness.R ## Fig. 2
Rscript R/ppp.R ## supplementary ? KEGG
# Rscript R/tissue_residency.R

### Plot specific genes
GENE="Nfat1c"
Rscript R/boxplot.R $GENE

#=============# python scripts #=============# 
echo "------------------Running downstream analysis in python3 scripts------------------"

## NKT + flavell + goldrath 
TESTDIR="data/test_data"

for file in "$TESTDIR"/*; do
    if [ -f "$file" ]; then # check if file
        echo "Processing: $file"
        # extract file name 
        fname=$(basename "$file" _counts.tsv)
        echo "output file: ${fname}_scores.tsv"
        # run python script 
        python3 python/innate_score.py -c "$file" \
              -b output/b_levels/results_mouse_beta_table_filt.csv \
              -o output/b_scores/ -n "${fname}_scores.tsv"
    fi
done

## immgen
python3 python/innate_score.py -c data/countmatrix/immgen_ULI_RNAseq.csv \
      -b output/b_levels/results_mouse_beta_table_filt.csv \
      -o output/b_scores/ -n immgen_scores.tsv

## single cell MAIT cells
TESTDIR="data/single_cell"

for file in "$TESTDIR"/*; do
    if [ -f "$file" ]; then 
        echo "Processing: $file"
        # extract file name 
        fname=$(basename "$file" _counts.tsv)
        echo "output file: ${fname}_scores.tsv"
        # run python script 
        python3 python/innate_score.py -c "$file" \
              -b output/b_levels/results_mouse_beta_table_filt.csv \
              -o output/b_scores/ -n "${fname}_sc_scores.tsv"
    fi
done

## run test graphs on all test datasets:
Rscript R/b_scores_test.R


# ATAC-seq
if [ $atac == 'y' ]
    then
    echo "Running Downstream ATAC-seq analysis"
    Rscript R/trackplot.R
    Rscript R/gviz_tracks.R
fi


### Run T helper scripts
echo "------------------Running downstream analysis in R scripts for T helper analysis ------------------"
Rscript R/rna_m_thelper_lmm.R
Rscript R/thelper.R
Rscript R/heatmap.R
Rscript R/thelper_comparison.R


## terminate
echo "All processes for innateness LMM have been completed. Check output files and stdout/stderr"

### options for scripts 
#
# bash bash/lmm_job_newcluster.sh -s mouse -p 0.05 -b 10 -c y
#
#
