#!/usr/bin/python3

# libraries
import numpy as np
import pandas as pd
import os
import yaml
import argparse
import re

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-c", "--countmatrix", help = "Transcriptome countmatrix to calculate innateness score")
#parser.add_argument("-m","--metadata", help = "Metadata for input samples to make groups")
parser.add_argument("-b", "--betatable", help = "Beta levels table")
parser.add_argument("-o", "--output", help = "Output directory")
parser.add_argument("-n", "--name", help="filename for output innateness scores")
#parser.add_argument("-s", "--singlecell", help="optional for reading single cell count matrices (type yes)")
# Read YAML file
with open("config.yaml", 'r') as stream:
    config_list = yaml.safe_load(stream)

# parameters
if os.path.exists(config_list['innate_dir']['local']):
    print("working locally")
    innate_path = config_list['innate_dir']['local']
else:
    print("working on server")
    innate_path = config_list['innate_dir']['server']

# parse arguments
args = parser.parse_args()

# read tables into pandas
if __name__ == '__main__':
    cells_path = args.countmatrix
    beta_path = args.betatable
    innate = pd.read_table(beta_path, index_col=0, delimiter=",")
    #single_cell = args.singlecell

if not os.path.exists(cells_path):
    raise FileNotFoundError(f"The file '{cells_path}' does not exist.")

file_match = re.match(r'^(.+)\.(gz|tsv|csv)$', cells_path)

if file_match:
    base_name, extension = file_match.groups()
    if extension == 'gz':
        if base_name.endswith('.tsv'):
            cell_types = pd.read_csv(cells_path, index_col=0, compression='gzip', sep='\t')
        elif base_name.endswith('.csv'):
            cell_types = pd.read_csv(cells_path, index_col=0, compression='gzip', sep=',')
        else:
            raise ValueError("Unsupported file extension. Please provide a .tsv, .csv, or .gz file.")
    elif extension == 'tsv':
        cell_types = pd.read_csv(cells_path, index_col=0, sep='\t')
    elif extension == 'csv':
        cell_types = pd.read_table(cells_path, index_col=0, delimiter=',')
else:
    raise ValueError("Invalid file name format. Please provide a valid file path.")

## verbose
print("reading files from:",cells_path)
print("The dimension of your expression matrix are:", cell_types.shape)

## Arrange columns and indexes
if cell_types.index.name != 'gene':
    cell_types = cell_types.set_index('gene')

print(cell_types.head())
print("The dimension of your beta matrix are:", innate.shape)

# filtering
vc = innate['gene'].value_counts()
vc[vc > 1].shape
innate_filt = innate[innate['gene'].isin(vc[vc ==1].index.tolist())]
tdf = innate_filt.set_index("gene")[['beta_norm']]
print(tdf.head())
cell_types_filt = cell_types[cell_types.index.isin(tdf.index.tolist())]
print("filtered genes x samples:", cell_types_filt.shape)
cell_types_filt = cell_types_filt.join(tdf, how = 'left')
# multiplication for all genes: woosh!
#‌　 ∧＿∧　
#（。·ω·。)つ━☆·*。
#⊂　　 ノ 　　　·゜+.
#　しーＪ　　　°。+ *´¨)
#　　　　　　　　　.· ´¸.·*´¨) ¸.·*¨)
#　　　　　　　　　　(¸.·´ (¸.·'*  ☆
dfs = []
for c in cell_types_filt.columns.tolist():
    if c != "beta_norm":
        df = np.multiply(cell_types_filt[c],cell_types_filt["beta_norm"])
        df = df.to_frame(c)
        dfs.append(df)
cell_types_filt_trans = pd.concat(dfs, axis = 1)
# aggregate
aggregate_score_cell_types = cell_types_filt_trans.sum(axis = 0)
aggregate_score_cell_types = aggregate_score_cell_types.to_frame("innateness_score")
aggregate_score_cell_types['category'] = "cell_types"
aggregate_score_cell_types["cell_type"] = aggregate_score_cell_types.index
aggregate_score_cell_types = aggregate_score_cell_types.reset_index(drop = True)
#
print("This is what the output looks like: ")
print(aggregate_score_cell_types.head())
# export data
if args.name is not None:
    fn_out = args.output + args.name
else:
    fn_out = args.output + 'innateness_scores_generic.tsv'
## write data
aggregate_score_cell_types.to_csv(fn_out, sep = '\t', index = False)
# python3 python/innate_score.py -c data/countmatrix/immgen_ULI_RNAseq.csv \ 
#       -b output/b_levels/results_mouse_beta_table.tsv \
#       -o output/b_scores/
# python3 python/innate_score.py -c data/countmatrix/immgen_ULI_RNAseq.csv \
#       -b output/b_levels/results_mouse_beta_atac.tsv \
#       -o output/b_scores/
#

# (⁠ノ⁠ಠ⁠益⁠ಠ⁠)⁠ノ⁠彡⁠┻⁠━⁠┻
# Test data
# python3 python/innate_score.py -c data/test_data/GSE157225_counts.tsv \ 
#       -b output/b_levels/results_mouse_beta.tsv \
#       -o output/b_scores/ -n GSE157225_scores.tsv
