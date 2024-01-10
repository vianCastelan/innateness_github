#!/home/gascui/miniconda3/bin/python3
# libraries
import pandas as pd
import sys
import os
import yaml
import argparse

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-c", "--countmatrix", help = "Transcriptome countmatrix to calculate innateness score")
#parser.add_argument("-m","--metadata", help = "Metadata for input samples to make groups")
parser.add_argument("-b", "--betatable", help = "Beta levels table")
parser.add_argument("-o", "--output", help = "Output directory")

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

# read arguments
args = parser.parse_args()

# read tables (These need to be .tsv files) into pandas
if __name__ == '__main__':
    cells_path = args.countmatrix
    cell_types = pd.read_table(cells_path, index_col = 0, delimiter=",")
    beta_path = args.betatable
    innate = pd.read_table(beta_path, index_col = 0, delimiter="\t")
#
print("reading files from:",cells_path)
print("The dimension of your expression matrix are:", cell_types.shape)
print("The dimension of your beta matrix are:", innate.shape)

#
#cell_types = cell_types.set_index('gene')
#cell_types = cell_types.drop(['Excel_gene'], axis=1)
# filtering
vc = innate['gene'].value_counts()
vc[vc > 1].shape
innate_filt = innate[innate['gene'].isin(vc[vc ==1].index.tolist())]
tdf = innate_filt.set_index("gene")[['beta']]
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
    if c != "beta":
        df = cell_types_filt[c] * cell_types_filt["beta"]
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
fn_out = args.output + 'innateness_scores_human.tsv'
aggregate_score_cell_types.to_csv(fn_out, sep = '\t', index = False)
## mouse 
# python3 python/innate_score.py -c data/countmatrix/immgen_ULI_RNAseq.csv \ 
#       -b output/b_levels/results_mouse_beta_table.tsv \
#       -o output/b_scores/
## human
# python3 python/innate_score_human.py -c data/countmatrix/immgen_ULI_RNAseq.csv \
#       -b output/b_levels/human_beta_homologs.tsv \
#       -o output/b_scores/
#
#