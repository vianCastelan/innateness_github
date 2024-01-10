#!/usr/bin/python3

# libraries
import numpy as np
import pandas as pd

def run_innate(input_file):
    ## Process input file and generate output
    #output_file_path = "data/output.csv" 
    ## read tables into pandas
    cells = input_file
    if not cells.empty:
        cell_types = pd.DataFrame(cells)
    else:
        raise ValueError('cannot read file for some reason')
    beta_path = "data/results_mouse_beta_filt.tsv"
    innate = pd.read_table(beta_path, index_col=0, delimiter='\t')
    ## Arrange columns and indexes 
    cols = cell_types.columns
    if cols[0] != 'gene':
        cols.values[0] = 'gene'
    for i in range(len(cols)):
        if cols[i] in cols[:1]:
            count = 1
            while f"{cols[i]}_{count}" in cols:
                count += 1
            cols.values[i] = f"{cols[i]}_{count}"
    cell_types.columns = cols
    if cell_types.index.name != 'gene':
        cell_types = cell_types.set_index('gene')
    ## filtering 
    vc = innate['gene'].value_counts()
    #vc[vc > 1].shape
    innate_filt = innate[innate['gene'].isin(vc[vc ==1].index.tolist())]
    tdf = innate_filt.set_index("gene")[['beta']]
    cell_types_filt = cell_types[cell_types.index.isin(tdf.index.tolist())]
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
    # aggregate results
    aggregate_score_cell_types = cell_types_filt_trans.sum(axis = 0)
    aggregate_score_cell_types = aggregate_score_cell_types.to_frame("innateness_score")
    aggregate_score_cell_types['category'] = "cell_types"
    aggregate_score_cell_types["cell_type"] = aggregate_score_cell_types.index
    aggregate_score_cell_types = aggregate_score_cell_types.reset_index(drop = True)
    #aggregate_score_cell_types.to_csv(output_file_path, sep = ',', index = False)
    return aggregate_score_cell_types
