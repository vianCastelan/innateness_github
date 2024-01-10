import pandas as pd

def find_different_rows(file1, file2):
    df_file2 = pd.read_csv(file2)
    genes_in_file2 = set(df_file2.iloc[:, 0])

    df_file1 = pd.read_csv(file1)
    different_rows = df_file1[~df_file1.iloc[:, 0].isin(genes_in_file2)]

    return different_rows

def write_filtered_file(input_file, output_file, filtered_rows):
    filtered_rows.to_csv(output_file, index=False)

def merge_files(file1, file2, output_file):
    df_file1 = pd.read_csv(file1)
    df_file2 = pd.read_csv(file2)

    # Identify columns in file2 but not in file1
    extra_columns = [col for col in df_file2.columns if col not in df_file1.columns]

    # Merge with only unique columns
    merged_df = pd.merge(df_file1, df_file2[extra_columns + ['GeneSymbol']], on='GeneSymbol', how='left')
    merged_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    file1_name = '../data/countmatrix/GSE109125_Gene_count_table.csv'  # replace with your actual file name
    file2_name = '../data/countmatrix/immgen_ULI_RNAseq.csv'  # replace with your actual file name
    output_file_name = '../data/countmatrix/GSE109125_filtered_genes.csv'  # replace with your desired output file name

    different_rows = find_different_rows(file1_name, file2_name)
    write_filtered_file(file2_name, output_file_name, different_rows)

    merge_files(file1_name, file2_name, output_file_name)