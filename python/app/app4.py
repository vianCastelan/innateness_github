#!/usr/bin/env python3

import dash

from dash.dependencies import Input, Output, State
from dash import dcc, html, dash_table

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

## parameters
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

## input data 
innate_levels = pd.read_csv('output/b_levels/results_mouse_beta_filt.tsv', sep = '\t')
immgen_table = pd.read_csv("data/countmatrix/immgen_ULI_RNAseq.csv")
immgen_meta = pd.read_csv("data/metadata/immgen_ULI_RNAseq_metadata.csv")

## fix gene name 
immgen_table.rename(columns={"Unnamed: 0":"gene_id"}, inplace=True)

## filter metadata
immgen_meta = immgen_meta[immgen_meta.eval("tissue == 'Spleen' & activation == 'Naive'")]
d_columns = immgen_meta['names'].to_list()
d_columns = [value.strip('gene_idgene_id') for value in d_columns]
d_columns.append('gene_id')

## filter ImmGEN dataset
immgen_table = immgen_table.loc[:, immgen_table.columns.isin(d_columns)]



## create app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "Innateness Dashboard 1.0.0"

server = app.server

## define sections 
def description_card():
    """
    
    :return: A Div contaning dashboard title & descripions.
    """
    return html.Div(
        id="description-card",
        children=[
            html.H5("Innateness Dashboard"),
            html.H3("Welcome to the Mouse Innateness Dashboard"),
            html.Div(
                id="intro",
                children="Explore beta levels of mouse innateness linear mixed models and use to analyze your own transcriptomic datasets."
            ),
        ],
    )


def generate_control_card():
    """
    
    :return: A Div contaning controls for graphs 
    """
    return html.Div(
        id="control-table",
        children=[
            html.H3("Select gene:"),
            dash_table.DataTable(
                id='datatable',
                columns=[{"name":i, "id": i} for i in innate_levels.columns],
                data=innate_levels.to_dict('records'),
                row_selectable="single",
                sort_action="native",
                page_size=10
            ),
        ]
    )

# def generate_graph(gene, reset):
#     """
#     :param: gene: which gene is to be plotted by selecting on the table.


#     :return: ImmGEN dataset results for each cell type.
#     """
    
#     gene = "Gzma"

#     ## indexed vectors
#     filtered_df = immgen_table[immgen_table['gene_id']== gene].drop(columns='gene_id')
#     df_melt = filtered_df.melt(var_name='names', value_name='expression')
#     df_merged = pd.merge(df_melt, immgen_meta, on='names')

#     ## plot 
#     boxplot = px.box(df_merged, x='cell_type',y='expression')
#     return boxplot


#===============# app container #===============# 
app.layout = html.Div(
    id="app-container",
    children=[
        # Store original data and selected rows 
        dcc.Store(id='data-store', data=immgen_table.to_dict('records')),
        dcc.Store(id='selected-rows-store', data=[]),
        # Banner
        html.Div(
            id="banner",
            className="banner",
            children=[html.Img(src=app.get_asset_url("plotly_logo.png"))],
        ),
        # Left column
        html.Div(
            id="left-column",
            className="four columns",
            children=[description_card(), generate_control_card()]
            + [
                html.Div(
                    ["initial child"], id="output-clientside", style={"display": "none"}
                )
            ],
        ),
        # right column
        html.Div(
            id="right-column",
            className="eight columns",
            children=[
                ## boxplot 
                html.Div(
                    id="immgen_boxplot",
                    children=[
                        html.B("Gene expression"),
                        dcc.Graph(id="innateness_immgen")
                    ]
                )
            ]

        )
    ]

)

## callback for stored rows 
@app.callback(
        Output('selected-rows-store','data'),
        Input('datatable','derived_virtual_selected_rows'),
        State('data-store','data')
)

def store_selected_rows(derived_virtual_selected_rows, data):
    ## store the selected rows indices
    return derived_virtual_selected_rows or []


## app callbacks for table and figure 
@app.callback(
    Output("innateness_immgen", "figure"),
    Input("selected-rows-store","data"),
    State('data-store','data')
)
def update_graph(selected_rows,data):
    if selected_rows is None or len(selected_rows) == 0:
        return {'data': []}
    else:
        specific_row = pd.DataFrame.from_records(data).iloc[selected_rows]
        
        # indexed vectors
        gene = specific_row['gene_id']
        df_melt = specific_row.melt(var_name='names', value_name='expression')
        df_merged = pd.merge(df_melt, immgen_meta, on='names')

        ## figure 
        box_fig = px.box(
            data_frame=df_merged,
            title=f'Boxplot of selected values for {gene}',
            x='cell_type',
            y='expression',
            points='all'
        )

        # Update layout of the box plot
        box_fig.update_layout(
            boxmode='group',  # Group mode for box plots
            showlegend=True
        )

        ## filter per gene 
        return box_fig


# Run the server
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')