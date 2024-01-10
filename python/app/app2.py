import dash 
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from datetime import datetime
import numpy as np

app = dash.Dash(__name__)

## components


## input data 
innate_levels = pd.read_csv('output/b_levels/results_mouse_beta_filt.tsv', sep = '\t')
#exp_table = pd.read_csv('')

## app layout
app.layout = html.Div(
    children=[
        html.Img(src='https://www.lji.org/wp-content/uploads/2020/08/building_south_hero-e1663176149609.jpg',
                style={'width':'50%','height':'50%'}),
        html.H1("New Dashboard for Innateness"),
        html.Span(children=[
            f"Prepared: {datetime.now().date()}", 
            html.Br(),
            " by ", html.B("Gabriel Ascui and Viankail Cedillo "),
            html.Br(),
            html.I("Kronenberg Lab")
        ]),
        ## Dropdown callback 
        dcc.Dropdown(id='title_dd',
             options=[{'label':'Title 1','value':'Title 1'},
                       {'label':'Title 2','value':'Title 2'}]),
        dcc.Graph(id='my_graph'),
        html.Div(
            children=[
                html.H2("This box will be Input data "),
                html.H3("User upload expression dataset button"),
                html.Br(),
                html.H3("User upload metadata per sample button")],
            style={'background-color': 'rgb(224,255,252)',
                   'height': 250, 'width': 500}),
        html.Div(
            children=[
                html.H2("This box will be for some graphs maybe "),
                html.H3("Select specific genes to plot")],
            style={'background-color': 'lightblue',
                   'height': 250, 'width': 250,
                   'border':'5px dotted red',
                   'padding':'100px'})
    ]
)
## callbacks
@app.callback(
    Output(component_id='my_graph',
           component_property='figure'),
    Input(component_id='title_dd',
          component_property='value')
)
def update_plot(selection):
    title = "None Selected"
    if selection:
        title = selection

    # figure 
    bar_fig = px.scatter(
        data_frame=innate_levels,
        title=f'{title}',
        x='beta',y='pval',
        log_y=True,
        hover_data='gene'
    )
    bar_fig.update_layout(
        yaxis=dict(autorange="reversed")
    )
    return bar_fig



if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')