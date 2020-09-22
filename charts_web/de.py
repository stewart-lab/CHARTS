import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output
from dash_table import DataTable

from app import app
import load_data
import common

DEFAULT_TUMOR_1 = "PJ025"
DEFAULT_TUMOR_2 = "PJ035"
PAGE_SIZE = 50

@app.callback(
    Output(component_id='de-table-1', component_property='data'),
    [
        Input(component_id='select-tumor-de-1', component_property='value'),
        Input(component_id='select-cluster-de-1', component_property='value')
    ]
)
def update_de_table_1(tumor, cluster):
    return build_table_data(tumor, cluster, 1)

@app.callback(
    Output(component_id='de-table-2', component_property='data'),
    [
        Input(component_id='select-tumor-de-2', component_property='value'),
        Input(component_id='select-cluster-de-2', component_property='value')
    ]
)
def update_de_table_2(tumor, cluster):
    return build_table_data(tumor, cluster, 2)


@app.callback(
    [
        Output(component_id='select-cluster-de-1', component_property='options'),
        Output(component_id='select-cluster-de-1', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-1', component_property='value')
    ]
)
def update_select_cluster_de_1(tumor):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted(set(load_data.load_tumor_clusters(tumor)))
        ],
        '0'
    )


@app.callback(
    [
        Output(component_id='select-cluster-de-2', component_property='options'),
        Output(component_id='select-cluster-de-2', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-2', component_property='value')
    ]
)
def update_select_cluster_de_2(tumor):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted(set(load_data.load_tumor_clusters(tumor)))
        ],
        '0'
    )


def build_table_data(tumor, cluster, table_id):
    genes, log_fc, pval_adj = load_data.load_de(tumor, cluster)
    abs_log_fc = map(abs, log_fc)
    df = pd.DataFrame(
        {
            'de-genes-col-{}'.format(table_id): genes,
            'de-lfc-col-{}'.format(table_id): log_fc,
            'de-alfc-col-{}'.format(table_id): abs_log_fc,
            'de-pval-col-{}'.format(table_id): pval_adj
        }
    )
    return df.to_dict('records')
    

def build_layout():
    layout = dcc.Tab(
        label='Differential Expression',
        children=[
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.Div(['Compare differentially expressed genes between tumor sub-clusters of cells.']),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            dbc.Row([
                dbc.Col([
                    dbc.CardHeader(
                        "Table 1",
                        style={
                            "background-color":"#e3e3e3",
                            "font-weight":"bold",
                            "font-size":"Large",
                            "text-align": "center"
                        }
                    ),
                    dbc.CardBody([
                        html.H6("Select a tumor:"),
                        common.build_tumor_dropdown('select-tumor-de-1'),
                        html.H6("Select a cluster:"),
                        dcc.Dropdown(
                            options=[
                                {'label': 'Cluster {}'.format(clust), 'value': clust}
                                for clust in sorted(set(load_data.load_tumor_clusters(DEFAULT_TUMOR_1)))
                            ],  
                            id='select-cluster-de-1',
                            value='0'
                        ),
                        dbc.Row(html.Hr(), style={'height': '0.5%'}),
                        DataTable(
                            columns=[
                                {
                                    "id": "de-genes-col-1", 
                                    "name": "DE Genes"
                                },
                                {
                                    "id": "de-lfc-col-1",
                                    "name": "log-FC"
                                },
                                {
                                    "id": "de-alfc-col-1",
                                    "name": "Abs. log-FC"
                                },
                                {
                                    "id": "de-pval-col-1",
                                    "name": "FDR"
                                }
                            ],
                            data=build_table_data("PJ025", "0", 1),
                            style_cell_conditional=[
                                {
                                    'textAlign': 'center'
                                }   
                            ],  
                            id='de-table-1',
                            filter_action="native",
                            sort_action="native",
                            page_size=PAGE_SIZE
                        )
                    ])
                ]), 
                dbc.Col([], width=100, style={"width": "15px"}),
                dbc.Col([
                    dbc.CardHeader(
                        "Table 2",
                        style={
                            "background-color":"#e3e3e3",
                            "font-weight":"bold",
                            "font-size":"Large",
                            "text-align": "center"
                        }
                    ),
                    dbc.CardBody([
                        html.H6("Select a tumor:"),
                        common.build_tumor_dropdown('select-tumor-de-2'),
                        html.H6("Select a cluster:"),
                        dcc.Dropdown(
                            options=[
                                {'label': 'Cluster {}'.format(clust), 'value': clust}
                                for clust in sorted(set(load_data.load_tumor_clusters(DEFAULT_TUMOR_2)))
                            ],  
                            id='select-cluster-de-2',
                            value='0'
                        ),
                        dbc.Row(html.Hr(), style={'height': '0.5%'}),
                        DataTable(
                            columns=[
                                {
                                    "id": "de-genes-col-2",
                                    "name": "DE Genes"
                                },
                                {
                                    "id": "de-lfc-col-2",
                                    "name": "log-FC"
                                },
                                {
                                    "id": "de-alfc-col-2",
                                    "name": "Abs. log-FC"
                                },
                                {
                                    "id": "de-pval-col-2",
                                    "name": "FDR"
                                }
                            ],
                            data=build_table_data("PJ025", "0", 1),
                            style_cell_conditional=[
                                {
                                    'textAlign': 'center'
                                }   
                            ],  
                            id='de-table-2',
                            filter_action="native",
                            sort_action="native",
                            page_size=PAGE_SIZE
                        )
                    ])
                ])  
            ])  
        ]   
    )
    return layout
