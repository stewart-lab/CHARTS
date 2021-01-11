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
        Input(component_id='select-cluster-res-de-1', component_property='value'),
        Input(component_id='select-cluster-de-1', component_property='value')
    ]
)
def update_de_table_1(tumor, res, cluster):
    return build_table_data(res, tumor, cluster, 1)

@app.callback(
    Output(component_id='de-table-2', component_property='data'),
    [
        Input(component_id='select-tumor-de-2', component_property='value'),
        Input(component_id='select-cluster-res-de-2', component_property='value'),
        Input(component_id='select-cluster-de-2', component_property='value')
    ]
)
def update_de_table_2(tumor, res, cluster):
    return build_table_data(res, tumor, cluster, 2)


@app.callback(
    [
        Output(component_id='select-cluster-de-1', component_property='options'),
        Output(component_id='select-cluster-de-1', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-1', component_property='value'),
        Input(component_id='select-cluster-res-de-1', component_property='value')
    ]
)
def update_select_cluster_de_1(tumor, res):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted(set(load_data.load_tumor_clusters(tumor, res)))
        ],
        '0'
    )


@app.callback(
    [
        Output(component_id='select-cluster-de-2', component_property='options'),
        Output(component_id='select-cluster-de-2', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-2', component_property='value'),
        Input(component_id='select-cluster-res-de-2', component_property='value')
    ]
)
def update_select_cluster_de_2(tumor, res):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted(set(load_data.load_tumor_clusters(tumor, res)))
        ],
        '0'
    )


def build_table_data(res, tumor, cluster, table_id):
    genes, log_fc, sig_rank = load_data.load_de(res, tumor, cluster)
    abs_log_fc = map(abs, log_fc)
    df = pd.DataFrame(
        {
            'de-genes-col-{}'.format(table_id): genes,
            'de-lfc-col-{}'.format(table_id): log_fc,
            'de-alfc-col-{}'.format(table_id): abs_log_fc,
            'de-sig-col-{}'.format(table_id): sig_rank
        }
    )
    return df.to_dict('records')
    

def build_layout():
    layout = dcc.Tab(
        label='Cluster Markers',
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
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select a tumor:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    common.build_tumor_dropdown('select-tumor-de-1')
                                ])
                            ])
                        ]),
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select cluster resolution:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    common.build_cluster_res_dropdown('select-cluster-res-de-1')
                                ])
                            ])
                        ]),
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select a cluster:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    dcc.Dropdown(
                                        options=[
                                            {'label': 'Cluster {}'.format(clust), 'value': clust}
                                            for clust in sorted(set(load_data.load_tumor_clusters(DEFAULT_TUMOR_1, 1)))
                                        ],  
                                        id='select-cluster-de-1',
                                        value='0'
                                    )
                                ])
                            ])
                        ]),
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
                                    "id": "de-sig-col-1",
                                    "name": "Significance Rank"
                                }
                            ],
                            data=build_table_data("1", "PJ025", "0", 1), # TODO add resolution
                            style_cell_conditional=[
                                {
                                    'textAlign': 'center'
                                }   
                            ],  
                            id='de-table-1',
                            filter_action="native",
                            sort_action="native",
                            page_size=PAGE_SIZE,
                            export_format='csv',
                            export_headers='display'
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
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select a tumor:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    common.build_tumor_dropdown('select-tumor-de-2')
                                ])
                            ])
                        ]),
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select cluster resolution:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    common.build_cluster_res_dropdown('select-cluster-res-de-2')
                                ])
                            ])
                        ]),
                        html.Div([
                            dbc.Row(children=[
                                dbc.Col([
                                    html.H6("Select a cluster:")
                                ], width=100, style={"width": "40%"}),
                                dbc.Col([
                                    dcc.Dropdown(
                                        options=[
                                            {'label': 'Cluster {}'.format(clust), 'value': clust}
                                            for clust in sorted(set(load_data.load_tumor_clusters(DEFAULT_TUMOR_1, 1)))
                                        ],
                                        id='select-cluster-de-2',
                                        value='0'
                                    )
                                ])
                            ])
                        ]),
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
                                    "id": "de-sig-col-2",
                                    "name": "Significance Rank"
                                }
                            ],
                            data=build_table_data("1", "PJ025", "0", 1), # TODO add resolution
                            style_cell_conditional=[
                                {
                                    'textAlign': 'center'
                                }   
                            ],  
                            id='de-table-2',
                            filter_action="native",
                            sort_action="native",
                            page_size=PAGE_SIZE,
                            export_format='csv',
                            export_headers='display'
                        )
                    ])
                ])  
            ])  
        ]   
    )
    return layout
