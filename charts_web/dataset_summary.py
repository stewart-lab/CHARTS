import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import h5py
import numpy as np
import json
from dash.dependencies import Input, Output
from dash_table import DataTable

from app import app
import common
import load_data
import nav
import footer

COLNAME_TO_ID = {
    'Tumor': 'data-summary-tumor-col',   
    'Cancer Type': 'data-summary-cancer-type-col',
    'Publication': 'data-summary-pub-col',
    'No. Cells': 'data-summary-num_cells-col'
}

PAGE_SIZE = 50

def build_layout():
    df = load_data.dataset_summary()
    df_new_col = df.copy()
    df_new_col.columns = [COLNAME_TO_ID[x] for x in df.columns]
    layout = html.Div(
        children=[
            nav.LAYOUT,
            html.Div(
                children=[
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H2(children='Datasets Summary', style={"text-align": "center"}),
                    html.Div(
                        children=[
                            """
                            A summary of the tumor datasets exposed by CHARTS.
                            """
                        ]
                    ),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    DataTable(
                        columns=[
                            {
                                "id": COLNAME_TO_ID[col],
                                "name": col
                            }
                            for col in df.columns
                        ],
                        data=df_new_col.to_dict('records'),
                        style_cell_conditional=[
                            {
                                'textAlign': 'center'
                            }
                        ],
                        id='data-summary-table',
                        filter_action="native",
                        sort_action="native",
                        page_size=PAGE_SIZE
                    )
                ]
            ),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            footer.LAYOUT
        ],
        style={"padding-left": "10%", "padding-right": "10%", "width": "100%"}
    )
    return layout


