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
import dim_reduc
import clust_compare
import nav
import de
import footer

def build_layout():
    layout = html.Div(children=[
        nav.LAYOUT,
        html.Div(children=[
            dbc.Row(html.Hr(), style={'height': '3%'}),
            dbc.Container(fluid=True, children=[
                dcc.Tabs([
                    dim_reduc.build_layout(),
                    de.build_layout()
                    #clust_compare.build_layout()
                ]),
            ], style={"padding-right": "0px"}),
            dbc.Row(html.Hr(), style={'height': '3%'})
        ]),
        footer.LAYOUT
    ], style={"padding-left": "10%", "padding-right": "10%", "width": "100%"}, className="fullpage")
    return layout
