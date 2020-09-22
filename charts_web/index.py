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
#import load_data
import dim_reduc
#import clust_compare
import nav
import about
import charts
import dataset_summary
import faq

@app.callback(Output('page-content', 'children'),
            [Input('url', 'pathname')])
def display_page(pathname):
    print(pathname)
    if pathname == '/about':
        return about.build_layout()
    elif pathname == '/data_set_summary':
        return dataset_summary.build_layout()
    elif pathname == '/faq':
        return faq.build_layout()
    else:
        return charts.build_layout()

app.layout = html.Div([
    dcc.Location(id = 'url', refresh = False),
    html.Div(id = 'page-content')
])


if __name__ == '__main__':
    app.run_server(debug=True)
