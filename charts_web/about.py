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
import nav
import footer

def build_layout():
    layout = html.Div(
        children=[
            nav.LAYOUT,
            html.Div(
                children=[
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H2(children='About', style={"text-align": "center"}),
                    html.Div(
                        children=[
                            """
                            CHARTS is a web application for exploring publicly available single-cell RNA-seq data
                            from human tumors.  CHARTS enables researchers to compare and characterize tumor
                            subpopulations between tumors across datasets and cancer types.
                            """
                        ]
                    ),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.Div(children=[
                        html.Span([
                            """
                            For more information, please see the manuscript:
                            """
                        ]),
                        html.A("Bernstein, M.N. et al. (2020). CHARTS: A web application for characterizing and comparing tumor subpopulations in publicly available single-cell RNA-seq datasets. bioRxiv.", href='https://plot.ly', target="_blank")
                    ]),
                    html.H2(children='Team', style={"text-align": "center"}),
                    #dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.Div([html.Span([
                        "CHARTS was developed at the ", 
                        html.A("Morgridge Institute for Research", href="https://morgridge.org/", target="_blank"),
                        " and the ",
                        html.A("University of Wisconsin - Madison", href="https://www.wisc.edu/", target="_blank"),
                        " by the following team:"
                    ])]),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.Div(children=[
                        html.A("Matthew N. Bernstein", href="https://mbernste.github.io", target="_blank"),
                        html.Br(),
                        html.A("Zijian Ni", href="https://stat.wisc.edu/staff/ni-zijian/", target="_blank"),
                        html.Br(),
                        html.A("Mike Collins", href="https://morgridge.org/profile/michael-collins/", target="_blank"),
                        html.Br(),
                        html.A("Mark E. Burkard", href="https://www.medicine.wisc.edu/hematology-oncology/welcome-burkard-research-group", target="_blank"),
                        html.Br(),
                        html.A("Christina Kendziorski", href="https://www.biostat.wisc.edu/~kendzior/", target="_blank"),
                        html.Br(),
                        html.A("Ron Stewart", href="https://morgridge.org/research/regenerative-biology/bioinformatics/", target="_blank")
                    ], style={"text-align": "center"})
                ],
                style={"padding-left": "5%", "padding-right": "5%", "width": "100%"}
            ),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            footer.LAYOUT
        ],
        style={"padding-left": "10%", "padding-right": "10%", "width": "100%"}
    )

    return layout


