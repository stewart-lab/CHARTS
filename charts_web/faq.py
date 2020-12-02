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
                    html.H2(children='Frequently asked questions', style={"text-align": "center"}),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H5(children='Can I upload my own single-cell data into CHARTS?'),
                    html.Div(children=[
                        html.Span(children=[
                            """
                            CHARTS enables
                            the exploration of publicly available, published data.  However, one can run the 
                            CHARTS application locally with a custom backend database. For documentation regarding
                            how one can go about processing data with the same pipeline used to process the
                            datasets exposed by CHARTS, please see the
                            """,
                            html.A("GitHub repository", href="https://github.com/stewart-lab/CHARTS"),
                            "."
                        ])
                    ]),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H5(children='Who do I contact to add a new dataset to CHARTS?'),
                    html.Div(children=[
                        html.Span(children=[
                            """
                            If there exists a public single-cell cancer dataset that you would like to see processed
                            and added to CHARTS, please 
                            """,
                            html.A("contact us", href="mailto:mbernstein@morgridge.org"),
                            " and we will work with you to add the dataset of interest."
                        ])
                    ]),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H5(children='How do I decode the cancer type abbreviations found throughout the website (e.g. "GBM")?'),
                    html.Div(children=[
                        html.Span(children=[
                            """
                            These cancer type abbreviations are those used by The Cancer Genome Atlas. To decode these
                            abbreviations, please reference
                            """,
                            html.A("this page", href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations"),
                            "."
                        ])
                    ])
                ],
                style={"padding-left": "5%", "padding-right": "5%", "width": "100%"}
            ),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            footer.LAYOUT
        ],
        style={"padding-left": "10%", "padding-right": "10%", "width": "100%"}
    )
    return layout


