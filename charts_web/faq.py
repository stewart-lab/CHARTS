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
                            Most cancer type abbreviations are those used by The Cancer Genome Atlas. To decode these
                            abbreviations, please reference
                            """,
                            html.A("this page", href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations"),
                            """. CHARTS also contains data for some cancer types not in the TCGA: high-grade glioma (HGG), 
                            gastrointestinal neuroendocrine cancer (GINET), basal cell carcinoma (BCC), and cutaneous squoumous cell carcinoma (SCC)."""
                        ])
                    ]),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H5(children='Where can I download data for further exploration?'),
                    html.Div(children=[
                        html.Span(children=[
                            """
                            Each tumor's dataset is available to download as tab-separated text files at
                            """,
                            html.A("https://charts.morgridge.org/download", href="https://charts.morgridge.org/download"),
                            ". See ",
                            html.A(" this Jupyter notebook ", href="https://github.com/stewart-lab/CHARTS/blob/master/explore_charts_results.ipynb"),
                            "for an example of how this data can be explored in Python with ",
                            html.A("Scanpy", href="https://scanpy.readthedocs.io/en/stable/"),
                            "."
                        ])
                    ]),
                    dbc.Row(html.Hr(), style={'height': '3%'}),
                    html.H5(children='How do I cite CHARTS?'),
                    html.Div(children=[
                        html.Span(children=[
                            "Please cite the original publication for any of the single-cell datasets used accessed via CHARTS in addition to the following publication:",
                            html.Br(),
                            html.Br(),
                            "Bernstein, M.N., Ni, Z., Collins, M., Burkard, M.E., Kendziorski, C., and Stewart, R. (2020). CHARTS: A web application for characterizing and comparing tumor subpopulations in publicly available single-cell RNA-seq datasets. bioRxiv.",
                            html.A(" https://doi.org/10.1101/2020.09.23.310441", href="https://doi.org/10.1101/2020.09.23.310441"),
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


