import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import plotly.figure_factory as ff
import numpy as np
from scipy.spatial.distance import pdist, squareform
import dash_bio as dashbio
from clustergrammer import Network
from collections import defaultdict

from app import app
import load_data
import common
import dash_clustergrammer

FIG_DIM = 600

import json
import pandas as df



#with open('example_clustergrammer.json', 'r') as f:
#    j = json.load(f)
#d = np.array(np.random.rand(100,100))
#l = [str(x) for x in range(100)]
#df = pd.DataFrame(data=d, index=l, columns=l)


def build_layout():
    df = load_data.load_gsva_compare_cluster('hallmark')

    # TODO THIS NEEDS TO BE CLEANED UP!!!!!!
    cat_to_true = defaultdict(lambda: [])
    for clust in df.index:
        if 'PJ030' in clust:
            cat_to_true['LGG'].append(clust)
        elif 'PJ' in clust:
            cat_to_true['GBM'].append(clust)
        elif 'LX' in clust:
            cat_to_true['LUAD'].append(clust)
        elif 'GSE146026' in clust:
            cat_to_true['OV'].append(clust)
        elif 'GSE72056' in clust:
            cat_to_true['SKCM'].append(clust)
        elif 'GSE103322' in clust:
            cat_to_true['HNSC'].append(clust)
        elif 'GSE111672' in clust:
            cat_to_true['PAAD'].append(clust)

    cats = [{
        'title': 'Cancer Type',
        'cats': {
            k: v
            for k,v in cat_to_true.items()
        }
    }]

    net = Network()
    net.load_df(df)
    net.add_cats('row', cats)
    net.make_clust()

    layout = dcc.Tab(
        label='Cluster Comparison',
        children=[dbc.Container(
            fluid=True, 
            children=[
                html.Link(
                    rel='stylesheet',
                    href='./static/custom.css'
                ),
                dash_clustergrammer.DashClustergrammer(
                    id='cgram-component',
                    label='',
                    network_data=net.viz
                )
            ]
        )]
    )
    return layout

