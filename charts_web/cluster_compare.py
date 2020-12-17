import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output

from app import app
import load_data
import common

SCATTER_HEIGHT = '550px'

FIG_DIM = 500

# Default settings
DEFAULT_NUM_DIM = 2
DEFAULT_ALGO = 'umap'
DEFAULT_GENE = 'SOX10'
DEFAULT_TUMOR_1 = "PJ025"
DEFAULT_TUMOR_2 = "PJ035"

# Color blind palette from:
# https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
PALETTE = [
    "#006ddb", # color blind
    "#db6d00", # color blind
    "#004949", # color blind
    "#920000", # color blind
    "#ffb6db", # color blind
    "#009292", # color blind
    "#ff6db6", # color blind
    "#490092", # color blind
    "#b66dff", # color blind
    "#6db6ff", # color blind
    "#b6dbff", # color blind
    "#924900", # color blind
    "#24ff24", # color blind
    "#ffff6d", # color blind
    "#000000", # color blind
    "#FE0D0D", # color blind
]

TUMOR_INFO_TEMPLATE = """
    <link type="text/css" rel="Stylesheet" href="https://codepen.io/chriddyp/pen/bWLwgP.css" />
    Cancer type: {typpe}<br>
    Sex: {sex}<br>
    Grade: {grade}<br>
    Stage: {stage}<br>
    Genomic alteration: {gen}<br>
    Publication: <a href={url} target="_blank">{name}</a><br>
"""


@app.callback(
    [
        Output(component_id='select-cluster-dist-1', component_property='options'),
        Output(component_id='select-cluster-dist-1', component_property='value')
    ],
    [
        Input(component_id='select-tumor-dist-1', component_property='value')
    ]
)
def update_select_cluster_de_1(tumor):
    return (
        [{'label': 'All', 'value': 'all'}] + [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted([int(x) for x in set(load_data.load_tumor_clusters(tumor))])
        ],
        'all'
    )

@app.callback(
    [
        Output(component_id='select-cluster-dist-2', component_property='options'),
        Output(component_id='select-cluster-dist-2', component_property='value')
    ],
    [
        Input(component_id='select-tumor-dist-2', component_property='value')
    ]
)
def update_select_cluster_de_2(tumor):
    return (
        [{'label': 'All', 'value': 'all'}] + [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in sorted([int(x) for x in set(load_data.load_tumor_clusters(tumor))])
        ],
        'all'
    )



@app.callback(
    Output(component_id='dist-plot-1', component_property='children'),
    [   
        Input(component_id='select-tumor-dist-1', component_property='value'),
        Input(component_id='select-cluster-dist-1', component_property='value'),
        Input(component_id='color-by-feature-dist-1', component_property='value'),
        Input(component_id='img-format-dist-1', component_property='value'),
    ]
)
def update_dim_reduc_1(tumor_id, clust, gene, img_format):
    return dcc.Graph(
        figure=_build_dist_plot(tumor_id, gene, clust),
        config=_build_plot_config(img_format),
    )


@app.callback(
    Output(component_id='msg-box-dist-1', component_property='srcDoc'),
    [
        Input(component_id='select-tumor-dist-1', component_property='value')
    ]
)
def _update_msg(tum):
    meta = load_data.get_tumor_meta(tum)
    return TUMOR_INFO_TEMPLATE.format(
        url=meta['pub_url'],
        name=meta['pub_name'],
        typpe=meta['cancer_type'],
        sex=meta['sex'],
        grade=meta['grade'],
        stage=meta['stage'],
        gen=meta['genomic_alteration']
    )


@app.callback(
    Output(component_id='msg-box-dist-2', component_property='srcDoc'),
    [
        Input(component_id='select-tumor-dist-2', component_property='value')
    ]
)
def _update_msg(tum):
    meta = load_data.get_tumor_meta(tum)
    return TUMOR_INFO_TEMPLATE.format(
        url=meta['pub_url'],
        name=meta['pub_name'],
        typpe=meta['cancer_type'],
        sex=meta['sex'],
        grade=meta['grade'],
        stage=meta['stage'],
        gen=meta['genomic_alteration']
    )


@app.callback(
    Output(component_id='dist-plot-2', component_property='children'),
    [
        Input(component_id='select-tumor-dist-2', component_property='value'),
        Input(component_id='select-cluster-dist-2', component_property='value'),
        Input(component_id='color-by-feature-dist-2', component_property='value'),
        Input(component_id='img-format-dist-2', component_property='value')
    ]
)
def update_dim_reduc_2(tumor_id, clust, gene, img_format):
    return dcc.Graph(
        figure=_build_dist_plot(tumor_id, gene, clust),
        config=_build_plot_config(img_format),
    )



def _build_dist_plot(tumor_id, gene, clust):
    df = load_data.load_tumor_gene_w_cluster(tumor_id, gene)
    if clust == 'all':
        fig = px.box(df, y="Expression log(TPM+1)", x="Cluster")
    else:
        df = df.loc[df['Cluster'] == int(clust)]
        fig = px.violin(df, y="Expression log(TPM+1)", x="Cluster", box=True, points="all")
    fig.update_layout(
        autosize=True,
        width=FIG_DIM+10,
        height=FIG_DIM
    )
    return fig


def _build_control_panel(plot_num):
    return [
        dbc.Row([
            dbc.Col([], width=100, style={"width": "15px"}),
            dbc.Col([
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Select a tumor to visualize:")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            common.build_tumor_dropdown('select-tumor-dist-{}'.format(plot_num))
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
                                options=[{'label': 'All', 'value': 'all'}] + [
                                    {'label': 'Cluster {}'.format(clust), 'value': clust}
                                    for clust in sorted([int(x) for x in set(load_data.load_tumor_clusters(DEFAULT_TUMOR_1))])
                                ],
                                #style={"width": '60%'},
                                id='select-cluster-dist-{}'.format(plot_num),
                                value='0'
                            )
                        ])
                    ])
                ]),
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Enter gene: ")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            html.Div([
                                dcc.Input(id='color-by-feature-dist-{}'.format(plot_num), value=DEFAULT_GENE)
                            ], id='color-by-feature-container-dist-{}'.format(plot_num))
                        ])
                    ])
                ])
            ])
        ])
    ]


@app.callback(Output("load-out-dist-1", "children"), [Input("loading-input-1", "value")])
def input_triggers_spinner(value):
    time.sleep(1)
    return value


@app.callback(Output("load-out-dist-2", "children"), [Input("loading-input-2", "value")])
def input_triggers_spinner(value):
    time.sleep(1)
    return value


def _build_plot_config(img_format):
    return {
        'displayModeBar': True,
        'toImageButtonOptions': {
            'format': img_format,
            'filename': 'dash_plot'
        },
        "displaylogo": False
    }

def build_layout():
    layout = dcc.Tab(
        label='Cluster Comparison',
        children=[
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.Div(['Compare tumors and subpopulations via violin plots. For each plot, select a tumor, a cluster within that tumor, and a gene.']),
            dbc.Container(fluid=True, children=[
                dbc.Row(html.Hr(), style={'height': '1%'}),
                dbc.Row(children=[
                    #dbc.Col(width=100, style={'width': '1%'}),
                    dbc.Row(
                        [
                            dbc.Card(
                                children=[
                                    dbc.CardHeader(
                                        "Plot 1",
                                        style={
                                            "background-color":"#e3e3e3",
                                            "font-weight":"bold",
                                            "font-size":"Large",
                                            "text-align": "center"
                                        }
                                    ),
                                    dbc.CardBody(
                                        _build_control_panel(1)
                                        + [
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            html.Div(["Tumor Information:"]),
                                            html.Div(
                                                [
                                                    html.Iframe(
                                                        srcDoc="", 
                                                        id='msg-box-dist-1', 
                                                        style={"width": "100%", "border": "0", "height": "80px"}
                                                    )
                                                ], 
                                                style={"margin": "3px", 'border-style': 'dotted', 'border-width': 'thin'}
                                            ),
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            dcc.Loading(
                                                id="loading-dist-1",
                                                type="default",
                                                color='black',
                                                children=html.Div([
                                                    html.Div(
                                                        id='dist-plot-1',
                                                        children=[
                                                            dcc.Graph(
                                                                figure=_build_dist_plot('PJ016', DEFAULT_GENE, 'all'),
                                                                config=_build_plot_config('svg'),
                                                            )
                                                        ]
                                                    )
                                                ], id='load-out-dist-1')
                                            ), 
                                            dbc.Row([
                                                dbc.Col([], width=100, style={"width": "15px"}),
                                                dbc.Col([
                                                    html.H6("Download format: ")
                                                ], width=100, style={"width": "40%"}),
                                                dbc.Col([
                                                    dcc.Dropdown(
                                                        options=[
                                                            {'label': 'SVG', 'value': 'svg'},
                                                            {'label': 'PNG', 'value': 'png'}
                                                        ],
                                                        value='svg',
                                                        id='img-format-dist-1'
                                                    )
                                                ])
                                            ])
                                        ]
                                    )
                                ]
                                #width={"size": "auto", "order": 1}
                            ),
                            dbc.Col([], width=100, style={"width": "15px"}),
                            dbc.Card(
                                children=[
                                    dbc.CardHeader(
                                        "Plot 2",
                                        style={
                                            "background-color":"#e3e3e3",
                                            "font-weight":"bold",
                                            "font-size":"Large",
                                            "text-align": "center"
                                        }
                                    ),
                                    dbc.CardBody(
                                        _build_control_panel(2)
                                        + [
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            html.Div(["Tumor Information:"]),
                                            html.Div(
                                                [
                                                    html.Iframe(
                                                        srcDoc="", 
                                                        id='msg-box-dist-2', 
                                                        style={"width": "100%", "border": "0", "height": "80px"}
                                                    )
                                                ], 
                                                style={"margin": "3px", 'border-style': 'dotted', 'border-width': 'thin'}
                                            ),
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            dcc.Loading(
                                                id="loading-dist-2",
                                                type="default",
                                                color='black',
                                                children=html.Div([
                                                    html.Div(
                                                        id='dist-plot-2',
                                                        children=[
                                                            dcc.Graph(
                                                                figure=_build_dist_plot('PJ016', DEFAULT_GENE, 'all'),
                                                                config=_build_plot_config('svg'),
                                                            )
                                                        ]
                                                    )
                                                ], id='load-out-dist-2')
                                            ),
                                            dbc.Row([
                                                dbc.Col([], width=100, style={"width": "15px"}),
                                                dbc.Col([
                                                    html.H6("Download format: ")
                                                ], width=100, style={"width": "40%"}),
                                                dbc.Col([
                                                    dcc.Dropdown(
                                                        options=[
                                                            {'label': 'SVG', 'value': 'svg'},
                                                            {'label': 'PNG', 'value': 'png'}
                                                        ],
                                                        value='svg',
                                                        id='img-format-dist-2'
                                                    )
                                                ])
                                            ])
                                        ]
                                    )
                                ]
                                #width={"size": "auto", "order": 1}
                            )
                        ], 
                        #style={"width": "80%"}
                    )
                ])
            ])
        ]
    )
    return layout
