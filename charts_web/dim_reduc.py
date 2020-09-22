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
    "#0652ff", #  electric blue
    "#e50000", #  red
    "#9a0eea", #  violet
    "#01b44c", #  shamrock
    "#fedf08", #  dandelion
    "#00ffff", #  cyan
    "#89fe05", #  lime green
    "#a2cffe", #  baby blue
    "#dbb40c", #  gold
    "#029386", #  teal
    "#ff9408", #  tangerine
    "#d8dcd6", #  light grey
    "#80f9ad", #  seafoam
    "#3d1c02", #  chocolate
    "#fffd74", #  butter yellow
    "#536267", #  gunmetal
    "#f6cefc", #  very light purple
    "#650021", #  maroon
    "#020035", #  midnight blue
    "#b0dd16", #  yellowish green
    "#9d7651", #  mocha
    "#c20078", #  magenta
    "#380282", #  indigo
    "#ff796c", #  salmon
    "#874c62", #  dark muave
    "#02ccfe", #  bright sky blue
    "#5fa052", #  muted green
    "#9a3001", #  auburn
    "#fc2647", #  pinky red
    "#d8863b", #  dull orange
    "#7b002c", #  bordeaux
    "#8e82fe", #  periwinkle
    "#ffff14", #  yellow
    "#ff073a", #  neon red
    "#6ecb3c", #  apple
    "#c45508", #  rust orange
    "#8756e4", #  purpley
    "#8756e4", #  diarrhea
    "#bcecac", #  light sage
    "#5d1451", #  grape purple
    "#028f1e", #  emerald green
    "#ffa62b", #  mango
    "#3a2efe", #  light royal blue
    "#c0022f", #  lipstick red
    "#0485d1", #  cerulean
    "#a57e52", #  puce
    "#380835", #  eggplant
    "#a9f971", #  spring green
    "#fe4b03", #  blood orange
    "#8cff9e", #  baby green
    "#86775f", #  brownish grey
    "#9d0759", #  dark fuchsia
    "#665fd1", #  dark periwinkle
    "#49759c", #  dullblue
    "#fffa86", #  manilla
    "#280137", #  midnight purple
    "#fa4224", #  orangey red
    "#d99b82", #  pinkish tan
    "#152eff", #  vivid blue
    "#f2ab15", #  squash
    "#70b23f", #  nasty green
    "#952e8f", #  warm purple
    "#bcf5a6", #  washed out green
    "#9a6200", #  raw sienna
    "#fb5ffc", #  violet pink
    "#ddd618", #  piss yellow
    "#fe420f", #  orangered
    "#c27e79", #  brownish pink
    "#adf802", #  lemon green
    "#29465b", #  dark grey blue
    "#48c072", #  dark mint
    "#edc8ff"  #  light lilac

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
    Output(component_id='dim-reduc-scatter-1', component_property='figure'),
    [   
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='dim-reduc-alg-1', component_property='value'),
        Input(component_id='num-dims-1', component_property='value'),
        Input(component_id='color-by-feature-1', component_property='value'),
        Input(component_id='select-feature-category-1', component_property='value'),
        Input(component_id='dot-size-1', component_property='value')
    ]
)
def update_dim_reduc_1(tumor, algo, num_dims, gene, category, dot_size):
    return _build_dim_reduc(tumor, algo, num_dims, gene, category, dot_size)


@app.callback(
    Output(component_id='color-by-feature-container-1', component_property='children'),
    [
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='select-feature-category-1', component_property='value')
    ]
)
def update_feature_category_selector_1(tumor, category):
    return build_features_selector('color-by-feature-1', tumor, category)

@app.callback(
    Output(component_id='msg-box-1', component_property='srcDoc'),
    [
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='select-feature-category-1', component_property='value')
    ]
)
def _update_msg(tum, feat):
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
    Output(component_id='msg-box-2', component_property='srcDoc'),
    [
        Input(component_id='select-tumor-2', component_property='value'),
        Input(component_id='select-feature-category-2', component_property='value')
    ]
)
def _update_msg(tum, feat):
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
    Output(component_id='dim-reduc-scatter-2', component_property='figure'),
    [
        Input(component_id='select-tumor-2', component_property='value'),
        Input(component_id='dim-reduc-alg-2', component_property='value'),
        Input(component_id='num-dims-2', component_property='value'),
        Input(component_id='color-by-feature-2', component_property='value'),
        Input(component_id='select-feature-category-2', component_property='value'),
        Input(component_id='dot-size-2', component_property='value')
    ]
)
def update_dim_reduc_2(tumor, algo, num_dims, gene, category, dot_size):
    return _build_dim_reduc(tumor, algo, num_dims, gene, category, dot_size)


@app.callback(
    Output(component_id='color-by-feature-container-2', component_property='children'),
    [
        Input(component_id='select-tumor-2', component_property='value'),
        Input(component_id='select-feature-category-2', component_property='value')
    ]
)
def update_feature_category_selector_2(tumor, category):
    return build_features_selector('color-by-feature-2', tumor, category)


def build_dim_reduc_selector(idd):
    return dcc.Dropdown(
        options=[
            {'label': 'UMAP', 'value': 'umap'},
            {'label': 'PHATE', 'value': 'phate'}
        ],
        value=DEFAULT_ALGO,
        id=idd
    )

def build_num_dims_selector(idd):
    return dcc.Dropdown(
        options=[
            {'label': '2', 'value': 2},
            {'label': '3', 'value': 3}
        ],
        value=DEFAULT_NUM_DIM,
        id=idd
    )

def build_feature_category_selector(idd):
    return dcc.Dropdown(
        options=[
            {'label': 'Gene', 'value': 'gene'},
            {'label': 'Cluster', 'value': 'cluster'},
            {'label': 'Cell Type Probability', 'value': 'cell_type_probability'},
            {'label': 'Predicted Cell Type', 'value': 'cell_type_classification'},
            {'label': 'Malignancy Score', 'value': 'malignancy_score'},
            {'label': 'CancerSEA Gene Set Score', 'value': 'cancersea_enrichment'},
            {'label': 'Hallmark Gene Set Score', 'value': 'hallmark_enrichment'}
            #{'label': 'Tumor', 'value': 'tumor'}
        ],
        value='gene',
        id=idd
    )


def build_features_selector(idd, tumor, category):
    if category == 'gene':
        return dcc.Input(id=idd, value=DEFAULT_GENE)
    elif category == 'cluster':
        return dcc.Dropdown(
            options=[{'label': 'Cluster', 'value': 'Cluster'}],
            value='Cluster',
            id=idd
        )
    elif category == 'cell_type_probability':
        if '&' in tumor:
            tum_1 = tumor.split('&')[0]
            tum_2 = tumor.split('&')[1]
            return common.build_cell_type_probability_dropdown_mult_tumors(
                tum_1,
                tum_2,
                idd     
            )
        else:
            return common.build_cell_type_probability_dropdown(
                tumor, 
                idd        
            )
    elif category == 'tumor':
        return dcc.Dropdown(
            options=[{'label': 'Tumor', 'value': 'Tumor'}],
            value='Tumor',
            id=idd
        )
    elif category == 'cell_type_classification': 
        return dcc.Dropdown(
            options=[{'label': 'Cell Type', 'value': 'Cell Type'}],
            value='Cell Type',
            id=idd
        )
    elif category == 'hallmark_enrichment':
        if '&' in tumor:
            tum_1 = tumor.split('&')[0]
            tum_2 = tumor.split('&')[1]
            return common.build_hallmark_enrichment_dropdown_mult_tumors(
                tum_1,
                tum_2,
                idd
            )
        else:
            return common.build_hallmark_enrichment_dropdown(
                tumor,
                idd
            )
    elif category == 'cancersea_enrichment':
        if '&' in tumor:
            tum_1 = tumor.split('&')[0]
            tum_2 = tumor.split('&')[1]
            return common.build_cancersea_enrichment_dropdown_mult_tumors(
                tum_1,
                tum_2,
                idd
            )
        else:
            return common.build_cancersea_enrichment_dropdown(
                tumor,
                idd
            )
    elif category == 'malignancy_score':
        return dcc.Dropdown(
            options=[{'label': 'Malignancy Score', 'value': 'Malignancy Score'}],
            value='Malignancy Score',
            id=idd
        )


def _build_dim_reduc(tumor_id, algo, num_dims, feat, category, dot_size):
    if algo == 'umap':
        df_dim_reduc = load_data.load_tumor_umap(tumor_id, num_dims)
    elif algo == 'phate':
        df_dim_reduc = load_data.load_tumor_phate(tumor_id, num_dims)

    if '&' in tumor_id:
        tums = tumor_id.split('&')
        tum_1 = tums[0]
        tum_2 = tums[1]
        hover_texts = load_data.load_hover_texts_mult_tumors( 
            tum_1,
            tum_2,
            df_dim_reduc.index
        )
    else:
        hover_texts = load_data.hover_texts(tumor_id, df_dim_reduc.index)

    if category in set(['gene', 'cell_type_probability', 'hallmark_enrichment', 'cancersea_enrichment', 'malignancy_score']):
        # Check if this is an aligned set of tumors
        if '&' in tumor_id:
            tums = tumor_id.split('&')
            tum_1 = tums[0]
            tum_2 = tums[1]
            if category == 'gene':
                df_color = load_data.load_color_by_real_value_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index,
                    'log1_tpm',
                    'gene_name',
                    feat
                )
            elif category == 'cell_type_probability':
                df_color = load_data.load_color_by_real_value_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index,
                    'cell_type_probability',
                    'cell_type_probability_columns',
                    feat
                )
            elif category == 'hallmark_enrichment':
                df_color = load_data.load_color_by_real_value_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index,
                    'hallmark_gsva',
                    'hallmark_gene_set_name',
                    feat
                )
            elif category == 'cancersea_enrichment':
                df_color = load_data.load_color_by_real_value_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index,
                    'cancersea_gsva',
                    'cancersea_gene_set_name',
                    feat
                )
            elif category == 'malignancy_score':
                df_color = load_data.load_malignancy_score_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index
                )
        else:
            if category == 'gene':
                df_color = load_data.load_tumor_gene(tumor_id, feat)
            elif category == 'cell_type_probability':
                df_color = load_data.load_tumor_cell_type_probabilities(
                    tumor_id, 
                    feat
                )
            elif category == 'hallmark_enrichment':
                df_color = load_data.load_tumor_hallmark_enrichment(
                    tumor_id, 
                    feat
                )
            elif category == 'cancersea_enrichment':
                df_color = load_data.load_tumor_cancersea_enrichment(
                    tumor_id,
                    feat
                )
            elif category == 'malignancy_score':
                df_color = load_data.load_malignancy_score(tumor_id)

        col = 'color_by'

        # Determine color range
        if category == 'cell_type_probability':
            cmin = 0.0
            cmax = 1.0
        elif category == 'hallmark_enrichment' or category == 'cancersea_enrichment':
            end = max(
                abs(min(df_color[col])), 
                abs(max(df_color[col]))
            )
            cmin = -1 * end
            cmax = end
        else:
            color_range = [
                min(df_color[col]),
                max(df_color[col])
            ]
            cmin = color_range[0]
            cmax = color_range[1]

        # Determine color map
        if category == 'hallmark_enrichment' or category == 'cancersea_enrichment':
            #palette = 'RdBu'
            palette = 'RdYlBu_r'
        else:
            palette = 'Viridis'

        # Create the full data frame
        df = df_dim_reduc.join(df_color)
        
        # Create figures
        if num_dims == 3:
            markers=dict(
                size=dot_size,
                color=df[col],
                colorscale=palette,
                opacity=0.0,
                cmin=cmin,
                cmax=cmax,
                colorbar=dict(
                    thickness=20
                )
            )
            fig = go.Figure(data=[go.Scatter3d(
                x=df[df_dim_reduc.columns[0]],
                y=df[df_dim_reduc.columns[1]],
                z=df[df_dim_reduc.columns[2]],
                mode='markers',
                marker=markers,
                showlegend=False,
                hovertemplate="%{hovertext}<extra></extra>",
                hovertext=hover_texts
            )])
        elif num_dims == 2:
            if len(df) > 5000:
                size=2.5
            else:
                size=5
            markers=dict(
                size=dot_size,
                color=df[col],
                colorscale=palette,
                opacity=1.0,
                cmin=cmin,
                cmax=cmax,
                colorbar=dict(
                    thickness=20
                )
            )
            fig = go.Figure(data=[go.Scatter(
                x=df[df_dim_reduc.columns[0]],
                y=df[df_dim_reduc.columns[1]],
                mode='markers',
                marker=markers,
                showlegend=False,
                hovertemplate="%{hovertext}<extra></extra>",
                hovertext=hover_texts
            )])
    elif category == 'cluster' or category == 'tumor' or category == 'cell_type_classification':
        if category == 'cluster':
            if '&' in tumor_id:
                tums = tumor_id.split('&')
                tum_1 = tums[0]
                tum_2 = tums[1]
                df_color = load_data.load_clusters_mult_tumors(
                    tum_1,
                    tum_2,
                    df_dim_reduc.index
                )
            else:
                df_color = load_data.load_tumor_clusters_for_cells(tumor_id)
        elif category == 'tumor':
            if '&' in tumor_id:
                tums = tumor_id.split('&')
                tum_1 = tums[0]
                tum_2 = tums[1]
                df_color = load_data.load_tumors_for_cells_mult_tumors(
                    tum_1, 
                    tum_2, 
                    df_dim_reduc.index
                )
            else:
                df_color = pd.DataFrame(
                    data={
                        'color_by': [tumor_id for x in df_dim_reduc.index]
                    },
                    index=df_dim_reduc.index
                )
        elif category == 'cell_type_classification':
            if '&' in tumor_id:
                tums = tumor_id.split('&')
                tum_1 = tums[0]
                tum_2 = tums[1]
                df_color = load_data.load_cell_type_classifications_mult_tumors(
                    tum_1,
                    tum_2, 
                    df_dim_reduc.index
                )
            else:
                df_color = load_data.load_tumor_cell_type_classifications(tumor_id)

        col = 'color_by'
        df = df_dim_reduc.join(df_color)
        # Convert hover-texts to a dataframe so that we can plot subsets at a time
        hover_df = pd.DataFrame(
            data={'hover_text': hover_texts},
            index=df.index
        )
        fig = go.Figure()
        for clust_i, group in enumerate(sorted(set(df[col]))):
            df_clust = df.loc[df[col] == group]
            if len(df) > 5000:
                size=1.5
            else:
                size=5
            markers=dict(
                size=dot_size,
                color=PALETTE[clust_i],
                opacity=1.0
            )
            if category == 'tumor' or category == 'cell_type_classification':
                group_name = group
            elif category == 'cluster':
                group_name = "Cluster {}".format(group)
            group_name
            if num_dims == 3:
                fig.add_trace(
                    go.Scatter3d(
                        x=df_clust[df_clust.columns[0]],
                        y=df_clust[df_clust.columns[1]],
                        z=df_clust[df_clust.columns[2]],
                        mode='markers',
                        marker=markers,
                        name=group_name,
                        hovertemplate="%{hovertext}<extra></extra>",
                        hovertext=list(hover_df.loc[df_clust.index]['hover_text'])
                    )
                )
            elif num_dims == 2:
                fig.add_trace(
                    go.Scatter(
                        x=df_clust[df_clust.columns[0]],
                        y=df_clust[df_clust.columns[1]],
                        mode='markers',
                        marker=markers,
                        name=group_name,
                        hovertemplate="%{hovertext}<extra></extra>",
                        hovertext=list(hover_df.loc[df_clust.index]['hover_text'])
                    )
                )

    # Avoid an exception when trying to set the z-axis title
    if num_dims == 2:
        z_axis_title = None
    elif num_dims == 3:
        z_axis_title = df_dim_reduc.columns[2]

    fig.update_layout(
        autosize=True,
        width=FIG_DIM+10,
        height=FIG_DIM,
        scene = dict(
            xaxis = dict(
                showticklabels=False,
                title=df_dim_reduc.columns[0]
            ),
            yaxis = dict(
                showticklabels=False,
                title=df_dim_reduc.columns[1]
            ),
            zaxis = dict(
                showticklabels=False,
                title=z_axis_title
            )
        ),
        xaxis=dict(
            showticklabels=False
        ),
        yaxis=dict(
            showticklabels=False
        ),
        xaxis_title=df_dim_reduc.columns[0],
        yaxis_title=df_dim_reduc.columns[1],
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
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
                            common.build_tumor_dropdown('select-tumor-{}'.format(plot_num))
                        ])
                    ])
                ]),
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Select algorithm:")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            build_dim_reduc_selector('dim-reduc-alg-{}'.format(plot_num))
                        ])
                    ])
                ]),
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Select dimensions:")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            build_num_dims_selector('num-dims-{}'.format(plot_num))
                        ])
                    ])
                ]),
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Color by:")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            build_feature_category_selector('select-feature-category-{}'.format(plot_num))
                        ])
                    ])
                ]),
                html.Div([
                    dbc.Row(children=[
                        dbc.Col([
                            html.H6("Enter feature: ")
                        ], width=100, style={"width": "40%"}),
                        dbc.Col([
                            html.Div([
                                build_features_selector('color-by-feature-{}'.format(plot_num), 'PJ016', 'gene')
                            ], id='color-by-feature-container-{}'.format(plot_num))
                        ])
                    ])
                ])
            ])
        ])
    ]

@app.callback(Output("louad-out-1", "children"), [Input("loading-input-1", "value")])
def input_triggers_spinner(value):
    time.sleep(1)
    return value

@app.callback(Output("louad-out-2", "children"), [Input("loading-input-2", "value")])
def input_triggers_spinner(value):
    time.sleep(1)
    return value

def build_layout():
    layout = dcc.Tab(
        label='Dimension Reduction',
        children=[
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.Div(['Compare tumors and subpopulations via dimension reduction scatterplots. For each plot, select a tumor as well as an attribute used to color each tumor.']),
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
                                                        id='msg-box-1', 
                                                        style={"width": "100%", "border": "0", "height": "80px"}
                                                    )
                                                ], 
                                                style={"margin": "3px", 'border-style': 'dotted', 'border-width': 'thin'}
                                            ),
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            dcc.Loading(
                                                id="loading-1",
                                                type="default",
                                                color='black',
                                                children=html.Div([
                                                    dcc.Graph(
                                                        id='dim-reduc-scatter-1',
                                                        figure=_build_dim_reduc('PJ016', DEFAULT_ALGO, DEFAULT_NUM_DIM, DEFAULT_GENE, 'gene', 2),
                                                        config={
                                                            'displayModeBar': True,
                                                            'toImageButtonOptions': {
                                                                'format':'svg',
                                                                'filename': 'dash_plot'
                                                            },
                                                            "displaylogo": False
                                                        }
                                                    )
                                                ], id='louad-out-1')
                                            ), 
                                            dbc.Row([
                                                dbc.Col([], width=100, style={"width": "15px"}),
                                                dbc.Col([
                                                    html.H6("Dot size: ")
                                                ], width=100, style={"width": "40%"}),
                                                dbc.Col([
                                                    dcc.Slider(
                                                        min=1,
                                                        max=6,
                                                        #step=None,
                                                        value=2,
                                                        id='dot-size-1'
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
                                                        id='msg-box-2', 
                                                        style={"width": "100%", "border": "0", "height": "80px"}
                                                    )
                                                ], 
                                                style={"margin": "3px", 'border-style': 'dotted', 'border-width': 'thin'}
                                            ),
                                            dbc.Row(html.Hr(), style={'height': '3%'}),
                                            dcc.Loading(
                                                id="loading-2",
                                                type="default",
                                                color='black',
                                                children=html.Div([
                                                    dcc.Graph(
                                                        id='dim-reduc-scatter-2',
                                                        figure=_build_dim_reduc('PJ016', DEFAULT_ALGO, DEFAULT_NUM_DIM, DEFAULT_GENE, 'gene', 2),
                                                        config={
                                                            'displayModeBar': True,
                                                            'toImageButtonOptions': {
                                                                'format':'svg',
                                                                'filename': 'dash_plot'
                                                            },
                                                            "displaylogo": False
                                                        }
                                                    )
                                                ], id='louad-out-2')
                                            ),
                                            dbc.Row([
                                                dbc.Col([], width=100, style={"width": "15px"}),
                                                dbc.Col([
                                                    html.H6("Dot size: ")
                                                ], width=100, style={"width": "40%"}),
                                                dbc.Col([
                                                    dcc.Slider(
                                                        min=1,
                                                        max=6,
                                                        #step=None,
                                                        value=2,
                                                        id='dot-size-2'
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
