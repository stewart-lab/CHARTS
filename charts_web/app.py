import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

external_scripts = [
    {'src': 'https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js'},
    {'src': 'https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/js/bootstrap.min.js'}
]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css', './custom.css']
#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
print(dbc.themes.BOOTSTRAP)
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP, {'src': './custom.css'}], external_scripts=external_scripts, suppress_callback_exceptions=True)
app.title = 'CHARTS'
