import os
import dash_html_components as html
import dash_bootstrap_components as dbc
from flask import send_from_directory
import app
from urllib.parse import quote as urlquote
import json

import nav
import footer


with open('config.json', 'r') as f:
    config = json.load(f)
    downloads_dir = config['downloads_dir']


@app.server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(downloads_dir, path, as_attachment=True)


def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)


def get_files():
    files = []
    for filename in os.listdir(downloads_dir):
        path = os.path.join(downloads_dir, filename)
        if os.path.isfile(path):
            files.append(filename)
    return sorted(files)


def build_layout():
    return html.Div(
        [
            nav.LAYOUT,
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.H2(children='Downloads', style={"text-align": "center"}),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.Div(['Download tumor datasets as tab-separated files.']),
            dbc.Row(html.Hr(), style={'height': '3%'}),
            html.Ul(
                id="file-list", 
                children=[html.Li(file_download_link(filename)) for filename in get_files()]
            ),
            footer.LAYOUT
        ],
        style={"padding-left": "10%", "padding-right": "10%", "width": "100%"}
    )

