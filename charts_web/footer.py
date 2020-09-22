import dash_bootstrap_components as dbc
import dash_html_components as html


MORGRIDGE_LOGO = "https://mk0morgridgeorgkfioq.kinstacdn.com/wp-content/uploads/MorgridgeLogo_vertical_no-artifact.png"
UW_LOGO = "https://raw.githubusercontent.com/deweylab/MetaSRA-website-frontend/master/src/assets/uw-logo.svg"

LAYOUT = html.Div(
    children=[
        html.Span([
            html.A(
                [
                    html.Img(src=MORGRIDGE_LOGO, style={'width':'200px'}),
                ], 
                href='https://morgridge.org',
                target='_blank'
            ),
            html.A(
                [
                    html.Img(src=UW_LOGO, style={'width':'70px'}),
                ],
                href='https://www.wisc.edu',
                target='_blank',
                style={'margin-left': '30px'}
            ),
            html.P(u"\u00A9" + " 2020", style={'float': 'right', 'margin-right': '100px'})
        ])
    ]
)

