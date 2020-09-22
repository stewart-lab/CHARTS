import dash_bootstrap_components as dbc
import dash_html_components as html


CHARTS_LOGO = "https://github.com/mbernste/cancer-single-cell-biomarker/raw/master/img/charts_logo.png"

LAYOUT = dbc.NavbarSimple(
    children=[
      dbc.NavItem(dbc.NavLink("About", href="/about")),
      dbc.DropdownMenu(
         nav=True,
         in_navbar=True,
         label="Menu",
         children=[
            dbc.DropdownMenuItem("FAQ",  href="/faq"),
            dbc.DropdownMenuItem("Datasets",  href="/data_set_summary")
                  ],
              ),
            ],
    brand="CHARTS: Characterizing Tumor Subpopulations",
    brand_style={"font-size": "200%"},
    brand_href="/charts",
    sticky="top",
)

