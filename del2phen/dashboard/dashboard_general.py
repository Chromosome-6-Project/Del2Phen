
from dash import dcc, html
import dash_bootstrap_components as dbc


nav_links = [
    dbc.NavItem(children=[
        dbc.NavLink(id="option_pane", children="Options")
        ])
    ]

nav_title = dcc.Link(href="home",
                     style={"textDecoration": "none"},
                     children=[
                         dbc.Row(align="center",
                                 className="g-0",
                                 children=[
                                     dbc.Col(children=[
                                         html.Img(src="assets/c6_logo_white.png", height="50px")
                                         ]),
                                     dbc.Col(children=[
                                         dbc.NavbarBrand("The Chromosome 6 Project", className="ms-2")
                                         ])
                                     ])
                         ])

nav_collapse = dbc.Collapse(id="navbar-collapse", navbar=True, children=[
    dbc.Nav(children=nav_links,
            className="ms-auto",
            navbar=True)
    ])

navbar = dbc.Navbar(children=[
    dbc.Container(children=[
        nav_title,
        dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
        nav_collapse
        ])
    ],
    color="#663399",
    # color="primary",
    dark=True,
    className="mb-5",
    )


head_layout = html.Div([navbar])
