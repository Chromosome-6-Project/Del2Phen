
from dash import dcc, html
import dash_bootstrap_components as dbc


login_modal = dbc.Modal(id="login_modal",
                        is_open=False,
                        children=[
                            dbc.ModalHeader(children=[
                                dbc.ModalTitle("Login")
                                ]),
                            dbc.ModalBody(children=[
                                dbc.Label("Username:"),
                                dbc.Input(id="username_input", type="text", placeholder="Enter username"),
                                dbc.Label("Password:"),
                                dbc.Input(id="password_input", type="password", placeholder="Enter password",
                                          n_submit=0)
                                ]),
                            dbc.ModalFooter(children=[
                                dbc.Button("Login", id="login_button", className="ms-auto", n_clicks=0)
                                ])
                            ])

nav_links = [
    dbc.NavItem(children=[
        dbc.NavLink("Case Explorer", href="case-explorer")
        ]),
    dbc.NavItem(children=[
        dbc.NavLink("Historical Trends", href="historical-trends")
        ]),
    html.Div(id="login_div", style={"display": "block"}, children=[
        dbc.NavItem(children=[
            dbc.NavLink(id="open_login", children="Sign in", n_clicks=0, href="#!")
            ])
        ]),
    html.Div(id="logout_div", style={"display": "none"}, children=[
        dbc.NavItem(children=[
            dbc.NavLink(id="logout_button", children="Sign out", n_clicks=0, href="#!")
            ])
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

nav_collapse = dbc.Collapse(
            children=[dbc.Nav(
                children=nav_links,
                className="ms-auto",
                navbar=True
                )],
            id="navbar-collapse",
            navbar=True
            )

navbar = dbc.Navbar(
    children=[dbc.Container(
        children=[nav_title,
                  dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                  nav_collapse]
        )],
    color="#663399",
    # color="primary",
    dark=True,
    className="mb-5",
    )

memory = html.Div(children=[dcc.Store(id="user_credentials", data={}, storage_type="session"),
                            dcc.Store(id="historical_coverage", data={}),
                            dcc.Store(id="historical_insert_len", data={})])

head_layout = html.Div([memory,
                        dcc.Location(id="url", refresh=False),
                        navbar,
                        login_modal])

home_layout = html.Div(html.H2("Home"))
