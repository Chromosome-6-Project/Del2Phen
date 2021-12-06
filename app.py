# -*- coding: utf-8 -*-
"""Test Dash App."""

import pickle

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

import network

app = dash.Dash(__name__)

# df = pd.read_csv("https://plotly.github.io/datasets/country_indicators.csv")

with open("/home/tyler/Documents/Chr6_docs/Plots/nx_graph.pickle", "rb") as infile:
    nx_graph = pickle.load(infile)

# available_indicators = df["Indicator Name"].unique()

# =============================================================================
# app.layout = html.Div([
#     html.Div([
#
#         html.Div([
#             dcc.Dropdown(
#                 id="xaxis-column",
#                 options=[{"label": i, "value": i} for i in available_indicators],
#                 value="Fertility rate, total (births per woman)"
#             ),
#             dcc.RadioItems(
#                 id="xaxis-type",
#                 options=[{"label": i, "value": i} for i in ["Linear", "Log"]],
#                 value="Linear",
#                 labelStyle={"display": "inline-block"}
#             )
#         ], style={"width": "48%", "display": "inline-block"}),
#
#         html.Div([
#             dcc.Dropdown(
#                 id="yaxis-column",
#                 options=[{"label": i, "value": i} for i in available_indicators],
#                 value="Life expectancy at birth, total (years)"
#             ),
#             dcc.RadioItems(
#                 id="yaxis-type",
#                 options=[{"label": i, "value": i} for i in ["Linear", "Log"]],
#                 value="Linear",
#                 labelStyle={"display": "inline-block"}
#             )
#         ], style={"width": "48%", "float": "right", "display": "inline-block"})
#     ]),
#
#     dcc.Graph(id="indicator-graphic"),
#
#     dcc.Slider(
#         id="year--slider",
#         min=df["Year"].min(),
#         max=df["Year"].max(),
#         value=df["Year"].max(),
#         marks={str(year): str(year) for year in df["Year"].unique()},
#         step=None
#     )
# ])
# =============================================================================
class DefaultSlider(dcc.Slider):
    def __init__(self, id):
        super().__init__(
            id, min=0, max=100, value=0,
            marks={str(sim): str(sim) for sim in range(0, 101, 10)},
            tooltip={"placement": "bottom"}
            )

app.layout = html.Div([
    dcc.Graph(id="graph-with-slider"),
    html.Div(id="five-output"),
    DefaultSlider("size-slider"),
    DefaultSlider("overlap-slider"),
    DefaultSlider("gene-slider"),
    DefaultSlider("hi-slider"),
])


@app.callback(
    Output("graph-with-slider", "figure"),
    Output("five-output", "children"),
    Input("size-slider", "value"),
    Input("overlap-slider", "value"),
    Input("gene-slider", "value"),
    Input("hi-slider", "value"))
def update_figure(size_threshold, overlap_threshold, gene_threshold, hi_threshold):
    filtered_graph = network.filter_graph_edges(
        nx_graph,
        size_threshold/100,
        overlap_threshold/100,
        gene_threshold/100,
        hi_threshold/100
        )
    sizes = network.get_subnet_sizes(filtered_graph)
    five_count = str(sum((1 for size in sizes if size >= 5)))
    five_count = f"Clusters with at least 5 individuals: {five_count}"

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=sizes,
        xbins=dict(
            start=0,
            end=500,
            size=2),
        autobinx=False,
        name="",
        hovertemplate="Size Range: %{x}<br>Count: %{y}<extra></extra>",
        ))

    fig.update_layout(
        transition_duration=500,
        title_text="Distribution of Subgroup Sizes",
        xaxis_title_text="Subgroup Size",
        yaxis_title_text="Subgroup Count"
        )

    return fig, five_count


if __name__ == "__main__":
    app.run_server(debug=True)
