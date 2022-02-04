# -*- coding: utf-8 -*-
"""Test Dash App."""

import pickle
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import network

app = dash.Dash(__name__)

with open("/home/tyler/Documents/Chr6_docs/Plots/nx_graph.pickle", "rb") as infile:
    nx_graph = pickle.load(infile)


class DefaultSlider(dcc.Slider):
    def __init__(self, id):
        super().__init__(
            id, min=0, max=100, value=0,
            marks={str(sim): str(sim) for sim in range(0, 101, 10)},
            tooltip={"placement": "bottom"}
            )

app.layout = html.Div([
    dcc.Graph(id="graph-with-slider"),
    html.Div(id="subsize-5-count"),
    html.Div(id="link-5-count"),
    DefaultSlider("size-slider"),
    DefaultSlider("overlap-slider"),
    DefaultSlider("gene-slider"),
    DefaultSlider("hi-slider"),
])


@app.callback(
    Output("graph-with-slider", "figure"),
    Output("subsize-5-count", "children"),
    Output("link-5-count", "children"),
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
    subsizes = network.get_subnet_sizes(filtered_graph)
    subsize_gt_5 = network.count_subnets_over_size_n(filtered_graph, 5)
    subsize_gt_5 = f"Clusters with at least 5 individuals: {subsize_gt_5}"

    links = sorted(list(network.get_node_degrees(filtered_graph).values()))
    links_gt_5 = network.count_nodes_over_n_degree(filtered_graph, 5)
    links_gt_5 = f"Individuals with at least 5 links: {links_gt_5}"

    fig = make_subplots(rows=1, cols=2)
    fig.add_trace(
        row=1,
        col=1,
        trace=go.Histogram(
            x=subsizes,
            xbins=dict(
            start=0,
            end=500,
            size=1),
            autobinx=False,
            name="",
            hovertemplate="Subnet Size: %{x}<br>Count: %{y}<extra></extra>"
            )
        )

    fig.add_trace(
        row=1,
        col=2,
        trace=go.Histogram(
            x=links,
            xbins=dict(
            start=0,
            end=500,
            size=1),
            autobinx=False,
            name="",
            hovertemplate="Links: %{x}<br>Count: %{y}<extra></extra>"
            )
        )

    fig.update_layout(
        transition_duration=500,
        title_text="Distribution of Subnet Sizes and Node Degrees",
        xaxis_title_text="Size",
        yaxis_title_text="Count"
        )

    return fig, subsize_gt_5, links_gt_5


if __name__ == "__main__":
    app.run_server(debug=True)
