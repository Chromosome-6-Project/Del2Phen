# -*- coding: utf-8 -*-
"""Test Dash App."""

from math import log10
import pickle

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import chr6
import network


_, _, _, comparison, _, ontology, _ = chr6.test(
    genotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_array_2022-02-21_18_28_06.csv",
    phenotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_questionnaire_2022-02-21_10_18_09.csv",
    patient_hpo="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_research_patients_2022-02-21_10_20_53.csv",
    drop_list_file="/home/tyler/Documents/Chr6_docs/PatientData/drop_list.txt",
    expand_hpos=True
    )

with open("/home/tyler/Documents/Chr6_docs/Phenotype_Homogeneity/selected_phenotypes_hpos.txt") as infile:
    selected_hpos = infile.readlines()
selected_hpos = [ontology[x.strip().split("\t")[1]] for x in selected_hpos]

nx_graph = network.make_nx_graph_from_comparisons(comparison)


class DefaultSlider(dcc.Slider):
    def __init__(self, id):
        super().__init__(
            id, min=0, max=100, value=0,
            marks={str(sim): str(sim) for sim in range(0, 101, 10)},
            tooltip={"placement": "bottom"}
            )

app = dash.Dash(__name__)

app.layout = html.Div([
    # Panes
    html.Div([
        # Left Pane
        html.Div([
            dcc.Graph(id="graph-with-slider"),
            html.Div(id="subsize-5-count"),
            html.Div(id="link-5-count"),
            html.Hr(),
            html.H2("Similarity Settings"),
            html.Br(),
            html.Label("Size Similarity"),
            DefaultSlider("length-slider"),
            html.Br(),
            html.Label("Overlap"),
            DefaultSlider("overlap-slider"),
            html.Br(),
            html.Label("Gene Similarity"),
            DefaultSlider("gene-slider"),
            html.Br(),
            html.Label("HI Gene Similarity"),
            DefaultSlider("hi-slider"),
            ], style={"padding": 10, "flex": 1}),
        # Right Pane
        html.Div([
            dcc.Graph(id="homogeneity-graph"),
            html.Label("Prevalence Level Threshold"),
            dcc.Slider("homogeneity-slider",
                       min=0, max=100, value=20,
                       marks={str(sim): str(sim) for sim in range(0, 101, 10)},
                       tooltip={"placement": "bottom"}),
            html.Br(),
            html.Label("Group Size Threshold"),
            dcc.Slider("group-size-slider",
                       min=0, max=100, value=5,
                       marks={str(sim): str(sim) for sim in range(0, 101, 10)},
                       tooltip={"placement": "bottom"}),
            ], style={"padding": 10, "flex": 1}),
        ], style={"display": "flex", "flex-direction": "row"}),
])

@app.callback(Output("homogeneity-graph", "figure"),
              Input("homogeneity-slider", "value"),
              Input("group-size-slider", "value"),
              Input("length-slider", "value"),
              Input("overlap-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-slider", "value"))
def update_homogeneity_figure(homogeneity_threshold, group_size_threshold,
                              length_threshold, overlap_threshold,
                              gene_threshold, hi_threshold):
    _, upper_homogens, lower_homogens = comparison.test_all_homogeneities(
        selected_hpos, length_threshold/100, overlap_threshold/100,
        gene_threshold/100, hi_threshold/100, hpo_similarity=0,
        group_size_threshold=group_size_threshold
        )
    upper_homogens = [100*homogen.calculate_homogeneity(homogeneity_threshold/100)
                      for homogen in upper_homogens.values()]
    lower_homogens = [100*homogen.calculate_homogeneity(homogeneity_threshold/100)
                      for homogen in lower_homogens.values()]
    fig = go.Figure()
    fig.add_trace(
        trace=go.Histogram(
            x=upper_homogens,
            name=f"Group size >= {group_size_threshold}",
            xbins=dict(start=0,
                       end=105,
                       size=5),
            autobinx=False,
            hovertemplate="Homogeneity Score: %{x}%<br>Count: %{y}<extra></extra>"
            )
        )
    fig.add_trace(
        trace=go.Histogram(
            x=lower_homogens,
            name=f"Group size < {group_size_threshold}",
            xbins=dict(start=0,
                       end=105,
                       size=5),
            autobinx=False,
            hovertemplate="Homogeneity Score: %{x}%<br>Count: %{y}<extra></extra>"
            )
        )
    fig.update_layout(
        transition_duration=500,
        title_text="Distribution of Homogeneity Scores",
        xaxis_title_text="Score",
        yaxis_title_text="Count",
        barmode="stack"
        )
    return fig

@app.callback(Output("graph-with-slider", "figure"),
              Output("subsize-5-count", "children"),
              Output("link-5-count", "children"),
              Input("length-slider", "value"),
              Input("overlap-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-slider", "value"))
def update_figure(length_threshold, overlap_threshold, gene_threshold, hi_threshold):
    filtered_graph = network.filter_graph_edges(nx_graph, length_threshold/100,
                                                overlap_threshold/100, gene_threshold/100,
                                                hi_threshold/100)
    subsizes = network.get_subnet_sizes(filtered_graph)
    subsize_gt_5 = network.count_subnets_over_size_n(filtered_graph, 5)
    subsize_gt_5 = f"Clusters with at least 5 individuals: {subsize_gt_5}"

    links = sorted(list(network.get_node_degrees(filtered_graph).values()))
    links_gt_5 = network.count_nodes_over_n_degree(filtered_graph, 5)
    links_gt_5 = f"Individuals with at least 5 links: {links_gt_5}"

    # fig1 = make_subplots(rows=1, cols=2)
    fig1 = go.Figure()
    fig1.add_trace(
        # row=1,
        # col=1,
        trace=go.Histogram(
            x=subsizes,
            xbins=dict(start=0,
                       end=500,
                       size=1),
            autobinx=False,
            name="Subnet Sizes",
            hovertemplate="Subnet Size: %{x}<br>Count: %{y}<extra></extra>"
            )
        )
    fig1.add_trace(
        # row=1,
        # col=2,
        trace=go.Histogram(
            x=links,
            xbins=dict(start=0,
                       end=500,
                       size=1),
            autobinx=False,
            name="Node Degrees",
            hovertemplate="Links: %{x}<br>Count: %{y}<extra></extra>"
            )
        )
    fig1.update_layout(
        transition_duration=500,
        title_text="Distribution of Subnet Sizes and Node Degrees",
        xaxis_title_text="Size",
        yaxis_title_text="Count"
        )
    fig1.update_yaxes(type="log", dtick=log10(2))

    return fig1, subsize_gt_5, links_gt_5


if __name__ == "__main__":
    app.run_server(debug=True, use_reloader=False)
