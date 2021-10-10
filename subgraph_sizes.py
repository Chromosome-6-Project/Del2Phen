# -*- coding: utf-8 -*-
"""Chromosome 6 Project - Filter Effects on Subgraph Sizes.

Created on Fri Sep 17 11:54:00 2021

@author: Ty
"""

import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"


def make_nx_nodes(nodes):
    nodes = [(node.id, node._asdict()) for node in nodes.values()]
    return nodes


def make_nx_edges(edges):
    edges = [(edge.id1, edge.id2, edge._asdict()) for edge in edges]
    return edges


def make_nx_graph(nodes, edges):
    nodes = make_nx_nodes(nodes)
    edges = make_nx_edges(edges)
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    return graph


def filter_graph_edges(graph, overlap_threshold=0, gene_threshold=0, hi_threshold=0):
    edges = list(graph.edges.items())
    edges = [(edge[0][0], edge[0][1], edge[1]) for edge in edges
             if edge[1]["sim"] >= gene_threshold]
    filtered_graph = nx.Graph()
    filtered_graph.add_nodes_from(graph.nodes)
    filtered_graph.add_edges_from(edges)
    return filtered_graph


def get_subnet_sizes(graph):
    components = nx.connected_components(graph)
    sizes = sorted([len(component) for component in components])
    return sizes


# Create figure
def plot_histograms_with_slider(graph):
    fig = go.Figure()

    for threshold in np.arange(0, 1.05, .05):

        sizes = get_subnet_sizes(filter_graph_edges(graph, gene_threshold=threshold))

        fig.add_trace(
            go.Histogram(
                x=sizes,
                visible=False,
                xbins=dict(
                    start=0,
                    end=500,
                    size=2),
                autobinx=False,
                alignmentgroup=0,
                name="",
                hovertemplate="Size Range: %{x}<br>Count: %{y}<extra></extra>",
                )
            )

    # Make 10th trace visible
    fig.data[0].visible = True

    # Create and add slider
    steps = []
    for i in range(len(fig.data)):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(fig.data)},
                  {"title": f"Number of groups of size n with gene similarity threshold of {i*.05:.0%}"}],  # layout attribute
            label=f"{i*.05:.0%}"
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Threshold: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        sliders=sliders,
        bargap=0.1
    )

    fig.show()
