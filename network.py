#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Network analysis.

@author: T.D. Medina
"""

from collections import namedtuple
from datetime import datetime
from math import log
import os

import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

# %% Nodes and Edges
Node = namedtuple("Node", ["id", "label", "group", "color",
                           "value", "title", "HI", "ranges"])
Edge = namedtuple("Edge", ["edge_id", "id1", "id2", "width", "color", "title",
                           "length_sim", "overlap_sim", "gene_sim", "hi_sim"])

def build_network_nodes(comparison_table):
    """Build patient node objects from comparison table."""
    colors = {"Literature case report": "blue",
              "Parental uploaded array report": "pink",
              None: "black",
              "HI Gene": "yellow"}

    nodes = {}
    for patient in comparison_table.index:
        lookup = comparison_table.lookup(patient)
        group = comparison_table.patient_db[patient].origin
        ranges = []
        for cnv in comparison_table.patient_db[patient].cnvs:
            ranges.append(f"{cnv.chromosome}:"
                          f"{cnv.range.start}:"
                          f"{cnv.range.stop}")
        ranges = ";".join(ranges)
        color = colors[group]
        if group == "HI Gene":
            hi = comparison_table.patient_db[patient].score
            value = 2
            title = (f"<p>{patient}<br>"
                     f"Group: {group}<br>"
                     f"HI score: {hi}<br>"
                     f"Genes: {value}<br></p>")
        else:
            hi = len(comparison_table.patient_db[patient].all_HI_genes())
            value = lookup.gene_count
            title = (f"<p>{patient}<br>"
                     f"Group: {group}<br>"
                     f"Affected genes: {value}<br>"
                     f"Affected HI genes: {hi}<br>"
                     f"HPO terms: {lookup.hpo_count}</p>")
        node = Node(patient, patient, group, color, value, title, hi, ranges)
        nodes[patient] = node
    return nodes


def write_network_nodes(nodes, out, normalize=True):
    """Write patient nodes to CSV file."""
    writer = ["id,label,group,color,value,title,hi,ranges\n"]
    for node in sorted(list(nodes.values()), key=lambda x: x[0]):
        if normalize:
            node = list(node)
            node[4] = log(node[4] + 1, 2)
        writer.append(",".join([str(x) for x in node]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)


def build_network_edges(comparison_table):
    """Build patient-patient edge objects from comparison table."""
    edges = []
    for intersect in comparison_table:
        if not intersect.gene_count or intersect.patients[0] == intersect.patients[1]:
            continue

        id1 = intersect.patients[0].id
        id2 = intersect.patients[1].id
        edge_id = f"{id1}_{id2}"

        length_sim = intersect.length_similarity

        overlap = intersect.loci_shared_size
        overlap_sim = intersect.loci_similarity

        gene_count = intersect.gene_count
        gene_sim = intersect.gene_similarity

        hpo_count = intersect.hpo_count
        hpo_sim = intersect.hpo_similarity

        hi_gene_count = intersect.hi_gene_count
        hi_gene_sim = intersect.hi_gene_similarity

        color = "gray"
        title = (f"<p>{id1}---{id2}:<br>"
                 f"Overlap: {overlap} ({overlap_sim:.2%})<br>"
                 f"Shared genes: {gene_count} ({gene_sim:.2%})<br>"
                 f"Shared HI genes: {hi_gene_count} ({hi_gene_sim:.2%})<br>"
                 f"Shared HPO terms: {hpo_count} ({hpo_sim:.2%})<br></p>")
        edges.append(Edge(edge_id, id1, id2, gene_count, color, title,
                          length_sim, overlap_sim, gene_sim, hi_gene_sim))
    return edges


def write_network_edges(edges, out, normalize=True):
    """Write patient-patient edges to CSV file."""
    writer = ["id,from,to,width,color,title,length_sim,overlap_sim,gene_sim,hi_gene_sim\n"]
    for edge in edges:
        if normalize:
            edge = list(edge)
            edge[3] = log(edge[3], 2)
        writer.append(",".join([str(x) for x in edge]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)


def write_network_files(comparison_table, out_dir=None, normalize=True):
    if out_dir is None:
        out_dir = ("/home/tyler/Documents/Chr6_docs/Network/"
                   + datetime.today().strftime("%Y_%m_%d"))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    node_path = f"{out_dir}/nodes_"
    edge_path = f"{out_dir}/edges_"
    file_no = 1
    while os.path.isfile(f"{node_path}{file_no}.csv"):
        file_no += 1
    node_path = f"{node_path}{file_no}.csv"
    edge_path = f"{edge_path}{file_no}.csv"

    nodes = build_network_nodes(comparison_table,)
    edges = build_network_edges(comparison_table)
    write_network_nodes(nodes, node_path, normalize)
    write_network_edges(edges, edge_path, normalize)
    return nodes, edges


# %% NetworkX
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


def filter_graph_edges(graph, length_sim_threshold=0, overlap_sim_threshold=0,
                       gene_sim_threshold=0, hi_gene_sim_threshold=0):
    edges = list(graph.edges.items())
    edges = [(edge[0][0], edge[0][1], edge[1]) for edge in edges
             if edge[1]["length_sim"] >= length_sim_threshold
             and edge[1]["overlap_sim"] >= overlap_sim_threshold
             and edge[1]["gene_sim"] >= gene_sim_threshold
             and edge[1]["hi_sim"] >= hi_gene_sim_threshold]
    filtered_graph = nx.Graph()
    filtered_graph.add_nodes_from(graph.nodes)
    filtered_graph.add_edges_from(edges)
    return filtered_graph


def get_subnet_sizes(graph):
    components = nx.connected_components(graph)
    sizes = sorted([len(component) for component in components])
    return sizes

def get_subnet_size_counts_over_n(graph, n):
    sizes = get_subnet_sizes(graph)
    five_count = sum((1 for size in sizes if size >= n))
    return five_count


def get_subnet_sizes_and_locations(graph, patient_db):
    subgroups = nx.connected_components(graph)
    middles = []
    for subgroup in subgroups:
        subgroup_genes = {gene for patient in subgroup
                          for gene in patient_db[patient].all_genes()
                          if gene.seqname == "6"}
        if not subgroup_genes:
            continue
        center_gene_locus = int(
            (min([gene.start for gene in subgroup_genes])
             + max([gene.end for gene in subgroup_genes]))
            / 2)
        middles.append([len(subgroup), center_gene_locus])
    middles.sort(key=lambda x: x[0], reverse=True)
    return middles

# %% Histograms
def plot_gene_sim_histograms_with_slider(graph):
    fig = go.Figure()

    for threshold in np.arange(0, 1.05, .05):

        sizes = get_subnet_sizes(
            filter_graph_edges(graph, gene_sim_threshold=threshold)
            )

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


def plot_hi_gene_sim_histograms_with_slider(graph):
    fig = go.Figure()

    for threshold in np.arange(0, 1.05, .05):

        sizes = get_subnet_sizes(
            filter_graph_edges(graph, hi_gene_sim_threshold=threshold)
            )

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
                  {"title": f"Number of groups of size n with HI gene similarity threshold of {i*.05:.0%}"}],  # layout attribute
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


def plot_overlap_sim_histograms_with_slider(graph):
    fig = go.Figure()

    for threshold in np.arange(0, 1.05, .05):

        sizes = get_subnet_sizes(
            filter_graph_edges(graph, overlap_threshold=threshold)
            )

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
                  {"title": f"Number of groups of size n with overlap similarity threshold of {i*.05:.0%}"}],  # layout attribute
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

def plot_length_sim_histograms_with_slider(graph):
    fig = go.Figure()

    for threshold in np.arange(0, 1.05, .05):

        sizes = get_subnet_sizes(
            filter_graph_edges(graph, length_threshold=threshold)
            )

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
                  {"title": f"Number of groups of size n with length similarity threshold of {i*.05:.0%}"}],  # layout attribute
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
