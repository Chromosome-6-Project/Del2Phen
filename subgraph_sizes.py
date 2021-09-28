# -*- coding: utf-8 -*-
"""Chromosome 6 Project - Filter Effects on Subgraph Sizes.

Created on Fri Sep 17 11:54:00 2021

@author: Ty
"""

import networkx as nx


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
