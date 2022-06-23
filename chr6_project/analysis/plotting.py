#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Plotting.

@author: T.D. Medina
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

from chr6_project.analysis import network
from chr6_project.analysis.utilities import overlap


pio.renderers.default = "browser"

GA_COLORS = ["#F0A3FF", "#FF0010", "#2BCE48", "#FFCC99", "#FFFF80",
             "#100AFF", "#5EF1F2", "#990000", "#C20088", "#003380",
             "#426600", "#19A405", "#FFA8BB", "#0075DC", "#808080",
             "#FFE100", "#8F7C00", "#94FFB5", "#4C005C", "#00998F",
             "#FF5000", "#993F00", "#005C31", "#9DCC00", "#191919",
             "#E0FF66"]


def hex_to_rgb(hex):
    rgb = tuple(int(hex[i:i+2], 16) for i in range(1, 6, 2))
    return rgb


def hex_to_rgba(hex, alpha):
    rgba = tuple([int(hex[i:i+2], 16) for i in range(1, 6, 2)] + [alpha])
    return rgba


# XXX: Deprecated.
def make_array(table):
    """Make array from comparisons."""
    array = []
    for pid in table:
        values = []
        for p2 in table:
            if pid == p2:
                values.append(1)
            elif p2 in table[pid]:
                values.append(table[pid][p2])
            else:
                values.append(table[p2][pid])
        array.append(values)
    array = np.array(array)
    return array


# XXX: Deprecated.
def gene_comparison_heatmap(table):
    """Plot heatmap of gene comparisons."""
    data_array = make_array(table)

    fig, ax = plt.subplots()
    ax.imshow(data_array)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(table)))
    ax.set_yticks(np.arange(len(table)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(list(table.keys()))
    ax.set_yticklabels(list(table.keys()))

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    # for i in range(len(vegetables)):
    #     for j in range(len(farmers)):
    #         text = ax.text(j, i, harvest[i, j],
    #                        ha="center", va="center", color="w")

    ax.set_title("Test Heat Map")
    fig.tight_layout()
    plt.show()


# def similarity_vs_hpo_scatter(patient_comparison):
#     points = []
#     for intersect in patient_comparison:
#         p1 = patient_comparison.patient_db[intersect.patients[0]]
#         p2 = patient_comparison.patient_db[intersect.patients[1]]
#         if p1.id == p2.id:
#             continue
#         points.append((intersect.gene_similarity, intersect.hpo_count))
#     points = list(zip(*points))
#     plt.scatter(*points)
#     plt.ylabel("Shared HPO terms", fontsize=24)
#     plt.xlabel("Shared genes (percent similarity)", fontsize=24)
#     plt.yticks(fontsize=20)
#     plt.xticks(fontsize=20)


def plot_individual_factors(comparison_table, percentage=True):
    """Plot scatterplots for similarity vs. shared HPO terms."""
    plotters = [intersect for intersect in comparison_table
                if intersect.patients[0] != intersect.patients[1]
                and intersect.patients[0].hpo and intersect.patients[1].hpo
                and intersect.patients[0].cnvs and intersect.patients[1].cnvs]

    # if percentage:
    #     hpos = [x.hpo_similarity for x in plotters]
    # else:
    #     hpos = [x.hpo_count for x in plotters]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    ax1.scatter(*list(zip(*[(x.length_similarity, x.hpo_similarity) for x in plotters if x.length_similarity > 0])))
    ax2.scatter(*list(zip(*[(x.loci_similarity, x.hpo_similarity) for x in plotters if x.loci_similarity > 0])))
    ax3.scatter(*list(zip(*[(x.gene_similarity, x.hpo_similarity) for x in plotters if x.gene_similarity > 0])))

    # ax4.scatter(*list(zip(*[(x.length_similarity, x.hpo_count) for x in plotters if x.length_similarity > 0])))
    # ax5.scatter(*list(zip(*[(x.loci_similarity, x.hpo_count) for x in plotters if x.loci_similarity > 0])))
    # ax6.scatter(*list(zip(*[(x.gene_similarity, x.hpo_count) for x in plotters if x.gene_similarity > 0])))

    # ax1.scatter([x.length_similarity for x in plotters], hpos)
    # ax2.scatter([x.loci_similarity for x in plotters if x.loci_similarity > 0], hpos)
    # ax3.scatter([x.gene_similarity for x in plotters if x.gene_similarity > 0], hpos)

    fig.suptitle("CNV Similarity vs. Shared HPO Terms", fontsize=30)
    y_label = "HPO term Jaccard similarity"

    ax1.set_ylabel(y_label, fontsize=24)

    ax1.set_xlabel("Length Similarity", fontsize=24)
    ax2.set_xlabel("Loci Similarity", fontsize=24)
    ax3.set_xlabel("Gene Similarity", fontsize=24)


def make_cnv_histogram_info(patient_db, chromosome, genome_dict):
    """Gather data on CNVs per megabase for a histogram."""
    cnvs = patient_db.cnvs[chromosome]

    bin_starts = range(1, genome_dict["6"].length, 1000000)
    bin_counts = {range(bin_start, bin_start + 1000000): 0
                  for bin_start in bin_starts[:-1]}
    bin_counts[range(bin_starts[-1], genome_dict["6"].length + 1)] = 0

    for cnv in cnvs:
        for bin_range in bin_counts:
            if overlap(cnv.range, bin_range):
                bin_counts[bin_range] += 1
    return bin_counts


def plot_cnv_histogram(hist_info):
    """Plot histogram of CNV coverage per megabase."""
    bins = [(x.start - 1)/1000000 for x in hist_info]
    heights = list(hist_info.values())

    # plt.bar(y_pos, performance, align='center', alpha=0.5)
    plt.bar(bins, heights, align="center", alpha=0.5, color="seagreen")
    plt.xticks(fontsize=30, fontname="Tahoma")
    plt.yticks(fontsize=30, fontname="Tahoma")
    plt.xlabel("Megabase", fontsize=36, fontname="Tahoma")
    plt.ylabel("Count", fontsize=36, fontname="Tahoma")
    plt.title("CNV coverage per chromosome 6 megabase",
              fontsize=36, fontname="Tahoma")

    plt.show()


def test_plot(hist_info, thing):
    """Plot histogram of CNV coverage per megabase."""
    bins = [(x.start - 1)/1000000 for x in hist_info]
    heights = list(hist_info.values())

    fig = go.Figure()
    fig.add_trace(go.Bar(x=bins, y=heights))
    fig.add_trace(go.Scatter(x=[x[1]/1e6 for x in thing],
                             y=[y[0] for y in thing],
                             mode="markers"))
    fig.show()


def plot_phenotype_homogeneity_heatmap(phenohomo_data):
    fig = px.imshow(phenohomo_data, aspect="auto")
    fig.show()
    return


def plot_median_degree_vs_hi_similarity(comparison_table):
    nx_graph = network.make_nx_graph_from_comparisons(comparison_table)
    degrees = []
    xs = np.linspace(0, 1, 101)
    for x in xs:
        filtered_graph = network.filter_graph_edges(nx_graph, hi_gene_sim_threshold=x)
        median_degrees = list(network.get_node_degrees(filtered_graph).values())
        median_degrees = np.median(median_degrees)
        degrees.append(median_degrees)
    df = pd.DataFrame({"HI Gene Similarity": xs, "Median Node Degree": degrees})
    fig = px.line(df, x="HI Gene Similarity", y="Median Node Degree",
                  title="Median Node Degree vs. HI Similarity")
    # fig.show()
    # return
    return fig

# def plot_phenotype_homogeneity_vs_hi_similarity(comparison_table, hi_sim_values):


def plot_min_degree_count_vs_hi_score(comparison, hi_scores, min_degrees):
    fig = go.Figure()
    for min_degree in min_degrees:
        graph = network.make_nx_graph_from_comparisons(comparison)
        points = []
        for hi_score in hi_scores:
            filtered = network.filter_graph_edges(graph, hi_gene_sim_threshold=hi_score)
            degree_count = network.count_nodes_over_n_degree(filtered, min_degree)
            points.append(degree_count)
        fig.add_trace(go.Scatter(x=hi_scores, y=points, name=f"n ≥ {min_degree}"))
    fig.update_layout(title="Number of patients with minimum connections vs. minHIGSS",
                      xaxis_title="Minimum HI gene similarity score",
                      yaxis_title="Number of patients with ≥ n connections",
                      legend_title="Connection Threshold (n)",
                      hovermode="x unified")
    return fig


def ph_score_vs_group_size(comparison, phenotypes, hi_scores,
                           rel_threshold=0.2, abs_threshold=2, min_size=5):
    thresh = (rel_threshold, abs_threshold)
    fig = go.Figure()
    for hi_score in sorted(hi_scores, reverse=True):
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_similarity=hi_score,
            group_size_threshold=min_size
            )
        data = [(group_ph.group_size, group_ph.calculate_homogeneity(*thresh))
                for group_ph in group_phs.values()]
        data.sort(key=lambda x: x[0])
        xs, ys = zip(*data)
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", name=f"HI ≥ {hi_score}"))
    fig.update_layout(title=f"PH Score vs Group Size, for groups ≥ {min_size}",
                      xaxis_title="Group size",
                      yaxis_title="PH score",
                      legend_title="HI Threshold")
    return fig


def ph_score_histograms_by_hi_score(comparison, phenotypes, hi_scores,
                                    rel_threshold=0.2, abs_threshold=2, min_size=5):
    thresh = (rel_threshold, abs_threshold)
    fig = go.Figure()
    for hi_score in sorted(hi_scores, reverse=True):
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_similarity=hi_score,
            group_size_threshold=min_size
            )
        data = [group_ph.calculate_homogeneity(*thresh)
                for group_ph in group_phs.values()]
        data = [sum([x >= i for x in data]) for i in np.linspace(0, 1, 51)]
        data_len = len(group_phs)
        perc = [point/data_len for point in data]
        fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=data,
                                 name=f"HI ≥ {hi_score}", customdata=perc,
                                 hovertemplate="PH ≥ %{x}: %{y} (%{customdata:.2%})"))
    fig.update_layout(title=f"Number of Patients with Minimum PH Score,"
                            f" for groups ≥ {min_size}",
                      xaxis_title="Minimum PH score",
                      yaxis_title="Number of Patients",
                      legend_title="HI Threshold",
                      hovermode="x unified")
    return fig


def ph_score_histograms_by_hi_score_w_error(comparison, phenotypes, hi_scores,
                                            rel_threshold=0.2, abs_threshold=2,
                                            min_size=5, random_combinations=50):
    thresh = (rel_threshold, abs_threshold)
    fig = go.Figure()
    x_axis = np.linspace(0, 1, 51)

    for i, hi_score in enumerate(sorted(hi_scores, reverse=True)):
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_similarity=hi_score,
            group_size_threshold=min_size
            )
        data = [group_ph.calculate_homogeneity(*thresh)
                for group_ph in group_phs.values()]
        data = np.array([sum([x >= j for x in data]) for j in x_axis])
        data_len = len(group_phs)
        perc = data / data_len

        error_data = np.ndarray([random_combinations+1, 51])
        error_data[0] = data
        for perm in range(1, random_combinations+1):
            phenos= np.random.choice(phenotypes, 20, replace=False)
            _, group_phs, _ = comparison.test_all_homogeneities(
                phenotypes=phenos,
                hi_similarity=hi_score,
                group_size_threshold=min_size
                )
            error_perm_data = [group_ph.calculate_homogeneity(*thresh)
                               for group_ph in group_phs.values()]
            error_data[perm] = [sum([x >= j for x in error_perm_data])
                                for j in np.linspace(0, 1, 51)]

        error_min = data - error_data.min(axis=0)
        error_max = error_data.max(axis=0) - data
        error_vals = np.concatenate([error_data.max(axis=0),
                                     error_data.min(axis=0)[::-1]])
        custom_data = pd.DataFrame(dict(perc=perc,
                                        error_min=error_min,
                                        error_max=error_max))

        fig.add_trace(go.Scatter(x=x_axis, y=data,
                                 name=f"HI ≥ {hi_score}",
                                 customdata=custom_data,
                                 hovertemplate="PH ≥ %{x}: %{y} " 
                                               "(%{customdata[0]:.2%})"
                                               " +%{customdata[2]}/-%{customdata[1]}",
                                 legendgroup=f"{hi_score}",
                                 line_color=GA_COLORS[i],
                                 ))
        fig.add_trace(go.Scatter(x=np.concatenate([x_axis, x_axis[::-1]]),
                                 y=error_vals,
                                 fill='toself',
                                 line_color="rgba(255,255,255,0)",
                                 fillcolor=f"rgba{hex_to_rgba(GA_COLORS[i], 0.3)}",
                                 hoverinfo="skip",
                                 legendgroup=f"{hi_score}",
                                 showlegend=False
                                 ))
    fig.update_layout(title=f"Number of Patients with Minimum PH Score,"
                            f" for groups ≥ {min_size}",
                      xaxis_title="Minimum PH score",
                      yaxis_title="Number of Patients",
                      legend_title="HI Threshold",
                      hovermode="x unified")
    return fig


def ph_score_histograms_by_hi_score_by_area(comparison, phenotypes, hi_scores,
                                            rel_threshold=0.2, abs_threshold=2, min_size=5):
    thresh = (rel_threshold, abs_threshold)
    divider = 57038355
    fig = make_subplots(1, 3)

    for hi_score in sorted(hi_scores, reverse=True):
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_similarity=hi_score,
            group_size_threshold=min_size
            )

        region_data = [[], [], []]
        for group_ph in group_phs.values():
            median_locus = group_ph.patients.get_median_cnv_position("6")
            region = (divider < median_locus) + (divider*2 < median_locus)
            region_data[region].append(group_ph)

        for sub, data in enumerate(region_data, start=1):
            data = [group_ph.calculate_homogeneity(*thresh) for group_ph in data]
            data = [sum([x >= i for x in data]) for i in np.linspace(0, 1, 51)]
            data_len = len(data)
            perc = [point/data_len for point in data]
            fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=data,
                                     name=f"HI ≥ {hi_score}", customdata=perc,
                                     hovertemplate="PH ≥ %{x}: %{y} (%{customdata:.2%})",
                                     legendgroup=f"{hi_score}"),
                          row=1, col=sub)

    fig.update_layout(title=f"Number of Patients with Minimum PH Score,"
                            f" for groups ≥ {min_size}",
                      xaxis_title="Minimum PH score",
                      yaxis_title="Number of Patients",
                      legend_title="HI Threshold",
                      hovermode="x unified")
    return fig
