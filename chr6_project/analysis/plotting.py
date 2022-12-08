#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Plotting.

@author: T.D. Medina
"""

from math import log10

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


# %% Group Size Plots
def node_degree_histogram(comparison_table, hi_sim=0.75):
    nx_graph = network.make_nx_graph_from_comparisons(comparison_table)
    nx_graph = network.filter_graph_edges(nx_graph, hi_gene_similarity=hi_sim)
    degrees = sorted(network.get_node_degrees(nx_graph).values())
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(x=degrees,
                     xbins=dict(start=0, end=500, size=1),
                     autobinx=False,
                     hovertemplate="Connections: %{x}<br>Patients: %{y}<extra></extra>")
        )
    fig.update_layout(title=("Histogram of Patient Connection Counts, "
                             f"at {hi_sim:.2%} HI Similarity"),
                      xaxis_title="Number of Connections",
                      yaxis_title=f"Number of Patients",
                      legend_title="HI Similarity",
                      hovermode="x unified")
    return fig


def connectivity_histogram(comparison, length_similarity=0, loci_similarity=0,
                           gene_similarity=0, hi_gene_similarity=0, dom_gene_match=True,
                           hpo_similarity=0, min_size=5):
    nx_graph = network.make_nx_graph_from_comparisons(comparison)
    filtered_graph = network.filter_graph_edges(nx_graph,
                                                length_similarity,
                                                loci_similarity,
                                                gene_similarity,
                                                hi_gene_similarity,
                                                dom_gene_match,
                                                hpo_similarity)
    subsizes = network.get_subnet_sizes(filtered_graph)
    subsize_gt_5 = network.count_subnets_over_size_n(filtered_graph, min_size)
    subsize_gt_5 = f"Clusters with at least {min_size} individuals: {subsize_gt_5:>5}"

    links = sorted(list(network.get_node_degrees(filtered_graph).values()))
    links_gt_5 = network.count_nodes_over_n_degree(filtered_graph, min_size)
    links_gt_5 = f"Individuals with at least {min_size} links:    {links_gt_5:>5}"

    params = dict(xbins=dict(start=0, end=500, size=1), autobinx=False)

    fig = go.Figure()
    fig.add_trace(trace=go.Histogram(
        x=subsizes, name="Subnet Sizes",
        hovertemplate="Subnet Size: %{x}<br>Count: %{y}<extra></extra>", **params)
        )
    fig.add_trace(trace=go.Histogram(
        x=links, name="Node Degrees",
        hovertemplate="Links: %{x}<br>Count: %{y}<extra></extra>", **params)
        )

    fig.update_layout(
        transition_duration=500,
        title_text="Distribution of Subnet Sizes and Node Degrees",
        xaxis_title_text="Size",
        yaxis_title_text="Count",
        legend_orientation="h"
        )
    fig.update_yaxes(type="log", dtick=log10(2))
    fig.add_annotation(text=subsize_gt_5 + "<br>" + links_gt_5, showarrow=False,
                       align="right", xref="paper", yref="paper", x=0.95, y=0.95)
    return fig


def plot_median_degree_vs_hi_similarity(comparison_table):
    nx_graph = network.make_nx_graph_from_comparisons(comparison_table)
    degrees = []
    xs = np.linspace(0, 1, 101)
    for x in xs:
        filtered_graph = network.filter_graph_edges(nx_graph, hi_gene_similarity=x)
        median_degrees = list(network.get_node_degrees(filtered_graph).values())
        median_degrees = np.median(median_degrees)
        degrees.append(median_degrees)
    df = pd.DataFrame({"HI Gene Similarity": xs, "Median Node Degree": degrees})
    fig = px.line(df, x="HI Gene Similarity", y="Median Node Degree",
                  title="Median Node Degree vs. HI Similarity")
    # fig.show()
    # return
    return fig


def plot_min_degree_count_vs_hi_score(comparison, hi_scores, min_degrees,
                                      dom_gene_match=True):
    fig = go.Figure()
    for min_degree in min_degrees:
        graph = network.make_nx_graph_from_comparisons(comparison)
        points = []
        for hi_score in hi_scores:
            filtered = network.filter_graph_edges(graph, hi_gene_similarity=hi_score,
                                                  dom_gene_match=dom_gene_match)
            degree_count = network.count_nodes_over_n_degree(filtered, min_degree)
            points.append(degree_count)
        fig.add_trace(go.Scatter(x=hi_scores, y=points, name=f"n ≥ {min_degree}"))
    if len(min_degrees) == 1:
        n = min_degrees[0]
    else:
        n = "n"
    fig.update_layout(title="Number of patients with minimum connections vs. minHIGSS",
                      xaxis_title="Minimum HI gene similarity score",
                      yaxis_title=f"Number of patients with ≥ {n} connections",
                      legend_title="Connection Threshold (n)",
                      hovermode="x unified")
    return fig


# def plot_min_degree_count_vs_genetic_similarity(comparison, min_degrees,
#                                                 genetic_similarities):
#     fig = go.Figure()
#     for min_degree in min_degrees:
#         graph = network.make_nx_graph_from_comparisons(comparison)
#         points = []
#         for gen_sim in genetic_similarities:
#             filtered = network.filter_graph_edges(graph, gen_sim)
#             degree_count = network.count_nodes_over_n_degree(filtered, min_degree)
#             points.append(degree_count)
#         fig.add_trace(go.Scatter(x=hi_scores, y=points, name=f"n ≥ {min_degree}"))
#     fig.update_layout(title="Number of patients with minimum connections vs. minHIGSS",
#                       xaxis_title="Minimum HI gene similarity score",
#                       yaxis_title="Number of patients with ≥ n connections",
#                       legend_title="Connection Threshold (n)",
#                       hovermode="x unified")
#     return fig


def ph_score_vs_group_size(comparison, phenotypes, hi_scores, dom_gene_match=True,
                           rel_threshold=0.2, abs_threshold=2, min_size=5):
    thresh = (rel_threshold, abs_threshold)
    fig = go.Figure()
    for hi_score in sorted(hi_scores, reverse=True):
        group_phs = comparison.test_all_homogeneities(phenotypes=phenotypes,
                                                      hi_gene_similarity=hi_score,
                                                      dom_gene_match=dom_gene_match)
        data = [(group_ph.group_size, group_ph.calculate_homogeneity(*thresh))
                for group_ph in group_phs.values() if group_ph.group_size >= min_size]
        data.sort(key=lambda x: x[0])
        xs, ys = zip(*data)
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", name=f"HI ≥ {hi_score}"))
    title = f"PH Score vs Group Size, for groups ≥ {min_size}"
    if dom_gene_match:
        title += " with matching DE genes"
    fig.update_layout(title=title,
                      xaxis_title="Group size",
                      yaxis_title="PH score",
                      legend_title="HI Threshold")
    return fig


# %% Patients per PH score
def patients_per_ph(comparison, phenotypes,
                    hi_scores=(0.5, 0.6, 0.7, 0.75, 0.8, 0.9),
                    dom_gene_match=True,
                    rel_threshold=0.2, abs_threshold=2, min_size=5,
                    raw_patient_count=False):
    thresh = (rel_threshold, abs_threshold)
    count_name = "Fraction"
    hovertemplate = "PH ≥ %{x}: %{y:.2%} (%{customdata})"
    if raw_patient_count:
        count_name = "Number"
        hovertemplate = "PH ≥ %{x}: %{y} (%{customdata:.2%})"
    fig = go.Figure()
    for n, hi_score in enumerate(sorted(hi_scores, reverse=True)):
        group_phs = comparison.compare_all_patient_pheno_prevalences(
            phenotypes=phenotypes,
            hi_gene_similarity=hi_score,
            dom_gene_match=dom_gene_match
            )
        data, _ = group_phs.calculate_all_homogeneities(*thresh, min_size)
        data = [sum([x >= i for x in data.values()]) for i in np.linspace(0, 1, 51)]
        customdata = [point/group_phs.size for point in data]
        if not raw_patient_count:
            data, customdata = customdata, data
        fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=data,
                                 name=f"HI ≥ {hi_score}", customdata=customdata,
                                 line_color=GA_COLORS[n],
                                 hovertemplate=hovertemplate))
    fig.update_layout(title=f"{count_name} of Patients with Minimum PH Score,"
                            f" for groups ≥ {min_size}",
                      xaxis_title="Minimum PH score",
                      yaxis_title=f"{count_name} of Patients",
                      legend_title="HI Threshold",
                      hovermode="x unified")
    return fig


def ph_histogram(comparison=None, phenotypes=None, existing_ph_database=None,
                 length_similarity=0, loci_similarity=0, gene_similarity=0,
                 hi_gene_similarity=0, dom_gene_match=True, hpo_similarity=0,
                 rel_threshold=0.2, abs_threshold=2, min_size=5,
                 raw_patient_count=True):
    if existing_ph_database is not None:
        group_phs = existing_ph_database
    else:
        score_thresholds = dict(length_similarity=length_similarity,
                                loci_similarity=loci_similarity,
                                gene_similarity=gene_similarity,
                                hi_gene_similarity=hi_gene_similarity,
                                dom_gene_match=dom_gene_match,
                                hpo_similarity=hpo_similarity)
        group_phs = comparison.compare_all_patient_pheno_prevalences(phenotypes=phenotypes,
                                                                     **score_thresholds)

    upper, lower = group_phs.calculate_all_homogeneities(rel_threshold, abs_threshold, min_size)
    upper = [100*score for score in upper.values()]
    lower = [100*score for score in lower.values()]

    histnorm = "" if raw_patient_count else "percent"
    count_name = "Count" if raw_patient_count else "Percent"

    params = dict(autobinx=False, xbins=dict(start=0, size=5), histnorm=histnorm,
                  hovertemplate="Homogeneity Score: %{x}%<br>Count: %{y}")

    fig = go.Figure()
    fig.add_trace(go.Histogram(x=upper, name=f"Group >= {min_size}", **params))
    fig.add_trace(go.Histogram(x=lower, name=f"Group < {min_size}", **params))

    fig.update_layout(title=f"Histogram of Patient PH Scores",
                      xaxis_title="PH score",
                      yaxis_title=f"{count_name} of Patients",
                      barmode="stack")
    return fig


def patients_per_ph_w_error(comparison, phenotypes,
                            hi_scores=(0.5, 0.6, 0.7, 0.75, 0.8, 0.9),
                            rel_threshold=0.2, abs_threshold=2,
                            min_size=5, random_combinations=50):
    thresh = (rel_threshold, abs_threshold)
    fig = go.Figure()
    x_axis = np.linspace(0, 1, 51)

    for i, hi_score in enumerate(sorted(hi_scores, reverse=True)):
        # TODO: Removed upper/lower homogen output, so this is broken.
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_gene_similarity=hi_score,
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
            # TODO: Removed upper/lower homogen output, so this is broken.
            _, group_phs, _ = comparison.test_all_homogeneities(
                phenotypes=phenos,
                hi_gene_similarity=hi_score,
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


def patients_per_ph_by_area(comparison, phenotypes,
                            hi_scores=(0.5, 0.6, 0.7, 0.75, 0.8, 0.9),
                            rel_threshold=0.2, abs_threshold=2, min_size=5,
                            raw_patient_count=False):
    thresh = (rel_threshold, abs_threshold)
    count_name = "Fraction"
    hovertemplate = "PH ≥ %{x}: %{y:.2%} (%{customdata})"
    if raw_patient_count:
        hovertemplate = "PH ≥ %{x}: %{y} (%{customdata:.2%})"
        count_name = "Number"
    divider = 57038355
    fig = make_subplots(1, 3)

    for n, hi_score in enumerate(sorted(hi_scores, reverse=True)):
        # TODO: Removed upper/lower homogen output, so this is broken.
        _, group_phs, _ = comparison.test_all_homogeneities(
            phenotypes=phenotypes,
            hi_gene_similarity=hi_score,
            group_size_threshold=min_size
            )

        region_data = [[], [], []]
        for group_ph in group_phs.values():
            median_locus = group_ph.patients.get_median_cnv_position("6")
            region = (divider < median_locus) + (divider*2 < median_locus)
            region_data[region].append(group_ph)

        for sub, data in enumerate(region_data, start=1):
            data = [group_ph.calculate_homogeneity(*thresh) for group_ph in data]
            data_len = len(data)
            data = [sum([x >= i for x in data]) for i in np.linspace(0, 1, 51)]
            perc = [point/data_len for point in data]

            if not raw_patient_count:
                data, perc = perc, data
            fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=data,
                                     name=f"HI ≥ {hi_score}", customdata=perc,
                                     hovertemplate=hovertemplate,
                                     legendgroup=f"{hi_score}",
                                     showlegend=(not bool(sub-1)),
                                     line_color=GA_COLORS[n]),
                          row=1, col=sub)

    fig.update_layout(title=f"{count_name} of Patients with Minimum PH Score,"
                            f" for groups ≥ {min_size}",
                      xaxis_title="Minimum PH score",
                      yaxis_title=f"{count_name} of Patients",
                      legend_title="HI Threshold",
                      hovermode="x unified")
    fig.update_traces()
    return fig


# %% Phenotype prediction plots
def plot_precision_stats(prediction_database, patient_database, phenotypes=None,
                         rel_threshold=0.2, abs_threshold=2, use_adjusted_frequency=True,
                         group_size_threshold=5, geneset=None, color_de_patients=False):

    params = (phenotypes, rel_threshold, abs_threshold, use_adjusted_frequency,
              group_size_threshold)

    prediction_stats = prediction_database.calculate_individual_precision(*params)
    positions = [patient_database[patient].get_median_cnv_position("6")
                 for patient in prediction_stats]
    sizes = [prediction_database.predictions[patient].patient_group.size
             for patient in prediction_stats]

    if color_de_patients:
        de_patients = [patient.id for patient in patient_database
                       if patient.all_dominant_genes()
                       and patient.id in prediction_stats]

    prediction_stats = pd.DataFrame.from_dict(prediction_stats, orient="index")
    prediction_stats["Position"] = positions
    prediction_stats["Size"] = sizes
    prediction_stats = prediction_stats.sort_values("Position")

    if geneset:
        fig = plot_genes(geneset, "6", 0, 1)
    else:
        fig = go.Figure()
    hovertemplate = "%{customdata}<br>Pos: %{x}<br>Size: %{marker.size}<br>Value: %{y}"
    for stat in prediction_stats.columns[:-2]:
        fig.add_trace(go.Scatter(x=prediction_stats["Position"],
                                 y=prediction_stats[stat],
                                 name=stat,
                                 mode="markers",
                                 marker=dict(size=prediction_stats["Size"]),
                                 customdata=list(prediction_stats.index),
                                 hovertemplate=hovertemplate,
                                 legendgroup=stat))
        if color_de_patients:
            fig.add_trace(go.Scatter(x=prediction_stats.loc[de_patients]["Position"],
                                     y=prediction_stats.loc[de_patients][stat],
                                     name=f"{stat}_DE",
                                     mode="markers",
                                     marker=dict(size=prediction_stats.loc[de_patients]["Size"],
                                                 color="rgba(0,0,0,0)",
                                                 line=dict(color="red", width=2)),
                                     showlegend=False,
                                     legendgroup=stat))

    fig.update_layout(title="Phenotype Prediction Precision Metrics",
                      xaxis_title="Chromosome 6 Position",
                      yaxis_title="Precision Value")
    return fig


def minimum_precision(comparison, phenotypes=None,
                      rel_threshold=0.2, abs_threshold=2, use_adjusted_frequency=True,
                      group_size_threshold=5, hi_scores=(0.5, 0.6, 0.7, 0.75, 0.8, 0.9),
                      dom_gene_match=True):
    params = (phenotypes, rel_threshold, abs_threshold, use_adjusted_frequency,
              group_size_threshold)

    fig = go.Figure()

    for i, hi_score in enumerate(hi_scores):
        prediction_database = comparison.test_all_phenotype_predictions(
            hi_gene_similarity=hi_score,
            dom_gene_match=dom_gene_match
            )
        stats = prediction_database.calculate_individual_precision(*params)
        stats = pd.DataFrame.from_dict(stats, orient="index")

        size = len(stats)

        ppv = [sum(stats["PPV"] >= j)/size for j in np.linspace(0, 1, 51)]
        sensitivity = [sum(stats["Sensitivity"] >= j)/size for j in np.linspace(0, 1, 51)]

        fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=ppv, mode="lines",
                                 line=dict(color=GA_COLORS[i]),
                                 name=f"minHIGGS ≥ {hi_score:.2}", legendgroup="PPV",
                                 legendgrouptitle_text="PPV"))
        fig.add_trace(go.Scatter(x=np.linspace(0, 1, 51), y=sensitivity, mode="lines",
                                 name=f"minHIGGS ≥ {hi_score:.2}",
                                 line=dict(dash="dot", color=GA_COLORS[i]),
                                 legendgroup="Sensitivity",
                                 legendgrouptitle_text="Sensitivity"))
    fig.update_layout(title=f"Fraction of Patients with Minimum Prediction Precision,"
                            f" for groups ≥ {group_size_threshold}",
                      xaxis_title="Statistic Value",
                      yaxis_title="Fraction of patients")
    fig.update_layout(legend=dict(groupclick="toggleitem"))
    return fig


def plot_genes(geneset, chromosome):
    cats = {0: ["blue", "Regular", "legendonly"],
            1: ["orange", "HI", True],
            2: ["red", "DE", True]}
    legend = [True, True, True]
    fig = go.Figure()
    for i, gene in enumerate(geneset.chromosomes[chromosome].values()):
        cat = max(gene.dominant*2, gene.is_haploinsufficient())
        color = cats[cat][0]
        fig.add_trace(go.Scatter(
            x=[gene.start, gene.end, gene.end, gene.start, gene.start],
            y=[0, 0, 1, 1, 0],
            fill="toself",
            mode="lines",
            line_color=color,
            legendgroup=cats[cat][1],
            name=f"{cats[cat][1]}",
            showlegend=legend[cat],
            visible=cats[cat][2]
            ))
        legend[cat] = False
    return fig


def plot_genes(geneset, chromosome, y_min=0.0, y_max=1.0):
    cats = {0: ["blue", "Regular", "legendonly"],
            1: ["orange", "HI", True],
            2: ["red", "DE", True]}

    dom_genes = [gene for gene in geneset.chromosomes[chromosome].values()
                 if gene.dominant]
    hi_genes = [gene for gene in geneset.chromosomes[chromosome].values()
                if gene.is_haploinsufficient()]
    other_genes = set(geneset.chromosomes[chromosome].values()) - set(dom_genes) - set(hi_genes)

    fig = go.Figure()
    for i, genes in enumerate([other_genes, hi_genes, dom_genes]):
        xs = [pos for gene in genes
              for pos in [gene.start, gene.end, gene.end, gene.start, gene.start, None]]
        ys = [pos for _ in genes
              for pos in [y_min, y_min, y_max, y_max, y_min, None]]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, fill="toself", mode="lines",
            line_color=cats[i][0],
            visible=cats[i][2],
            name=cats[i][1],
            opacity=0.5
            ))
    return fig
