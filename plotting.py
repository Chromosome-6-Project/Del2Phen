#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Plotting.

@author: Ty Medina
"""

import matplotlib.pyplot as plt
import numpy as np

from utilities import overlap, is_patient


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
#         if is_gene(p1) or is_gene(p2) or p1.id == p2.id:
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
                if is_patient(intersect.patients[0])
                and is_patient(intersect.patients[1])
                and intersect.patients[0] != intersect.patients[1]
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
