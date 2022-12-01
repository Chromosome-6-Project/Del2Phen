#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Phenotype Homogeneity.

@author: T.D. Medina
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

pio.renderers.default = "browser"


class PhenotypePrevalence:
    def __init__(self, group_id, patients, phenotype, group_size, prevalence):
        self.group_id = group_id
        self.patients = patients
        self.phenotype = phenotype
        self.group_size = group_size
        self.prevalence = prevalence
        self.relative_prevalence = 0
        if self.group_size > 0:
            self.relative_prevalence = self.prevalence/self.group_size

    def __repr__(self):
        string = (f"PhenotypePrevalence("
                  f"group_id={self.group_id}, "
                  f"phenotype='{self.phenotype.name}', "
                  f"group_size={self.group_size}, "
                  f"prevalence={self.prevalence} ({self.relative_prevalence:.2%}))")
        return string


class PatientGroupPrevalences:
    def __init__(self, group_id, patients, phenotype_prevalences):
        self.group_id = group_id
        self.prevalences = phenotype_prevalences
        self.patients = patients
        self.group_size = patients.size
        self.phenotypes = list(self.prevalences.keys())
        self.num_phenotypes = len(self.phenotypes)
        self.num_present = self.count_present_phenotypes()
        self.proportion_present = self.calculate_overall_proportion_present()

    def __str__(self):
        string = (f"PatientGroupPrevalences(num_phenotypes={self.num_phenotypes}, "
                  f"phenotypes_present={self.num_present} "
                  f"({self.proportion_present:.2%}))")
        return string

    def count_present_phenotypes(self):
        counts = sum([phenotype.prevalence >= 1
                      for phenotype in self.prevalences.values()])
        return counts

    def count_phenotypes_above_prevalence_thresholds(self, rel_threshold, abs_threshold):
        counts = sum([phenotype.relative_prevalence >= rel_threshold
                      and phenotype.prevalence >= abs_threshold
                      for phenotype in self.prevalences.values()])
        return counts

    def calculate_homogeneity(self, rel_threshold, abs_threshold):
        if self.num_present == 0:
            return 0
        counts = self.count_phenotypes_above_prevalence_thresholds(rel_threshold, abs_threshold)
        homogeneity = counts / self.num_present
        return homogeneity

    def calculate_overall_proportion_present(self):
        prevalence = self.num_present / self.num_phenotypes
        return prevalence

    # def set_group_size(self):
    #     group_size = {x.group_size for x in self.homogeneities.values()}
    #     if len(group_size) != 1:
    #         raise ValueError("PhenotypeHomogeneities have different group sizes.")
    #     group_size = list(group_size)[0]
    #     return group_size
    #
    # def set_group_id(self):
    #     group_id = {x.group_id for x in self.homogeneities.values()}
    #     if len(group_id) != 1:
    #         raise ValueError("PhenotypeHomogeneities have different group IDs.")
    #     group_id = list(group_id)[0]
    #     return group_id


class HomogeneityDatabase:
    def __init__(self, patient_group_prevalences):
        self.patient_group_prevalences = patient_group_prevalences
        self.phenotypes = {pheno for group_prev in self.patient_group_prevalences.values()
                           for pheno in group_prev.phenotypes}
        self.size = len(self.patient_group_prevalences)

    def sorted_by_chromosome(self, chromosome):
        sorted_prevs = sorted(self.patient_group_prevalences.values(),
                              key=lambda x: x.patients.get_median_cnv_position(chromosome))
        return sorted_prevs

    def calculate_all_homogeneities(self, rel_threshold=0.2, abs_threshold=2,
                                    split_by_size=None):
        if not split_by_size:
            scores = {pid: group_ph.calculate_homogeneity(rel_threshold, abs_threshold)
                      for pid, group_ph in self.patient_group_prevalences.items()}
            return scores
        upper = {pid: group_ph.calculate_homogeneity(rel_threshold, abs_threshold)
                 for pid, group_ph in self.patient_group_prevalences.items()
                 if group_ph.group_size >= split_by_size}
        lower = {pid: group_ph.calculate_homogeneity(rel_threshold, abs_threshold)
                 for pid, group_ph in self.patient_group_prevalences.items()
                 if group_ph.group_size < split_by_size}
        return upper, lower

    def make_phenotype_homogeneity_table(self, rel_threshold=0.2, abs_threshold=2,
                                         min_size=0):
        pheno_data = {pheno: [] for pheno in self.phenotypes}

        group_prevs = self.sorted_by_chromosome("6")
        group_prevs = [group_prev for group_prev in group_prevs
                       if group_prev.group_size >= min_size]
        for group_prev in group_prevs:
            for pheno, pheno_prev in group_prev.prevalences.items():
                pheno_data[pheno].append(pheno_prev)

        table = {phenotype.name:
                 {f"{pheno_prev.group_id} ({pheno_prev.group_size})": pheno_prev.relative_prevalence
                  for pheno_prev in pheno_prevs}
                 for phenotype, pheno_prevs in pheno_data.items()}
        table["Phenotype Homogeneity"] = {
            f"{group_prev.group_id} ({group_prev.group_size})":
            group_prev.calculate_homogeneity(rel_threshold, abs_threshold)
            for group_prev in group_prevs
            }
        table = pd.DataFrame.from_dict(table).transpose()
        return table


# def make_phenotype_homogeneity_table(group_homogeneities, selected_hpos, rel_threshold, abs_threshold):
#     pheno_data = {hpo: [] for hpo in selected_hpos}
#     group_homogeneities = sorted(group_homogeneities.values(),
#                                  key=lambda x: x.patients.get_median_cnv_position("6"))
#     for group_homo in group_homogeneities:
#         for pheno, homo in group_homo.homogeneities.items():
#             pheno_data[pheno].append(homo)
#     table = {phenotype.name:
#                  {f"{homo.group_id} ({homo.group_size})": homo.relative_prevalence
#                   for homo in homos}
#              for phenotype, homos in pheno_data.items()}
#     table["Phenotype Homogeneity"] = {
#         f"{homo.group_id} ({homo.group_size})": homo.calculate_homogeneity(rel_threshold, abs_threshold)
#         for homo in group_homogeneities
#         }
#     table = pd.DataFrame.from_dict(table).transpose()
#     return table


# def phenotype_homo_test(comparison, selected_hpos,
#                         length_similarity, loci_similarity, gene_similarity,
#                         hi_gene_similarity, hpo_similarity, group_size_threshold,
#                         rel_threshold, abs_threshold):
#     # if phenotypes is None:
#     #     with open("/home/tyler/Documents/Chr6_docs/Phenotype_Homogeneity/selected_phenotypes_hpos.txt") as infile:
#     #         selected_hpos = infile.readlines()
#     #     selected_hpos = [ontology[x.strip().split("\t")[1]] for x in selected_hpos]
#     # else:
#     #     selected_hpos = phenotypes
#     homos = comparison.calculate_all_patient_homogeneities(
#         phenotypes=selected_hpos,
#         length_similarity=length_similarity,
#         loci_similarity=loci_similarity, gene_similarity=gene_similarity,
#         hi_gene_similarity=hi_gene_similarity, hpo_similarity=hpo_similarity,
#         group_size_threshold=group_size_threshold
#         )
#     table = make_phenotype_homogeneity_table(homos, selected_hpos, rel_threshold, abs_threshold)
#     return table, homos


def plot_phenotype_homogeneities(comparison, selected_hpos, hi_gene_similarity, abs_threshold):
    homos = comparison.calculate_all_patient_homogeneities(selected_hpos, hi_gene_similarity=hi_gene_similarity)
    homos = sorted(homos.values(),
                   key=lambda x: x.patients.get_median_cnv_position("6"))
    table = {}
    for i in reversed(np.linspace(0.01, 1, 100)):
        table[f"{i:.2}"] = {
            f"{homo.group_id} ({homo.group_size})": homo.calculate_homogeneity(i, abs_threshold)
            for homo in homos
            }
    table = pd.DataFrame(table).transpose()
    fig = px.imshow(table, aspect="auto")
    fig.show()
    return table


def plot_phenotype_homogeneities_vs_hi_sim(comparison, selected_hpos, hi_similarities, abs_threshold):
    size = len(comparison)
    homo_scores = []
    cnv_positions = []
    ids = []
    sizes = []
    for hi_similarity in hi_similarities:
        homos = comparison.test_all_homogeneities(selected_hpos, hi_gene_similarity=hi_similarity)
        homo_scores += [homo.calculate_homogeneity(0.2, abs_threshold)
                        for homo in homos.values()]
        cnv_positions += [homo.patients.get_median_cnv_position("6")
                          for homo in homos.values()]
        ids += [homo.group_id for homo in homos.values()]
        sizes += [homo.group_size for homo in homos.values()]
    series = [n for i in hi_similarities for n in [i]*size]
    x_pos = [n for _ in hi_similarities for n in range(1, size+1)]
    table = pd.DataFrame(zip(series, cnv_positions, homo_scores, ids, sizes),
                         columns=["HI Sim Score", "CNV Pos", "PH Score", "ID", "Size"])
    table.sort_values(by=["HI Sim Score", "PH Score", "CNV Pos"], inplace=True,
                      ascending=False)
    fig = px.line(table, x=x_pos, y="PH Score", color="HI Sim Score",
                  markers=True, hover_name="ID", hover_data=["Size"])
    fig.show()
    return fig


def plot_phenotype_homogeneities_vs_hi_sim2(comparison, selected_hpos, hi_similarities, abs_threshold, min_size):
    mid_hi_sim = hi_similarities[len(hi_similarities) // 2]
    size = len(comparison)
    homo_scores = []
    cnv_positions = []
    ids = []
    sizes = []

    homos = comparison.test_all_homogeneities(selected_hpos, hi_gene_similarity=mid_hi_sim)
    homos = [(homo.group_id, homo.group_size, homo.calculate_homogeneity(0.2, abs_threshold))
             for homo in homos.values()]
    homos.sort(key=lambda x: x[1], reverse=True)
    homos.sort(key=lambda x: x[2], reverse=True)
    sort_order = [y[0] for y in homos]

    for hi_similarity in hi_similarities:
        homos = comparison.test_all_homogeneities(selected_hpos, hi_gene_similarity=hi_similarity)
        homo_scores += [homo.calculate_homogeneity(0.2, abs_threshold)
                        for homo in homos.values()]
        cnv_positions += [homo.patients.get_median_cnv_position("6")
                          for homo in homos.values()]
        ids += [homo.group_id for homo in homos.values()]
        sizes += [homo.group_size for homo in homos.values()]
    series = [n for i in hi_similarities for n in [i]*size]
    x_pos = [n for _ in hi_similarities for n in range(1, size+1)]
    data = list(zip(series, cnv_positions, homo_scores, ids, sizes))
    data.sort(key=lambda x: x[-1])
    data.sort(key=lambda x: x[1])
    data.sort(key=lambda x: sort_order.index(x[3]))
    data.sort(key=lambda x: x[0], reverse=True)
    # data = pd.DataFrame(data)
    # return data
    # data["x_pos"] = x_pos
    table = pd.DataFrame(data, columns=["HI Sim Score", "CNV Pos", "PH Score", "ID", "Size"])
    table["x_pos"] = x_pos
    table = table[table["Size"] >= min_size]
    fig = px.line(table, x="x_pos", y="PH Score", color="HI Sim Score",
                  markers=True, hover_name="ID", hover_data=["Size"])
    fig.update_layout(hovermode="x unified")
    fig.show()
    return fig
