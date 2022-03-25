#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Phenotype Homogeneity.

@author: T.D. Medina
"""

import pandas as pd


class PhenotypeHomogeneity:
    def __init__(self, group_id, patients, phenotype, group_size, prevalence,):
        self.group_id = group_id
        self.patients = patients
        self.phenotype = phenotype
        self.group_size = group_size
        self.prevalence = prevalence
        self.relative_prevalence = 0
        if self.group_size > 0:
            self.relative_prevalence = self.prevalence/self.group_size

    def __repr__(self):
        string = (f"PhenotypeHomogeneity("
                  f"group_id={self.group_id}, "
                  f"phenotype='{self.phenotype.name}', "
                  f"group_size={self.group_size}, "
                  f"prevalence={self.prevalence} ({self.relative_prevalence:.2%}))")
        return string


class GroupPhenotypeHomogeneity:
    def __init__(self, group_id, patients, phenotype_homogeneities):
        self.group_id = group_id
        self.homogeneities = phenotype_homogeneities
        self.patients = patients
        self.group_size = patients.size
        self.phenotypes = list(self.homogeneities.keys())
        self.num_phenotypes = len(self.phenotypes)
        self.num_present = self.count_present_phenotypes()
        self.proportion_present = self.calculate_overall_proportion_present()

    def __str__(self):
        string = (f"GroupPhenotypeHomogeneity(num_phenotypes={self.num_phenotypes}, "
                  f"phenotypes_present={self.num_present} "
                  f"({self.proportion_present:.2%}))")
        return string

    def count_present_phenotypes(self):
        counts = sum([phenotype.prevalence >= 1
                      for phenotype in self.homogeneities.values()])
        return counts

    def count_phenotypes_above_relative_prevalence(self, threshold):
        counts = sum([phenotype.relative_prevalence >= threshold
                      for phenotype in self.homogeneities.values()])
        return counts

    def calculate_homogeneity(self, threshold):
        if self.num_present == 0:
            return 0
        counts = self.count_phenotypes_above_relative_prevalence(threshold)
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


def make_phenotype_homogeneity_table(group_homogeneities, selected_hpos, homogeneity_threshold):
    pheno_data = {hpo: [] for hpo in selected_hpos}
    for group_homo in group_homogeneities.values():
        for pheno, homo in group_homo.homogeneities.items():
            pheno_data[pheno].append(homo)
    table = {phenotype.name:
                 {f"{homo.group_id} ({homo.group_size})": homo.relative_prevalence
                  for homo in sorted(group_homogeneities.values(), key=lambda x: x.group_size, reverse=True)}
             for phenotype, homos in pheno_data.items()}
    table["Phenotype Homogeneity"] = {
        f"{homo.group_id} ({homo.group_size})": homo.calculate_homogeneity(homogeneity_threshold)
        for homo in sorted(group_homogeneities.values(), key=lambda x: x.group_size)
        }
    table = pd.DataFrame.from_dict(table).transpose()
    return table


def phenotype_homo_test(comparison, ontology, hi_similarity, homogeneity_threshold, phenotypes=None):
    if phenotypes is None:
        with open("/home/tyler/Documents/Chr6_docs/Phenotype_Homogeneity/selected_phenotypes_hpos.txt") as infile:
            selected_hpos = infile.readlines()
        selected_hpos = [ontology[x.strip().split("\t")[1]] for x in selected_hpos]
    else:
        selected_hpos = phenotypes
    homos = comparison.test_all_homogeneities(selected_hpos, hi_similarity=hi_similarity)[0]
    table = make_phenotype_homogeneity_table(homos, selected_hpos, homogeneity_threshold)
    return table, selected_hpos, homos
