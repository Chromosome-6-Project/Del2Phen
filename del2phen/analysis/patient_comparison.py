#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Patient Comparisons.

@author: T.D. Medina
"""

from collections import Counter
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from chr6_project.analysis.data_objects import Patient, PatientDatabase, CNV
from chr6_project.analysis.gene_set import GeneSet
from chr6_project.analysis.phenotype_homogeneity import (
    PhenotypePrevalence,
    PatientGroupPrevalences,
    HomogeneityDatabase
    )
from chr6_project.analysis.phenotype_prediction import (
    TraitPrediction,
    PatientPredictions,
    PredictionDatabase
    )
from chr6_project.analysis.utilities import (
    jaccard,
    length_of_range_intersects,
    merge_range_list
    )


# TODO: Add a method to add a patient.
class ComparisonTable:
    """Data object holding all patient vs. patient comparisons."""

    def __init__(self, patient_db=None,
                 chromosomes: Optional[Union[str, List[str]]] = None,
                 copy_numbers: Optional[Union[int, List[int]]] = None,
                 pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=0.86,
                 comparison_table=None, mode="confirm"):
        if comparison_table is not None:
            self.read_from_existing(comparison_table)
            return
        self.patient_db = patient_db
        self.index, self.array = self.compare_patients(chromosomes, copy_numbers,
                                                       pLI_threshold, HI_threshold,
                                                       phaplo_threshold, mode)
        self.size = len(self.index)

        self.__iteri__ = 0
        self.__iterj__ = 0

    def __len__(self):
        return self.size

    def __iter__(self):
        """Initialize iterable."""
        self.__iteri__ = 0
        self.__iterj__ = 0
        return self

    def __next__(self):
        """Iterate 2-axis iterable."""
        if self.__iteri__ == self.size:
            raise StopIteration
        result = self.array[self.__iteri__][self.__iterj__]
        self.__iterj__ += 1
        if self.__iterj__ == self.size:
            self.__iteri__ += 1
            self.__iterj__ = self.__iteri__
        return result

    def lookup(self, pid1, pid2=None):
        """Look up patient or patient-patient intersect in comparison array."""
        if not pid2:
            return self.array[self.index[pid1]][self.index[pid1]]
        if pid2.lower() == "all":
            return self.array[self.index[pid1]]
        if pid1 not in self.index:
            raise KeyError("ID not found.")
        return self.array[self.index[pid1]][self.index[pid2]]

    def read_from_existing(self, comparison_table):
        self.patient_db = comparison_table.patient_db
        self.array = comparison_table.array
        self.size = comparison_table.size

    # %% === Comparison methods ==========

    @staticmethod
    def _make_comparison_array(comparisons):
        """Convert comparison dictionary to numpy array."""
        index = {j: i for i, j in enumerate(comparisons)}
        array = []
        for patient1 in index:
            values = []
            for patient2 in index:
                if patient2 in comparisons[patient1]:
                    values.append(comparisons[patient1][patient2])
                else:
                    values.append(comparisons[patient2][patient1])
            array.append(values)
        array = np.array(array)
        return index, array

    def compare_patients(self,
                         chromosomes: Optional[Union[str, List[str]]] = None,
                         copy_numbers: Optional[Union[int, List[int]]] = None,
                         pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=0.86,
                         mode="confirm"):
        """Compare all patients to each other."""
        ids = list(self.patient_db.patients.keys())
        comparisons = {}
        while ids:
            id_i = ids.pop()
            patient_i = self.patient_db[id_i]
            patient_comparison = dict()
            patient_comparison[id_i] = self.compare_pair(patient_i, patient_i,
                                                         chromosomes=chromosomes,
                                                         copy_numbers=copy_numbers,
                                                         HI_threshold=HI_threshold,
                                                         phaplo_threshold=phaplo_threshold,
                                                         pLI_threshold=pLI_threshold,
                                                         mode=mode)

            for id_j in ids:
                patient_j = self.patient_db[id_j]
                patient_comparison[id_j] = self.compare_pair(patient_i, patient_j,
                                                             chromosomes=chromosomes,
                                                             copy_numbers=copy_numbers,
                                                             HI_threshold=HI_threshold,
                                                             phaplo_threshold=phaplo_threshold,
                                                             pLI_threshold=pLI_threshold,
                                                             mode=mode)
            comparisons[id_i] = patient_comparison
        index, comparisons = self._make_comparison_array(comparisons)
        return index, comparisons

    def compare_patient_vs_others(self, patient: Patient, save_results=False,
                                  chromosomes: Optional[Union[str, List[str]]] = None,
                                  copy_numbers: Optional[Union[int, List[int]]] = None,
                                  pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=0.86,
                                  mode="confirm", **kwargs):
        """Compare one Patient object against all patients in the PatientDatabase."""
        ids = self.index.keys()
        comparisons = []
        for pid in ids:
            patient_2 = self.patient_db[pid]
            intersect = self.compare_pair(patient, patient_2,
                                          chromosomes=chromosomes,
                                          copy_numbers=copy_numbers,
                                          HI_threshold=HI_threshold,
                                          phaplo_threshold=phaplo_threshold,
                                          pLI_threshold=pLI_threshold,
                                          mode=mode)
            comparisons.append(intersect)
        if save_results:
            self_compare = self.compare_pair(patient, patient,
                                             chromosomes=chromosomes, copy_numbers=copy_numbers,
                                             HI_threshold=HI_threshold,
                                             phaplo_threshold=phaplo_threshold,
                                             pLI_threshold=pLI_threshold,
                                             mode=mode)
            self.patient_db.add_patient(patient)
            self.index[patient.id] = max(self.index.values())+1
            self.array = np.append(self.array,
                                   np.array([comparisons]).T,
                                   axis=1)
            comparisons.append(self_compare)
            self.array = np.append(self.array,
                                   np.array([comparisons]),
                                   axis=0)
        return comparisons

    @classmethod
    def compare_pair(cls, patient_1: Patient, patient_2: Patient,
                     chromosomes: Optional[Union[str, List[str]]] = None,
                     copy_numbers: Optional[Union[int, List[int]]] = None,
                     HI_threshold=10,
                     phaplo_threshold=0.86,
                     pLI_threshold=0.9,
                     mode="confirm"):
        """Compare all metrics between two patients. """
        copy_number = {cnv.copy_number for cnv in patient_1.cnvs + patient_2.cnvs}
        if len(copy_number) == 0:
            copy_number = "N/A"
        elif len(copy_number) > 1:
            copy_number = "Mixed"
        else:
            copy_number = list(copy_number)[0]
        length_compare = cls.compare_pair_length(patient_1, patient_2, chromosomes, copy_numbers)
        loci_compare = cls.compare_pair_loci(patient_1, patient_2, chromosomes, copy_numbers)
        gene_compare = cls.compare_pair_genes(patient_1, patient_2, chromosomes, copy_numbers)
        hi_compare = cls.compare_pair_HI_genes(patient_1, patient_2,
                                               chromosomes=chromosomes,
                                               copy_numbers=copy_numbers,
                                               pLI_threshold=pLI_threshold,
                                               HI_threshold=HI_threshold,
                                               phaplo_threshold=phaplo_threshold,
                                               mode=mode)
        dom_compare = cls.compare_pair_dominant_genes(patient_1, patient_2,
                                                      chromosomes, copy_numbers)
        dom_match = cls.compare_pair_dominant_match(patient_1, patient_2,
                                                    chromosomes, copy_numbers)
        hpo_compare = cls.compare_pair_hpos(patient_1, patient_2)
        comparison = PatientIntersect(
            patient_1, patient_2, copy_number,
            length_compare, loci_compare,
            gene_compare, hi_compare,
            dom_compare, dom_match,
            hpo_compare
            )
        return comparison

    @staticmethod
    def compare_pair_length(patient_1: Patient, patient_2: Patient,
                            chromosomes: Optional[Union[str, List[str]]] = None,
                            copy_numbers: Optional[Union[int, List[int]]] = None):
        """Compare CNV length between two patients."""
        size_1 = sum([cnv.length for cnv in patient_1.filter_cnvs(chromosomes, copy_numbers)])
        size_2 = sum([cnv.length for cnv in patient_2.filter_cnvs(chromosomes, copy_numbers)])
        if size_1 == 0 or size_2 == 0:
            return 0
        return jaccard(size_1, size_2)

    @staticmethod
    def compare_pair_loci(patient_1, patient_2,
                          chromosomes: Optional[Union[str, List[str]]] = None,
                          copy_numbers: Optional[Union[int, List[int]]] = None):
        """Compare CNV loci between two patients."""
        ranges_1 = patient_1.get_affected_ranges(chromosomes, copy_numbers)
        ranges_2 = patient_2.get_affected_ranges(chromosomes, copy_numbers)

        total_union = 0
        total_intersect = 0

        chromosomes = set(ranges_1.keys()) | set(ranges_2.keys())
        for chrom in chromosomes:
            if chrom not in ranges_1:
                ranges_1[chrom] = []
            if chrom not in ranges_2:
                ranges_2[chrom] = []

            ranges = ranges_1[chrom] + ranges_2[chrom]
            total_union += sum([len(x) for x in merge_range_list(ranges)])
            total_intersect += length_of_range_intersects(ranges)

        jaccard_index = jaccard(total_intersect, total_union)
        return jaccard_index, total_intersect

    @staticmethod
    def compare_pair_genes(patient_1: Patient, patient_2: Patient,
                           chromosomes: Optional[Union[str, List[str]]] = None,
                           copy_numbers: Optional[Union[int, List[int]]] = None):
        """Compare affected genes between two patients."""
        gene_info = [patient.get_all_genes(chromosomes, copy_numbers)
                     for patient in [patient_1, patient_2]]
        jaccard_index, intersect = jaccard(*gene_info)
        return jaccard_index, intersect

    @staticmethod
    def compare_pair_HI_genes(patient_1: Patient, patient_2: Patient,
                              chromosomes: Optional[Union[str, List[str]]] = None,
                              copy_numbers: Optional[Union[int, List[int]]] = None,
                              pLI_threshold=0.9,
                              HI_threshold=10,
                              phaplo_threshold=0.86,
                              mode="confirm"):
        """Compare affected HI genes between two patients. """
        gene_info = [patient.get_all_HI_genes(chromosomes, copy_numbers, pLI_threshold,
                                              HI_threshold, phaplo_threshold, mode)
                     for patient in [patient_1, patient_2]]
        jaccard_index, intersect = jaccard(*gene_info)
        return jaccard_index, intersect

    @staticmethod
    def compare_pair_dominant_genes(patient_1: Patient, patient_2: Patient,
                                    chromosomes: Optional[Union[str, List[str]]] = None,
                                    copy_numbers: Optional[Union[int, List[int]]] = None):
        """Compare affected dominant-effect genes between two patients."""
        gene_info = [patient.get_all_dominant_genes(chromosomes, copy_numbers)
                     for patient in [patient_1, patient_2]]
        jaccard_index, intersect = jaccard(*gene_info)
        return jaccard_index, intersect

    @staticmethod
    def compare_pair_dominant_match(patient_1: Patient, patient_2: Patient,
                                    chromosomes: Optional[Union[str, List[str]]] = None,
                                    copy_numbers: Optional[Union[int, List[int]]] = None) -> bool:
        """Check if two patients share the same affected dominant-effect genes."""
        match = (patient_1.get_all_dominant_genes(chromosomes, copy_numbers)
                 == patient_2.get_all_dominant_genes(chromosomes, copy_numbers))
        return match

    @staticmethod
    def compare_pair_hpos(patient_1: Patient, patient_2: Patient):
        """Compare HPO terms between two patients."""
        hpo_set1 = {hpo for hpo, response in patient_1.phenotypes.items() if response is True}
        hpo_set2 = {hpo for hpo, response in patient_2.phenotypes.items() if response is True}
        jaccard_index, intersect = jaccard(hpo_set1, hpo_set2)
        return jaccard_index, intersect

    # %% Phenotypes prevalence and prediction methods ==========

    def filter_patient_intersects(self, patient_id,
                                  length_similarity=0, loci_similarity=0,
                                  gene_similarity=0, hi_gene_similarity=0,
                                  dom_gene_match=True, hpo_similarity=0,
                                  include_self=False, as_patient_database=False):
        """Retrieve patients that satisfy comparison criteria with a selected patient.

        All patients whose comparison values with the specified patient are at least
        as high as the specified criteria (or evaluate to True, for dom_gene_match)
        are returned. Length, loci, gene, hi_gene, and hpo similarity are Jaccard
        index values from 0 to 1. Patients can be returned as a dictionary of
        'patient.id': patient, or as a PatientDatabase object. Specifying
        'include_self=True' will include the original patient in the results.
        """
        intersects = [inter for inter in self.lookup(patient_id, "all") if all(
            [
                inter.length_similarity >= length_similarity,
                inter.loci_similarity >= loci_similarity,
                inter.gene_similarity >= gene_similarity,
                inter.hi_gene_similarity >= hi_gene_similarity,
                inter.dom_gene_match or not dom_gene_match,
                inter.hpo_similarity >= hpo_similarity,
                not inter.self_compare
                ]
            )]
        if include_self:
            intersects.append(self.lookup(patient_id))
        if as_patient_database:
            patients = [inter.get_other_patient(patient_id) for inter in intersects]
            patients = {patient.id: patient for patient in patients}
            patients = PatientDatabase(patients)
            return patients
        return intersects

    def calculate_patient_pheno_prevalence(self, patient_id, phenotypes,
                                           length_similarity=0, loci_similarity=0,
                                           gene_similarity=0, hi_gene_similarity=0,
                                           dom_gene_match=True, hpo_similarity=0):
        """Calculate phenotype prevalence in a patient subgroup."""
        group = self.filter_patient_intersects(
            patient_id, length_similarity, loci_similarity,
            gene_similarity, hi_gene_similarity, dom_gene_match, hpo_similarity,
            include_self=True, as_patient_database=True
            )
        size = group.size
        prevs = {
            phenotype: sum([patient.phenotypes[phenotype] == "T" for patient in group
                            if phenotype in patient.phenotypes])
            for phenotype in phenotypes
            }
        prevs = {phenotype: PhenotypePrevalence(patient_id, group, phenotype, size, prevalence)
                 for phenotype, prevalence in prevs.items()}
        prevs = PatientGroupPrevalences(patient_id, group, prevs)
        return prevs

    def calculate_all_patient_pheno_prevalences(self, phenotypes, length_similarity=0,
                                                loci_similarity=0, gene_similarity=0,
                                                hi_gene_similarity=0, dom_gene_match=True,
                                                hpo_similarity=0):
        """Calculate phenotype prevalences for all patient subgroups."""
        params = [phenotypes, length_similarity, loci_similarity, gene_similarity,
                  hi_gene_similarity, dom_gene_match, hpo_similarity]
        all_prevs = {patient: self.calculate_patient_pheno_prevalence(patient, *params)
                     for patient in self.index.keys()}
        all_prevs = HomogeneityDatabase(all_prevs)
        return all_prevs

    @staticmethod
    def predict_group_phenotypes(comparison_group, show=10,
                                 additional_hpos: set = None,
                                 additional_groupname="Added"):
        """Predict phenotypes for a patient subgroup."""
        population = len(comparison_group)
        if population == 0:
            return dict()

        if additional_hpos is None:
            additional_hpos = set()
        group_terms = {term for patient in comparison_group
                       for term in patient.phenotypes.keys()}
        group_terms -= additional_hpos

        predictions = dict()
        for group, term_set in zip(["predicted", additional_groupname],
                                   [group_terms, additional_hpos]):
            for term in term_set:
                count = Counter([patient.phenotypes[term]
                                 for patient in comparison_group
                                 if term in patient.phenotypes])
                if count[True] == 0 and group == "predicted":
                    continue
                na_count = sum(count.values()) - (count[True] + count[False])
                prediction = TraitPrediction(term, population,
                                             true=count[True], false=count[False],
                                             na=na_count,
                                             group=group)
                predictions[term] = prediction

        return predictions

    def test_patient_pheno_predictions(self, patient_id,
                                       length_similarity=0, loci_similarity=0,
                                       gene_similarity=0, hi_gene_similarity=0,
                                       dom_gene_match=True, hpo_similarity=0,
                                       skip_no_hpos=True):
        """Predict phenotypes for a patient and compare with their known phenotypes."""
        params = [patient_id, length_similarity, loci_similarity, gene_similarity,
                  hi_gene_similarity, dom_gene_match, hpo_similarity]

        index_hpos = {hpo for hpo, response in self.patient_db[patient_id].phenotypes.items()
                      if response is True}

        comparison_group = self.filter_patient_intersects(*params, include_self=False,
                                                          as_patient_database=True)
        comparison_group = list(comparison_group.patients.values())

        if skip_no_hpos:
            comparison_group = [patient for patient in comparison_group
                                if patient.phenotypes]

        predictions = self.predict_group_phenotypes(comparison_group, show=0,
                                                    additional_hpos=index_hpos,
                                                    additional_groupname="index")
        predictions = PatientPredictions(
            patient=patient_id,
            patient_group=PatientDatabase({patient.id: patient for patient in comparison_group}),
            predictions=predictions
            )
        return predictions

    def test_all_patient_pheno_predictions(self, length_similarity=0, loci_similarity=0,
                                           gene_similarity=0, hi_gene_similarity=0,
                                           dom_gene_match=True, hpo_similarity=0,
                                           skip_no_hpos=True, filter_unknowns=True,
                                           **kwargs):
        """Compare phenotypes for each patient and compare with their known phenotypes."""
        params = [length_similarity, loci_similarity, gene_similarity, hi_gene_similarity,
                  dom_gene_match, hpo_similarity, skip_no_hpos]

        all_predictions = {patient_id: self.test_patient_pheno_predictions(patient_id, *params)
                           for patient_id in self.index}

        if filter_unknowns:
            all_predictions = {patient_id: results
                               for patient_id, results in all_predictions.items()
                               if len(results) > 0}

        all_predictions = PredictionDatabase(all_predictions)
        return all_predictions

    def tabulate_summary(self):
        table = []
        for comparison in self:
            entry = list(comparison.patients.keys())
            entry.extend([
                comparison.length_similarity,
                comparison.loci_similarity,
                comparison.gene_similarity,
                comparison.hi_gene_similarity,
                comparison.dom_gene_match
                ])
            table.append(entry)
        table = pd.DataFrame(table)
        return table


class PatientIntersect:
    """Record of similarity between two patients."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, patient_1, patient_2, copy_number,
                 length_compare, loci_compare,
                 gene_compare, hi_gene_compare,
                 dom_gene_compare, dom_gene_match,
                 hpo_compare):
        self.patients = {patient_1.id: patient_1, patient_2.id: patient_2}
        self.ids = {patient_1.id, patient_2.id}
        self.copy_number = copy_number

        self.length_similarity = length_compare

        self.loci_similarity = loci_compare[0]
        self.loci_shared_size = loci_compare[1]

        self.gene_similarity = gene_compare[0]
        self.genes = gene_compare[1]
        self.gene_count = len(gene_compare[1])

        self.hi_gene_similarity = hi_gene_compare[0]
        self.hi_genes = hi_gene_compare[1]
        self.hi_gene_count = len(hi_gene_compare[1])

        self.dom_gene_similarity = dom_gene_compare[0]
        self.dom_genes = dom_gene_compare[1]
        self.dom_gene_count = len(dom_gene_compare[1])
        self.dom_gene_match = dom_gene_match

        self.hpo_similarity = hpo_compare[0]
        self.hpos = hpo_compare[1]
        self.hpo_count = len(hpo_compare[1])

        self.self_compare = False
        if len(self.ids) == 1:
            self.self_compare = True

    def __repr__(self):
        """Get official string representation."""
        string = (f"PatientIntersect(patients=[{', '.join(self.patients.keys())}], "
                  f"length_similarity={self.length_similarity}, "
                  f"loci_similarity={self.loci_similarity}, "
                  f"gene_similarity={self.gene_similarity}, "
                  f"hi_gene_similarity={self.hi_gene_similarity}, "
                  f"dom_gene_similarity={self.dom_gene_similarity}, "
                  f"hpo_similarity={self.hpo_similarity})")
        return string

    def __str__(self):
        """Get pretty-printing string representation."""
        string = (f"Similarities of {' vs. '.join(self.patients.keys())}:\n"
                  f"  Genes: {self.gene_similarity}\n"
                  f"  HPO terms: {self.hpo_similarity}\n"
                  f"  Length: {self.length_similarity}\n"
                  f"  Position: {self.loci_similarity}")
        return string

    def gene_info(self):
        """Get gene portion of intersect."""
        return self.genes, self.gene_count, self.gene_similarity

    def hpo_info(self):
        """Get HPO portion of intersect."""
        return self.hpos, self.hpo_count, self.hpo_similarity

    def get_other_id(self, ID):
        if ID not in self.patients:
            raise IndexError("This ID not present.")
        if set(self.patients.keys()) == {ID}:
            return ID
        other_id = list(set(self.patients.keys()) - {ID})[0]
        return other_id

    def get_other_patient(self, ID):
        other_id = self.get_other_id(ID)
        return self.patients[other_id]

    def get_similarities(self):
        similarities = [self.length_similarity, self.loci_similarity, self.gene_similarity,
                        self.hi_gene_similarity, self.dom_gene_match, self.hpo_similarity]
        return similarities


def predict_phenotypes_for_patient(patient: Patient, comparison_database: ComparisonTable,
                                   chromosomes: Optional[Union[str, List[str]]] = None,
                                   copy_numbers: Optional[Union[int, List[int]]] = None,
                                   length_similarity=0, loci_similarity=0,
                                   gene_similarity=0, hi_gene_similarity=0,
                                   dom_gene_match=True, hpo_similarity=0,
                                   pLI_threshold=0.9, HI_threshold=10,
                                   phaplo_threshold=0.86, mode="confirm",
                                   tabulate=False):
    """Predict phenotypes for a patient objects with CNVs."""
    params = dict(chromosomes=chromosomes, copy_numbers=copy_numbers,
                  pLI_threshold=pLI_threshold, HI_threshold=HI_threshold,
                  phaplo_threshold=phaplo_threshold, mode=mode, save_results=False)
    intersections = comparison_database.compare_patient_vs_others(patient, **params)
    intersections = [inter for inter in intersections
                     if all([inter.length_similarity >= length_similarity,
                             inter.loci_similarity >= loci_similarity,
                             inter.gene_similarity >= gene_similarity,
                             inter.hi_gene_similarity >= hi_gene_similarity,
                             inter.dom_gene_match or not dom_gene_match,
                             inter.hpo_similarity >= hpo_similarity])]
    group = [inter.get_other_patient(patient.id) for inter in intersections]
    group = PatientDatabase({group_patient.id: group_patient for group_patient in group})
    predictions = comparison_database.predict_group_phenotypes(group, show=0)
    predictions = PatientPredictions(patient.id, group, predictions)
    if tabulate:
        predictions = predictions.convert_patient_predictions_to_df()
        return predictions
    return predictions


def predict_phenotypes_for_cnvs(cnvs: list[CNV],
                                comparison_database: ComparisonTable,
                                chromosomes: Optional[Union[str, List[str]]] = None,
                                copy_numbers: Optional[Union[int, List[int]]] = None,
                                length_similarity=0, loci_similarity=0,
                                gene_similarity=0, hi_gene_similarity=0,
                                dom_gene_match=True, hpo_similarity=0,
                                pLI_threshold=0.9, HI_threshold=10,
                                phaplo_threshold=0.86, mode="confirm",
                                tabulate=False):
    """Predict phenotypes based on a list of CNV objects.

    CNVs must already have affected genes assigned."""
    comparison_params = dict(
        comparison_database=comparison_database,
        chromosomes=chromosomes,
        copy_numbers=copy_numbers,
        length_similarity=length_similarity, loci_similarity=loci_similarity,
        gene_similarity=gene_similarity, hi_gene_similarity=hi_gene_similarity,
        dom_gene_match=dom_gene_match, hpo_similarity=hpo_similarity,
        pLI_threshold=pLI_threshold, HI_threshold=HI_threshold,
        phaplo_threshold=phaplo_threshold, mode=mode,
        tabulate=tabulate
        )
    patient = Patient("query")
    patient.cnvs = cnvs
    predictions = predict_phenotypes_for_patient(patient, **comparison_params)
    return predictions


def _convert_cnv_str_to_cnv(cnv_string: str):
    """Convert string format CNV ('chr:start-stop:copy_number') to CNV object ."""
    cnv = cnv_string.split(":")
    cnv = cnv[:1] + [int(n) for n in cnv[1].split("-")] + [int(cnv[-1])]
    cnv = CNV(*cnv, ID="query")
    return cnv


def predict_phenotypes_for_cnv_strings(cnv_strings: list[str],
                                       geneset: GeneSet,
                                       comparison_database: ComparisonTable,
                                       chromosomes: Optional[Union[str, List[str]]] = None,
                                       copy_numbers: Optional[Union[int, List[int]]] = None,
                                       length_similarity=0, loci_similarity=0,
                                       gene_similarity=0, hi_gene_similarity=0,
                                       dom_gene_match=True, hpo_similarity=0,
                                       pLI_threshold=0.9, HI_threshold=10,
                                       phaplo_threshold=0.86, mode="confirm",
                                       tabulate=False, **kwargs):
    """Predict phenotypes based on a list of strings containing CNV loci and copy numbers.

    CNV string format is 'chr:start-stop:copy_number'. CNV objects will be created
    automatically, and affected genes will be assigned."""
    cnvs = [_convert_cnv_str_to_cnv(cnv_string) for cnv_string in cnv_strings]
    for cnv in cnvs:
        cnv.set_genes(geneset)
    comparison_params = dict(
        comparison_database=comparison_database,
        chromosomes=chromosomes,
        copy_numbers=copy_numbers,
        length_similarity=length_similarity, loci_similarity=loci_similarity,
        gene_similarity=gene_similarity, hi_gene_similarity=hi_gene_similarity,
        dom_gene_match=dom_gene_match, hpo_similarity=hpo_similarity,
        pLI_threshold=pLI_threshold, HI_threshold=HI_threshold,
        phaplo_threshold=phaplo_threshold, mode=mode,
        tabulate=tabulate
        )
    predictions = predict_phenotypes_for_cnvs(cnvs, **comparison_params)
    return predictions
