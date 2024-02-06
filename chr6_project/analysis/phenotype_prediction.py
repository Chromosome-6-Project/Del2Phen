"""The Chromosome 6 Project - Phenotype Prediction.

@author: T.D. Medina
"""

import numpy as np
import pandas as pd
from pronto import Term

from chr6_project.analysis.plotting import plot_precision_stats


class TraitPrediction:
    def __init__(self, trait, population, true, false,
                 unsure, na, group=None):
        self.trait = trait
        self.group = group
        self.population = population
        self.true_count = true
        self.false_count = false
        self.unsure_count = unsure
        self.na_count = na
        self._found = False

        if self.population == 0:
            self.freq = 0
        else:
            self.freq = self.true_count/self.population

        self.population_adjusted = self.true_count + self.false_count
        if self.population_adjusted == 0:
            self.freq_adjusted = 0
        else:
            self.freq_adjusted = self.true_count/self.population_adjusted

        if self.true_count > 0:
            self._found = True

    def __str__(self):
        string = (f"{self.trait.name}\t{self.population}\t{self.population_adjusted}\t"
                  f"{self.freq}\t{self.freq_adjusted}")
        return string

    def __repr__(self):
        string = "TraitPrediction("
        attrs = []
        for attr, value in self.__dict__.items():
            if isinstance(value, float):
                value = round(value, 3)
            attrs.append(f"{attr}={value}")
        string = string + ", ".join(attrs) + ")"
        return string

    def make_table_row(self):
        row = {}
        for attr, value in self.__dict__.items():
            if attr.startswith("_"):
                continue
            if isinstance(value, Term):
                row[attr] = value.id
                row["trait name"] = value.name
                continue
            row[attr] = value
        return row


class PatientPredictions:
    def __init__(self, patient, patient_group, predictions):
        self.patient = patient
        self.patient_group = patient_group
        self.predictions = predictions

    def __len__(self):
        return len(self.predictions)

    def make_confusion_matrix(self, phenotypes=None, rel_threshold=0.2, abs_threshold=2,
                              use_adjusted_frequency=True):
        predictions = self.predictions
        if phenotypes is not None:
            predictions = {term: prediction for term, prediction in predictions.items()
                           if term in phenotypes}
        confusion = np.zeros([2, 2])

        index = {term for term, prediction in predictions.items()
                 if prediction.group == "index"}
        if use_adjusted_frequency:
            above_threshold = {term for term, prediction in predictions.items()
                               if prediction.true_count >= abs_threshold
                               and prediction.freq_adjusted >= rel_threshold}
        else:
            above_threshold = {term for term, prediction in predictions.items()
                               if prediction.true_count >= abs_threshold
                               and prediction.freq >= rel_threshold}
        confusion[0][0] = len(index & above_threshold)
        confusion[1][0] = len(above_threshold - index)
        confusion[0][1] = len(index - above_threshold)
        confusion[1][1] = len(set(predictions.keys()) - (above_threshold | index))

        return confusion

    def convert_patient_predictions_to_df(self):
        table = [prediction.make_table_row() for prediction in self.predictions.values()]
        column_labels = [key.title() for key in table[0]]
        table = pd.DataFrame(table)
        table.columns = column_labels
        return table


class PredictionDatabase:
    def __init__(self, patient_predictions):
        self.predictions = patient_predictions

    def make_overall_confusion_matrix(self, phenotypes=None, rel_threshold=0.2,
                                      abs_threshold=2, use_adjusted_frequency=True,
                                      group_size_threshold=5):
        params = (phenotypes, rel_threshold, abs_threshold, use_adjusted_frequency)
        matrices = []
        for prediction in self.predictions.values():
            if prediction.patient_group.size < group_size_threshold:
                continue
            matrices.append(prediction.make_confusion_matrix(*params))
        matrix = sum(matrices)
        return matrix

    def make_individual_confusion_matrices(self, phenotypes=None, rel_threshold=0.2,
                                           abs_threshold=2, use_adjusted_frequency=True,
                                           group_size_threshold=5):
        params = (phenotypes, rel_threshold, abs_threshold, use_adjusted_frequency)
        matrices = {
            patient: predict.make_confusion_matrix(*params)
            for patient, predict in self.predictions.items()
            if predict.patient_group.size >= group_size_threshold}
        return matrices

    def calculate_individual_precision(self, phenotypes=None, rel_threshold=0.2,
                                       abs_threshold=2, use_adjusted_frequency=True,
                                       group_size_threshold=5):
        params = (phenotypes, rel_threshold, abs_threshold, use_adjusted_frequency,
                  group_size_threshold)
        stat_names = ("Sensitivity", "Specificity", "PPV", "NPV")
        precision_stats = {}
        for patient, mat in self.make_individual_confusion_matrices(*params).items():
            margins = list(mat.sum(1)) + list(mat.sum(0))
            # The mod here just makes it convenient to divide the relevant matrix cells
            # by their respective margins from the list above.
            patient_stats = [cell / margins[i] if (cell := mat[i % 2][i % 2]) > 0 else 0
                             for i in range(4)]
            precision_stats[patient] = dict(zip(stat_names, patient_stats))
        return precision_stats

    def make_individual_precision_table(self, phenotypes=None, rel_threshold=0.2,
                                        abs_threshold=2, use_adjusted_frequency=True,
                                        group_size_threshold=5):
        precisions = self.calculate_individual_precision(phenotypes, rel_threshold,
                                                         abs_threshold, use_adjusted_frequency,
                                                         group_size_threshold)
        for pid in precisions.keys():
            precisions[pid]["Group_Size"] = self.predictions[pid].group_size
        table = pd.DataFrame.from_dict(precisions, orient="index")
        return table

    def calculate_overall_precision(self, phenotypes=None, rel_threshold=0.2,
                                    abs_threshold=2, use_adjusted_frequency=True,
                                    group_size_threshold=5):
        params = (phenotypes, rel_threshold, abs_threshold, group_size_threshold,
                  use_adjusted_frequency)
        confusion = self.make_overall_confusion_matrix(*params)
        stat_names = ("Sensitivity", "Specificity", "PPV", "NPV")
        margins = list(confusion.sum(1)) + list(confusion.sum(0))
        stats = [cell / margins[i] if (cell := confusion[i % 2][i % 2]) > 0 else 0
                 for i in range(4)]
        stats = dict(zip(stat_names, stats))
        return stats

    def plot_precision_stats(self, patient_database, phenotypes=None, rel_threshold=0.2,
                             abs_threshold=2, use_adjusted_frequency=True,
                             group_size_threshold=5):
        precision_plot = plot_precision_stats(self, patient_database, phenotypes,
                                              rel_threshold, abs_threshold,
                                              use_adjusted_frequency,
                                              group_size_threshold)
        return precision_plot


def predict_test(comparison, patient_list="C:/Users/Ty/Documents/Chr6/Predict_tests/test_patients.txt"):
    with open(patient_list) as infile:
        aafkes_patients = infile.readlines()
    aafkes_patients = [x.strip() for x in aafkes_patients]
    all_tests = comparison.test_all_patient_pheno_predictions(gene_similarity=.7)
    aafkes_tests = {x: y for x, y in all_tests.items() if x in aafkes_patients}
    return all_tests, aafkes_tests
