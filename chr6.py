#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Automatic Phenotype Prediction.

@author: T.D. Medina
"""

import bz2
import csv
from collections import namedtuple
from math import sqrt
import pickle

import numpy as np
import pandas as pd
from pronto import Ontology, Term

from data_objects import Patient, PatientDatabase, GenomeDict
import gene_set
from phenotype_homogeneity import (
    PhenotypeHomogeneity,
    GroupPhenotypeHomogeneity,
    phenotype_homo_test
    )
from utilities import jaccard, length_of_range_intersects, merge_range_list


class DataManager:
    """Collection of Chromosome 6 Project data handling methods."""

    def __init__(self):
        pass

    @classmethod
    def read_data(cls, data_file):
        """Import general CSV data."""
        data = []
        with open(data_file, encoding="utf-8") as infile:
            reader = csv.DictReader(infile, delimiter=",", quotechar='"')
            for row in reader:
                for k, v in row.items():
                    if not v:
                        row[k] = None
                    else:
                        row[k] = v.replace("\n", " ")
                        row[k] = row[k].replace("\t", "")
                data.append(row)
        data = cls.coerce_booleans(data)
        return data

    @staticmethod
    def write_table(data, filename):
        """Write CSV data."""
        writer = ["\t".join(data[0].keys()) + "\n"]
        for entry in data:
            writer.append("\t".join([str(value) for value in entry.values()])
                          + "\n")
        with open(filename, "w") as outfile:
            outfile.writelines(writer)

    @staticmethod
    def coerce_booleans(data):
        """Replace True/False and Yes/No strings with bools."""
        coerced = data
        for entry in coerced:
            for k, v in entry.items():
                if not v:
                    continue
                if v.lower() in ["true", "yes"]:
                    entry[k] = True
                elif v.lower() in ["false", "no"]:
                    entry[k] = False
        return coerced

    @staticmethod
    def fix_hg19_positions(data):
        """Fill missing position data, if available, and convert to int."""
        fixed = data
        start_report = "Start position in report"
        start_convert = "Start positie in Hg19"
        stop_report = "Stop position in report"
        stop_convert = "Stop positie in Hg19"

        for entry in fixed:
            build = entry["Human genome build version"]
            for reported, converted in [(start_report, start_convert),
                                        (stop_report, stop_convert)]:
                reported_value = entry[reported]
                converted_value = entry[converted]

                # Try to fill missing if an hg19 position was already reported.
                if not converted_value:
                    # If both positions are missing, skip this record.
                    if not reported_value:
                        continue
                    # If not build hg19, skip this record.
                    if build != "hg19/NCBI37":
                        continue
                    entry[converted] = reported_value

                # Clean up value and convert to integer.
                try:
                    new_int = entry[converted].replace(" ", "")
                    new_int = new_int.replace(",", "")
                    new_int = int(new_int)
                    entry[converted] = new_int
                except ValueError as err:
                    print(f"Could not parse {err}")
        return fixed

    @staticmethod
    def trim_chromosome_names(data):
        """Remove 'chr' from contig names."""
        trimmed = data
        for entry in trimmed:
            if not entry["Chromosome"]:
                continue
            if entry["Chromosome"].lower().startswith("chr"):
                entry["Chromosome"] = entry["Chromosome"].strip("chrCHR")
        return trimmed

    @staticmethod
    def fix_patient_hpos(data):
        """Change Molgenis-styled HPO terms to OBO HPO terms."""
        # old_keys = list(data[0].keys())
        new_data = {}
        # new_data = {patient["label"]: patient for patient in data}
        for patient in data:
            new_patient = dict(patient)
            for key in patient:
                new_patient[key.replace("_", ":")] = new_patient.pop(key)
            new_data[new_patient["label"]] = new_patient
        return new_data

    @staticmethod
    def fix_patient_hpos2(data):
        """Change Molgenis-styled HPO terms to OBO HPO terms."""
        new_data = {
            entry["label"]: {x.replace("_", ":"): y
                             for x, y in list(entry.items())[1:]}
            for entry in data
            }
        return new_data

    @classmethod
    def fix_genotype_data(cls, data):
        """Apply all data fixes."""
        fixed = cls.trim_chromosome_names(data)
        fixed = cls.fix_hg19_positions(fixed)
        return fixed

    @staticmethod
    def make_patients(genotypes, phenotypes, geneset=None, hpos=None,
                      ontology=None, expand_hpos=True):
        """Construct Patient objects from data."""
        subjects = ({record["subject"] for record in genotypes}
                    | {record["owner"] for record in phenotypes})
        patients = {subject: Patient(subject) for subject in subjects}
        for genotype in genotypes:
            patients[genotype["subject"]].genotypes.append(genotype)
        for phenotype in phenotypes:
            patients[phenotype["owner"]].phenotypes.append(phenotype)
        if geneset:
            for patient in patients.values():
                patient.cnvs = patient.extract_cnvs()
                patient.identify_gene_overlaps(geneset)
        if ontology and hpos:
            for patient in patients.values():
                if patient.id not in hpos:
                    continue
                patient.hpo = {ontology[hpo_id]: response
                               for hpo_id, response in hpos[patient.id].items()}
                # patient_hpos = []
                # for hpo_id, value in hpos[patient.id].items():
                #     if hpo_id == "label":
                #         continue
                #     if value == "T":
                #         patient_hpos.append(ontology[hpo_id])
                # patient.hpo = patient_hpos
                if expand_hpos:
                    patient.expand_hpo_terms(ontology)
        for patient in patients.values():
            if patient.genotypes:
                patient.origin = patient.genotypes[0]["origin info"]
            patient.convert_birthday_to_datetime()
        return patients

    @staticmethod
    def print_summary_counts(patients):
        """Print a summary of genotype and phenotype counts."""
        size = len(patients.values())
        sums = [(len(patient.genotypes), len(patient.phenotypes))
                for patient in patients.values()]
        print(f"Patients: {size}")
        print("---------------------")
        print(f"{'Genotypes':13}{'Phenotypes':14}Counts")
        for geno, pheno in sorted(list(set(sums)), key=lambda x: f"{x[0]}_{x[1]}"):
            print(f"{geno:<13}{pheno:<14}{sums.count((geno,pheno)):3}")

    @staticmethod
    def make_ucsc_browser_tracks(patients, out, filter_chrom=None):
        """Write UCSC BED file of patient CNVs."""
        writer = ["browser hide all position chr6\n"]
        patient_list = sorted(
            [patient for patient in patients.values() if patient.cnvs],
            key=lambda x: x.cnvs[0].range.start
            )
        for i, patient in enumerate(patient_list):
            patient_writer = [f"track name=Patient_{i} visibility=2\n"
                              f"#chrom\tchromStart\tchromEnd\tname\n"]
            for cnv in patient.cnvs:
                patient_writer.append(
                    f"chr{cnv.chromosome}\t{cnv.range.start}\t"
                    f"{cnv.range.stop - 1}\t{str(cnv.change).replace(' ', '')}\n")
            writer += patient_writer
        with open(out, "w") as outfile:
            outfile.writelines(writer)

    @staticmethod
    def read_patient_drop_list(drop_list_file):
        with open(drop_list_file) as infile:
            drop_list = infile.readlines()
        drop_list = [name.strip() for name in drop_list]
        drop_list = {name for name in drop_list if not name.startswith("#")}
        return drop_list

    @staticmethod
    def filter_patients(patient_dict, drop_list, remove_ungenotyped=True,
                        remove_unphenotyped=True):
        patients_filtered = {name: patient for name, patient in patient_dict.items()
                             if name not in drop_list
                             and not (remove_ungenotyped and not patient.genotypes)
                             and not (remove_unphenotyped and not patient.phenotypes)}
        return patients_filtered


# %% Comparisons
# TODO: Add a method to add a patient.
class ComparisonTable:
    """Data object holding all patient vs. patient comparisons."""

    def __init__(self, patient_db=None, comparison_table=None):
        if comparison_table is not None:
            self.read_from_existing(comparison_table)
            return
        self.patient_db = patient_db
        self.raw = self.compare_patients()
        self.index = self.make_index()
        self.array = self.make_array()
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

    def make_array(self):
        """Convert raw comparison dictionary to numpy array."""
        array = []
        for patient1 in self.index:
            values = []
            for patient2 in self.index:
                if patient2 in self.raw[patient1]:
                    values.append(self.raw[patient1][patient2])
                else:
                    values.append(self.raw[patient2][patient1])
            array.append(values)
        array = np.array(array)
        return array

    def make_index(self):
        """Create name-to-number mapping to look up names in numpy array."""
        return {j: i for i, j in enumerate(self.raw)}

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
        self.raw = comparison_table.raw
        self.index = comparison_table.index
        self.array = comparison_table.array
        self.size = comparison_table.size

    def compare_patients(self):
        """Compare all patients to each other."""
        ids = list(self.patient_db.patients.keys())
        comparisons = {}
        while ids:
            id_i = ids.pop()
            patient_i = self.patient_db[id_i]
            patient_comparison = {}
            patient_comparison[id_i] = self.compare_all(patient_i, patient_i)

            for id_j in ids:
                patient_j = self.patient_db[id_j]
                patient_comparison[id_j] = self.compare_all(patient_i, patient_j)
            comparisons[id_i] = patient_comparison
        return comparisons

    @classmethod
    def compare_all(cls, patient_1, patient_2):
        """Compare all metrics between two patients."""
        cnv_type = {cnv.change for cnv in patient_1.cnvs + patient_2.cnvs}
        if len(cnv_type) == 0:
            cnv_type = "N/A"
        elif len(cnv_type) > 1:
            cnv_type = "Mixed"
        else:
            cnv_type = list(cnv_type)[0]
        length_compare = cls.compare_length(patient_1, patient_2)
        loci_compare = cls.compare_loci(patient_1, patient_2)
        gene_compare = cls.compare_genes(patient_1, patient_2)
        hi_compare = cls.compare_HI_genes(patient_1, patient_2)
        hpo_compare = cls.compare_hpos(patient_1, patient_2)
        comparison = PatientIntersect(
            patient_1, patient_2, cnv_type,
            length_compare, loci_compare, gene_compare, hi_compare, hpo_compare
            )
        return comparison

    @staticmethod
    def compare_length(patient_1, patient_2):
        """Compare CNV length between two patients."""
        size_1 = sum([cnv.length for cnv in patient_1.cnvs])
        size_2 = sum([cnv.length for cnv in patient_2.cnvs])
        if size_1 == 0 or size_2 == 0:
            return 0
        return jaccard(size_1, size_2)

    @staticmethod
    def compare_loci(patient_1, patient_2):
        """Compare CNV loci between two patients."""
        ranges_1 = patient_1.get_affected_ranges()
        ranges_2 = patient_2.get_affected_ranges()

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
    def compare_genes(patient_1, patient_2):
        """Compare affected genes between two patients."""
        jaccard_index, intersect = jaccard(patient_1.all_genes(), patient_2.all_genes())
        return jaccard_index, intersect

    @staticmethod
    def compare_HI_genes(patient_1, patient_2,
                         pLI_threshold=0.9, HI_threshold=10):
        """Compare affected HI genes between two patients."""
        jaccard_index, intersect = jaccard(
            patient_1.all_HI_genes(pLI_threshold, HI_threshold),
            patient_2.all_HI_genes(pLI_threshold, HI_threshold)
            )
        return jaccard_index, intersect

    @staticmethod
    def compare_hpos(patient_1, patient_2):
        """Compare HPO terms between two patients."""
        hpo_set1 = {hpo for hpo, response in patient_1.hpo.items() if response == "T"}
        hpo_set2 = {hpo for hpo, response in patient_2.hpo.items() if response == "T"}
        jaccard_index, intersect = jaccard(hpo_set1, hpo_set2)
        return jaccard_index, intersect

    # TODO: Needs to be updated, won't work currently.
    def write_all_comparisons(self, outfile):
        """Write comparison results to TSV file."""
        properties = ["length_similarity",
                      "loci_similarity", "loci_shared_size",
                      "gene_similarity", "gene_count",
                      "hpo_similarity", "hpo_count"]
        write_me = ["\t".join(["patient1", "patient2"] + properties) + "\n"]
        for intersect in self:
            p1, p2 = intersect.patients
            if p1 == p2:
                continue
            this_intersect = "\t".join(
                [f"{intersect.patients[0].id}", f"{intersect.patients[1].id}"]
                + [f"{intersect.__getattribute__(prop)}" for prop in properties]
                ) + "\n"
            write_me.append(this_intersect)
        with open(outfile, "w") as out:
            out.writelines(write_me)

    def filter_patient_comparisons(self, patient_id, length_similarity=0,
                                   loci_similarity=0, gene_similarity=0,
                                   hi_similarity=0, hpo_similarity=0):
        patient1 = self.patient_db[patient_id]
        intersections = self.lookup(patient_id, "all")
        filtered = []
        for intersect in intersections:
            if intersect.patients[0] == intersect.patients[1]:
                continue

            if patient1 != intersect.patients[0]:
                patient2 = intersect.patients[0]
            else:
                patient2 = intersect.patients[1]

            if not patient2.hpo:
                continue

            if not all([intersect.length_similarity >= length_similarity,
                        intersect.loci_similarity >= loci_similarity,
                        intersect.gene_similarity >= gene_similarity,
                        intersect.hi_gene_similarity >= hi_similarity,
                        intersect.hpo_similarity >= hpo_similarity]):
                continue

            filtered.append(patient2)
        return filtered

    def filter_patient_intersections(self, patient_id, length_similarity=0,
                                     loci_similarity=0, gene_similarity=0,
                                     hi_similarity=0, hpo_similarity=0,
                                     include_self=False):
        intersections = {intersect.get_other_id(patient_id): intersect
                         for intersect in self.lookup(patient_id, "all")
                         if all([intersect.length_similarity >= length_similarity,
                                 intersect.loci_similarity >= loci_similarity,
                                 intersect.gene_similarity >= gene_similarity,
                                 intersect.hi_gene_similarity >= hi_similarity,
                                 intersect.hpo_similarity >= hpo_similarity])}
        if include_self is False and patient_id in intersections:
            del intersections[patient_id]
        return intersections

    def make_patient_intersection_group(self, patient_id, length_similarity=0,
                                        loci_similarity=0, gene_similarity=0,
                                        hi_similarity=0, hpo_similarity=0,
                                        include_self=True, as_patient_database=True):
        intersections = self.filter_patient_intersections(
            patient_id, length_similarity, loci_similarity, gene_similarity,
            hi_similarity, hpo_similarity, include_self
            )
        patients = [intersection.get_other_patient(patient_id)
                    for intersection in intersections.values()]
        if as_patient_database:
            patients = PatientDatabase({patient.id: patient for patient in patients})
        return patients

    def test_phenotype_homogeneities(self, patient_id, phenotypes, length_similarity=0,
                                     loci_similarity=0, gene_similarity=0,
                                     hi_similarity=0, hpo_similarity=0):
        group = self.make_patient_intersection_group(
            patient_id, length_similarity, loci_similarity,
            gene_similarity, hi_similarity, hpo_similarity,
            as_patient_database=True)
        size = group.size
        homogeneities = {
            phenotype: sum([patient.hpo[phenotype] == "T" for patient in group
                            if phenotype in patient.hpo])
            for phenotype in phenotypes
            }
        homogeneities = {phenotype: PhenotypeHomogeneity(patient_id, group, phenotype, size, prevalence)
                         for phenotype, prevalence in homogeneities.items()}
        homogeneities = GroupPhenotypeHomogeneity(patient_id, group, homogeneities)
        return homogeneities

    def test_all_homogeneities(self, phenotypes, length_similarity=0,
                               loci_similarity=0, gene_similarity=0,
                               hi_similarity=0, hpo_similarity=0,
                               group_size_threshold=5, include_isolates=False):
        params = [phenotypes, length_similarity, loci_similarity,
                  gene_similarity, hi_similarity, hpo_similarity]
        all_homogeneities = {patient: self.test_phenotype_homogeneities(patient, *params)
                             for patient in self.index.keys()}
        upper_homogeneities = {patient: homogen for patient, homogen in all_homogeneities.items()
                               if homogen.group_size >= group_size_threshold}
        lower_homogeneities = {patient: homogen for patient, homogen in all_homogeneities.items()
                               if homogen.group_size < group_size_threshold}

        return all_homogeneities, upper_homogeneities, lower_homogeneities

# =============================================================================
#     @staticmethod
#     def predict_phenotypes(comparison_group, show=10):
#         size = len(comparison_group)
#         if size == 0:
#             return dict()
#
#         all_hpo = {hpo: 0 for patient in comparison_group for hpo in patient.hpo}
#         for hpo in all_hpo:
#             for patient in comparison_group:
#                 if hpo in patient.hpo:
#                     all_hpo[hpo] += 1
#         all_hpo = {hpo: (count, count/size) for hpo, count in all_hpo.items()}
#
#         if show == 0:
#             return all_hpo
#         if show == "all":
#             show = len(all_hpo)
#         string = "Top {show} phenotypes out of {len(all_hpo)}:\n"
#         for hpo, (count, freq) in sorted(all_hpo.items(), key=lambda x: x[1][0], reverse=True)[:show]:
#             string += f"{hpo.name}:    {count}/{size}    {freq:.2%}\n"
#         print(string)
#
#         return all_hpo
# =============================================================================

    @staticmethod
    def predict_phenotypes2(comparison_group, show=10,
                            additional_hpos=None, additional_groupname="Added"):
        population = len(comparison_group)
        if population == 0:
            return dict()

        all_hpo = {hpo: {"t": 0, "f": 0, "unsure": 0, "na": 0, "group": "predicted"}
                   for patient in comparison_group for hpo in patient.hpo.keys()}
        if additional_hpos is not None:
            all_hpo.update({hpo: {"t": 0, "f": 0, "unsure": 0, "na": 0, "group": additional_groupname}
                            for hpo in additional_hpos})

        for hpo, counts in all_hpo.items():
            for patient in comparison_group:
                # TODO: This is to make expanded HPOs work. But how should
                # they be counted? This counts them as NAs.
                # if hpo not in patient.hpo:
                #     counts["na"] += 1
                response = patient.hpo[hpo].lower()
                counts[response] += 1

        all_hpo = {hpo: TraitPrediction(hpo, hpo.name, population, *counts.values())
                   for hpo, counts in all_hpo.items()
                   if counts["group"] == additional_groupname or counts["t"] > 0}

        if show == 0:
            return all_hpo

        if show == "all":
            show = len(all_hpo)

        print(f"\nTop {show} phenotypes out of {len(all_hpo)}:\n\n"
              "Trait\tPop\tPop.Adjust\tFreq\tFreq.Adjust")
        for trait_freq in sorted(all_hpo.values(), key=lambda x: x.true_count, reverse=True)[:show]:
            print(trait_freq)

        return all_hpo

# =============================================================================
#     def test_predictions(self, patient_id, freq_threshold=0,
#                          length_similarity=0, loci_similarity=0,
#                          gene_similarity=0, hpo_similarity=0):
#
#         comparison_group = self.filter_patient_comparisons(
#             patient_id,
#             length_similarity,
#             loci_similarity,
#             gene_similarity,
#             hpo_similarity
#             )
#
#         predictions = self.predict_phenotypes(comparison_group, show=0)
#         test_hpos = self.patient_db[patient_id].hpo
#
#         predicted_hpos = {hpo: PredictInfo(count, freq, hpo in test_hpos)
#                           for hpo, (count, freq) in predictions.items()
#                           if freq >= freq_threshold}
#         true_positives = sum([info.TP for info in predicted_hpos.values()])
#         # found_hpos = {hpo: (count, freq) for hpo, (count, freq) in predictions.items()
#         #               if hpo in test_hpos and freq >= freq_threshold}
#
#         if not test_hpos:
#             percent_found = None
#         else:
#             percent_found = true_positives / len(test_hpos)
#
#         return len(test_hpos), percent_found, len(comparison_group), predicted_hpos
# =============================================================================

    def test_phenotype_prediction(self, patient_id, freq_threshold=0,
                                  length_similarity=0, loci_similarity=0,
                                  gene_similarity=0, hpo_similarity=0):

        index_hpos = {hpo for hpo, response in self.patient_db[patient_id].hpo.items()
                      if response == "T"}

        comparison_group = self.filter_patient_comparisons(patient_id, length_similarity,
                                                           loci_similarity, gene_similarity,
                                                           hpo_similarity)
        predictions = self.predict_phenotypes2(comparison_group, show=0,
                                               additional_hpos=index_hpos,
                                               additional_groupname="index")
        return predictions

    def test_all_phenotype_predictions(self, filter_unknowns=True, freq_threshold=0,
                                       length_similarity=0, loci_similarity=0,
                                       gene_similarity=0, hpo_similarity=0):
        all_predictions = {}

        for patient_id in self.index:
            prediction = self.test_phenotype_prediction(
                patient_id, freq_threshold, length_similarity, loci_similarity,
                gene_similarity, hpo_similarity
                )
            all_predictions[patient_id] = prediction

        if filter_unknowns:
            all_predictions = {patient_id: results
                               for patient_id, results in all_predictions.items()
                               if len(results) > 0}

        return all_predictions


# =============================================================================
#     def test_all_predictions(self, filter_unknowns=True, freq_threshold=0,
#                              length_similarity=0, loci_similarity=0,
#                              gene_similarity=0, hpo_similarity=0):
#         all_predictions = {patient_id: self.test_predictions(
#             patient_id,
#             freq_threshold, length_similarity, loci_similarity,
#             gene_similarity, hpo_similarity
#             ) for patient_id in self.index}
#         if filter_unknowns:
#             all_predictions = {patient_id: results for patient_id, results in all_predictions.items()
#                                if results[0] > 0}
#         return all_predictions
#
#     def predictions_as_df(self, patient_id, predictions, threshold=0):
#         found_count, found_perc, group_size, found_hpos =  predictions[patient_id]
#         table = []
#         added = set()
#
#         for hpo, info in found_hpos.items():
#             added.add(hpo)
#             table.append([hpo.id, hpo.name, info.TP, info.count, info.frequency])
#
#         for hpo in self.patient_db[patient_id].hpo:
#             if hpo in added:
#                 continue
#             table.append([hpo.id, hpo.name, True, 0, 0])
#             # if hpo in found_hpos:
#             #     table.append([hpo.id, hpo.name, found_hpos[hpo][1] >= threshold,
#             #                   found_hpos[hpo][0], found_hpos[hpo][1]])
#             # else:
#             #     table.append([hpo.id, hpo.name, False, 0, 0])
#
#         table.sort(key=lambda line: line[-1], reverse=True)
#         table = pd.DataFrame(table, columns=["HPO_ID", "HPO_Name", "Known", "Count", "Frequency"])
#         return table
# =============================================================================


    def convert_patient_predictions_to_df(self, predictions):
        table = []
        for prediction in sorted(predictions.values(), key=lambda x: x.freq_adjusted, reverse=True):
            table.append(prediction.make_table_row())
        column_labels = [attr.title()
                         for attr in list(predictions.values())[0].__dict__.keys()
                         if not attr.startswith("_")]
        table = pd.DataFrame(table, columns=column_labels)
        return table

    # def make_summary_table(self, predictions):
    #     table = []
    #     for patient, results in predictions.items():
    #         known_count, found_perc, group_size, found_hpos = results
    #         new_hpos = sum([not hpo.TP for hpo in found_hpos.values()])
    #         # passing = sum([hpo.frequency >= threshold for hpo in found_hpos.values()
    #         #                if hpo.TP])
    #         # passing = passing / known_count
    #         # table.append([patient, known_count, found_perc, passing, group_size])
    #         table.append([patient, group_size, known_count,
    #                       found_perc, None,
    #                       new_hpos, None])
    #     table = pd.DataFrame(
    #         sorted(table, key=lambda x: x[0]),
    #         columns=["ID", "Group_Size", "Known_HPO_Terms",
    #                  "Known_Found", "Known_Above_Threshold",
    #                  "New_Found", "New_Above_Threshold"]
    #         )
    #     return table

    def make_summary_table(self, patient_predictions):
        table = []
        for patient, predictions in patient_predictions.items():
            predictions = list(predictions.values())

            population = predictions[0].population
            median_adjusted_population = np.median(
                [pred.population_adjusted for pred in predictions]
                )

            known = [pred for pred in predictions if pred.group == "index"]
            known_count = len(known)
            if known_count == 0:
                known_found = 0
            else:
                known_found = sum(1 for pred in known if pred._found) / known_count
            # known_above_thresh = (sum(1 for pred in known
            #                           if pred.found and pred.freq_adjusted >= threshold)
            #                       / known_found)

            new_count = sum(1 for pred in predictions if pred._found and pred.group == "predicted")

            table.append([patient, population, known_count, known_found, None,
                          new_count, None])
        table = pd.DataFrame(
            sorted(table, key=lambda x: x[0]),
            columns=["ID", "Group_Size", "Known_HPO_Terms",
                     "Known_Found", "Known_Above_Threshold",
                     "New_Found", "New_Above_Threshold"]
            )
        return table

    @staticmethod
    def write_excel_sheet(table, name, writer=None, path=None, threshold=0,
                          formats=None, widths=None, conditionals=None):
        close = False
        if writer is None:
            writer = pd.ExcelWriter(path, engine="xlsxwriter")
            close = True

        table.to_excel(writer, sheet_name=name, startrow=3, index=False)
        col_names = [{"header": col_name} for col_name in table.columns]
        if name == "Summary":
            formula = ('=COUNTIFS('
                       'INDIRECT("Table_" & [@ID] & "[Freq_Adjusted]"), ">=" & $B$2,'
                       'INDIRECT("Table_" & [@ID] & "[Group]"), "=index")'
                       '/[@[Known_HPO_Terms]]')
            formula2 = ('=COUNTIFS('
                       'INDIRECT("Table_" & [@ID] & "[Freq_Adjusted]"), ">=" & $B$2,'
                       'INDIRECT("Table_" & [@ID] & "[Group]"), "=predicted")')
            col_names[4]["formula"] = formula
            col_names[6]["formula"] = formula2

        workbook = writer.book
        perc_format = workbook.add_format({"num_format": 10})

        worksheet = writer.sheets[name]
        worksheet.write("A1", name)
        worksheet.write("A2", "Threshold:")
        worksheet.write("B2", threshold, perc_format)

        table_name = "Table_" + name.replace("-", "_")

        worksheet.add_table(3, 0, table.shape[0]+3, table.shape[1]-1,
                            {"columns": col_names,
                             "name": table_name,
                             "style": "Table Style Light 1"})

        empty = [None] * len(col_names)
        if formats is None:
            formats = empty
        if widths is None:
            widths = empty
        if conditionals is None:
            conditionals = []

        for i, (f, w) in enumerate(zip(formats, widths)):
            worksheet.set_column(i, i, w, f)

        for conditional in conditionals:
            worksheet.conditional_format(*conditional)

        if close:
            writer.close()

    def write_all_predictions_to_excel(self, predictions, path, threshold=0):
        writer = pd.ExcelWriter(path, engine="xlsxwriter")

        try:
            workbook = writer.book

            perc_format = workbook.add_format({"num_format": 10})
            green = workbook.add_format({'bg_color': '#C6EFCE',
                                         'font_color': '#006100'})
            red = workbook.add_format({"bg_color": "#E6B8B7",
                                       "font_color": "#C0504D"})
            center_url = workbook.add_format({"center_across": True,
                                              "font_color": "#0000FF",
                                              "underline": True})

            # Add summary sheet at the beginning of the workbook.
            table = self.make_summary_table(predictions)
            self.write_excel_sheet(table=table, name="Summary", writer=writer, threshold=threshold,
                                   formats=[None, None, None, perc_format, perc_format, None, None],
                                   widths=[20, 12, 20, 20, 25, 20, 25])

            # Make patient ID's clickable on summary page.
            worksheet = writer.sheets["Summary"]
            for i, patient in enumerate(table["ID"], start=4):
                worksheet.write_url(i, 0, f"internal:{patient}!A1", string=patient)

            # Write individual patient sheets.
            for patient_id in sorted(predictions):
                table = self.convert_patient_predictions_to_df(predictions[patient_id])

                # Setup conditional formatting.
                freq_col = table.columns.get_loc("Freq_Adjusted")
                group_col = table.columns.get_loc("Group")
                conditionals = [(4, freq_col, table.shape[0]+3, freq_col,
                                 {'type': 'cell',
                                  'criteria': ">=",
                                  'value': "'Summary'!$B$2",
                                  'format': green}),
                                (4, group_col, table.shape[0]+3, group_col,
                                 {'type': 'cell',
                                  'criteria': "==",
                                  'value': '"predicted"',
                                  'format': red})]

                self.write_excel_sheet(table=table, name=patient_id,
                                       writer=writer, threshold="='Summary'!$B$2",
                                       formats=[None]*8 + [perc_format, None, perc_format],
                                       widths=[12, 50] + [12]*9,
                                       conditionals=conditionals)
                worksheet = writer.sheets[patient_id]
                worksheet.write_url(0, 4, "internal:Summary!A1",
                                    string="Home", cell_format=center_url)

        except Exception as error:
            writer.close()
            raise error
        writer.close()


class PatientIntersect:
    """Record of similarity between two patients."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, patient_1, patient_2, cnv_type,
                 length_compare, loci_compare,
                 gene_compare, hi_gene_compare,
                 hpo_compare):
        self.patients = {patient_1.id: patient_1, patient_2.id: patient_2}
        self.ids = {patient_1.id, patient_2.id}
        self.change = cnv_type

        self.length_similarity = length_compare

        self.loci_similarity = loci_compare[0]
        self.loci_shared_size = loci_compare[1]

        self.gene_similarity = gene_compare[0]
        self.genes = gene_compare[1]
        self.gene_count = len(gene_compare[1])

        self.hi_gene_similarity = hi_gene_compare[0]
        self.hi_genes = hi_gene_compare[1]
        self.hi_gene_count = len(hi_gene_compare[1])

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

# XXX: This probably doesn't work anymore.
# def write_comparison_table(table, patients, out, self_match="size"):
#     """Write comparison table to file. Deprecated."""
#     header = "PatientID"
#     writer = []
#     for pid in table:
#         header += f",{pid}"
#         string = pid
#         for p2 in table:
#             if pid == p2:
#                 if self_match == "size":
#                     string += f",{len(patients[pid].all_genes())}"
#                 else:
#                     string += f",{self_match}"
#             elif p2 in table[pid]:
#                 string += f",{table[pid][p2]}"
#             else:
#                 string += f",{table[p2][pid]}"
#         string += "\n"
#         writer.append(string)
#     writer = [header + "\n"] + writer
#     with open(out, "w") as outfile:
#         outfile.writelines(writer)


# %% Predictions
class TraitPrediction:
    def __init__(self, trait_id, trait_name, population, true, false,
                 unsure, na, group=None):
        self.trait_id = trait_id
        self.trait_name = trait_name
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
        string = (f"{self.trait_name}\t{self.population}\t{self.population_adjusted}\t"
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
        row = []
        for attr, value in self.__dict__.items():
            if attr.startswith("_"):
                continue
            if isinstance(value, Term):
                row.append(value.id)
                continue
            row.append(value)
        return row


class PredictionDatabase:
    def __init__(self, predictions):
        self.predictions = predictions

    def __getitem__(self, key):
        """Retrieve prediction using patient ID."""
        return self.predictions[key]

    def get_sizes(self):
        sizes = {patient: list(predictions.values())[0].population
                 for patient, predictions in self.predictions.items()}
        return sizes


PredictInfo = namedtuple("PredictionInfo",
                         ["count", "frequency", "TP"])


# %% Main
def main(geno_file, pheno_file):
    """Run main."""
    genotypes = DataManager.read_data(geno_file)
    phenotypes = DataManager.read_data(pheno_file)
    patients = DataManager.make_patients(genotypes, phenotypes)
    DataManager.print_summary_counts(patients)


def test(genotypes, phenotypes, patient_hpo, drop_list_file, expand_hpos=False):
    """Test run."""
    # Read genome dictionary.
    print("Reading reference genome dictionary...")
    genomedict = GenomeDict("Data/human_g1k_v37_phiX.dict")

    # Read geneset.
    print("Loading gene set...")
    # Build GeneSet from source GTF file (gzipped) (slow, large file):
    # geneset = GeneSet("C:/Users/Ty/Documents/Chr6/hg19.ensGene.gtf.gz")

    # Build chr6 GeneSet only from source GTF file (gzipped) (slow, large file)
    # and add pLI and HI info automatically from default sources:
    geneset = gene_set.main()

    # Or, load pre-made GeneSet from pickle (faster, large file):
    # with open("GeneSets/hg19.ensGene.pkl", "rb") as infile:
    #     geneset = pickle.load(infile)

    # Or, load pre-made GeneSet from bz2 pickle (less fast, small file):
    # with bz2.BZ2File("GeneSets/hg19.ensGene.pkl.bz2", "rb") as infile:
    #     geneset = cPickle.load(infile)

    # Read HPO ontology.
    print("Loading Human Phenotype Ontology...")
    ontology = Ontology("Data/hpo.obo")

    # Read patient genotypes.
    print("Reading patient genotype data...")
    genotypes = DataManager.read_data(genotypes)
    genotypes = DataManager.fix_genotype_data(genotypes)
    # genotypes = trim_chromosome_names(genotypes)

    # Read patient phenotypes.
    print("Reading patient phenotype data...")
    phenotypes = DataManager.read_data(phenotypes)

    # Read patient HPO terms.
    print("Reading patient HPO data...")
    hpos = DataManager.read_data(patient_hpo)
    hpos = DataManager.fix_patient_hpos2(hpos)

    # Build patient objects.
    # !!!: This is where you can choose whether or not to expand HPO terms.
    print("Building patient objects...")
    patients = DataManager.make_patients(genotypes, phenotypes, geneset,
                                         hpos, ontology, expand_hpos=expand_hpos)

    print("Filtering patients...")
    drop_list = DataManager.read_patient_drop_list(drop_list_file)
    patients = DataManager.filter_patients(patients, drop_list)

    patients = PatientDatabase(patients)

    print("Running comparisons...")
    comparison = ComparisonTable(patients)

    print("Done.")
    return (
        genotypes,
        phenotypes,
        patients,
        comparison,
        geneset,
        ontology,
        # hi_genes,
        genomedict
        )


def predict_test(comparison, patient_list="C:/Users/Ty/Documents/Chr6/Predict_tests/test_patients.txt"):
    with open(patient_list) as infile:
        aafkes_patients = infile.readlines()
    aafkes_patients = [x.strip() for x in aafkes_patients]
    all_tests = comparison.test_all_phenotype_predictions(gene_similarity=.7)
    aafkes_tests = {x: y for x, y in all_tests.items() if x in aafkes_patients}
    return all_tests, aafkes_tests


if __name__ == "__main__":
    # import sys
    # main(sys.argv[1], sys.argv[2])
    # genotypes, phenotypes, patients, geneset = test()
    _, _, my_patients, my_comparison, my_geneset, my_ontology, my_genomedict = test(
        genotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_array_2022-02-21_18_28_06.csv",
        phenotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_questionnaire_2022-02-21_10_18_09.csv",
        patient_hpo="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_research_patients_2022-02-21_10_20_53.csv",
        drop_list_file="/home/tyler/Documents/Chr6_docs/PatientData/drop_list.txt",
        expand_hpos=False
        )
    table, selected_hpos, homos = phenotype_homo_test(my_comparison, my_ontology, .85, .2, 2)
