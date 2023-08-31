#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Patient Comparisons.

@author: T.D. Medina
"""

import numpy as np
import pandas as pd

from chr6_project.analysis.data_objects import PatientDatabase
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

    def __init__(self, patient_db=None, pLI_threshold=0.9, HI_threshold=10,
                 phaplo_threshold=0.86, comparison_table=None, mode="confirm"):
        if comparison_table is not None:
            self.read_from_existing(comparison_table)
            return
        self.patient_db = patient_db
        # self.raw = self.compare_patients(pLI_threshold, HI_threshold,
        #                                  phaplo_threshold, mode)
        # self.index = self.make_index()
        # self.array = self.make_array()
        self.index, self.array = self.compare_patients(pLI_threshold, HI_threshold,
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

    # def make_index(self):
    #     """Create name-to-number mapping to look up names in numpy array."""
    #     return {j: i for i, j in enumerate(self.raw)}

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
        # self.raw = comparison_table.raw
        # self.index = comparison_table.index
        self.array = comparison_table.array
        self.size = comparison_table.size

    def compare_patients(self, pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=0.86,
                         mode="any"):
        """Compare all patients to each other."""
        ids = list(self.patient_db.patients.keys())
        comparisons = {}
        while ids:
            id_i = ids.pop()
            patient_i = self.patient_db[id_i]
            patient_comparison = dict()
            patient_comparison[id_i] = self.compare_all(patient_i, patient_i,
                                                        pLI_threshold, HI_threshold,
                                                        phaplo_threshold, mode)

            for id_j in ids:
                patient_j = self.patient_db[id_j]
                patient_comparison[id_j] = self.compare_all(patient_i, patient_j,
                                                            pLI_threshold, HI_threshold,
                                                            phaplo_threshold, mode
                                                            )
            comparisons[id_i] = patient_comparison
        index, comparisons = self._make_comparison_array(comparisons)
        return index, comparisons

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

    # def recompare_patients(self, pLI_threshold=0.9, HI_threshold=10,
    #                        phaplo_threshold=0.86, mode="any"):
    #     self.index, self.array = self.compare_patients(pLI_threshold, HI_threshold,
    #                                                    phaplo_threshold, mode)
    #     # self.raw = self.compare_patients(pLI_threshold, HI_threshold, phaplo_threshold, mode)
    #     # self.array = self.make_array()

    @classmethod
    def compare_all(cls, patient_1, patient_2,
                    pLI_threshold=0.9, HI_threshold=10,
                    phaplo_threshold=0.86, mode="any"):
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
        hi_compare = cls.compare_HI_genes(patient_1, patient_2,
                                          pLI_threshold, HI_threshold,
                                          phaplo_threshold, mode)
        dom_compare = cls.compare_dominant_genes(patient_1, patient_2)
        dom_match = cls.check_dominant_match(patient_1, patient_2)
        hpo_compare = cls.compare_hpos(patient_1, patient_2)
        comparison = PatientIntersect(
            patient_1, patient_2, cnv_type,
            length_compare, loci_compare,
            gene_compare, hi_compare,
            dom_compare, dom_match,
            hpo_compare
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
                         pLI_threshold=0.9, HI_threshold=10,
                         phaplo_threshold=0.86, mode="any"):
        """Compare affected HI genes between two patients."""
        jaccard_index, intersect = jaccard(
            patient_1.all_HI_genes(pLI_threshold, HI_threshold, phaplo_threshold, mode),
            patient_2.all_HI_genes(pLI_threshold, HI_threshold, phaplo_threshold, mode)
            )
        return jaccard_index, intersect

    @staticmethod
    def compare_dominant_genes(patient_1, patient_2):
        """Compare affected dominant-effect genes between two patients."""
        jaccard_index, intersect = jaccard(patient_1.all_dominant_genes(),
                                           patient_2.all_dominant_genes())
        return jaccard_index, intersect

    @staticmethod
    def check_dominant_match(patient_1, patient_2):
        return patient_1.all_dominant_genes() == patient_2.all_dominant_genes()

    @staticmethod
    def compare_hpos(patient_1, patient_2):
        """Compare HPO terms between two patients."""
        hpo_set1 = {hpo for hpo, response in patient_1.hpo.items() if response == "T"}
        hpo_set2 = {hpo for hpo, response in patient_2.hpo.items() if response == "T"}
        jaccard_index, intersect = jaccard(hpo_set1, hpo_set2)
        return jaccard_index, intersect

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

    def convert_to_3d_array(self):
        array = np.ndarray([6, *self.array.shape])
        for intersect in self:
            if not intersect.self_compare:
                pid1, pid2 = intersect.ids
            else:
                pid1 = list(intersect.ids)[0]
                pid2 = pid1
            pid1, pid2 = self.index[pid1], self.index[pid2]
            for i, sim in enumerate(intersect.get_similarities()):
                array[i, pid1, pid2] = sim
                array[i, pid2, pid1] = sim
        return array
    # def export_tables_per_comparator(self):
    #     tables = dict()
    #     comparators = ["length_similarity", "loci_similarity", "gene_similarity",
    #                    "hi_gene_similarity", "dom_gene_match"]
    #     for comparator in comparators:
    #         table = {pid1: {pid2: intersect[]}}

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

    def filter_patient_intersects(self, patient_id,
                                  length_similarity=0, loci_similarity=0,
                                  gene_similarity=0, hi_gene_similarity=0,
                                  dom_gene_match=True, hpo_similarity=0,
                                  include_self=False, as_patient_database=False):
        intersects = [inter for inter in self.lookup(patient_id, "all") if all(
            [inter.length_similarity >= length_similarity,
             inter.loci_similarity >= loci_similarity,
             inter.gene_similarity >= gene_similarity,
             inter.hi_gene_similarity >= hi_gene_similarity,
             inter.dom_gene_match or not dom_gene_match,
             inter.hpo_similarity >= hpo_similarity,
             not inter.self_compare or include_self]
            )]
        if as_patient_database:
            patients = [inter.get_other_patient(patient_id) for inter in intersects]
            patients = {patient.id: patient for patient in patients}
            patients = PatientDatabase(patients)
            return patients
        return intersects

    # def filter_patient_comparisons(self, patient_id,
    #                                length_similarity=0, loci_similarity=0,
    #                                gene_similarity=0, hi_similarity=0,
    #                                hpo_similarity=0, split_by_dom=False):
    #     intersections = self.lookup(patient_id, "all")
    #     filtered = []
    #     for intersect in intersections:
    #         if intersect.self_compare:
    #             continue
    #         # if intersect.patients[0] == intersect.patients[1]:
    #         #     continue
    #         patient2 = intersect.get_other_patient(patient_id)
    #         # if patient1 != intersect.patients[0]:
    #         #     patient2 = intersect.patients[0]
    #         # else:
    #         #     patient2 = intersect.patients[1]
    #
    #         if not patient2.hpo:
    #             continue
    #
    #         if not all([intersect.length_similarity >= length_similarity,
    #                     intersect.loci_similarity >= loci_similarity,
    #                     intersect.gene_similarity >= gene_similarity,
    #                     intersect.hi_gene_similarity >= hi_similarity,
    #                     (not split_by_dom) or intersect.dom_gene_similarity in {0, 1},
    #                     intersect.hpo_similarity >= hpo_similarity]):
    #             continue
    #
    #         filtered.append(patient2)
    #     return filtered

    # def filter_patient_intersections(self, patient_id, length_similarity=0,
    #                                  loci_similarity=0, gene_similarity=0,
    #                                  hi_similarity=0, hpo_similarity=0,
    #                                  include_self=False):
    #     intersections = {intersect.get_other_id(patient_id): intersect
    #                      for intersect in self.lookup(patient_id, "all")
    #                      if all([intersect.length_similarity >= length_similarity,
    #                              intersect.loci_similarity >= loci_similarity,
    #                              intersect.gene_similarity >= gene_similarity,
    #                              intersect.hi_gene_similarity >= hi_similarity,
    #                              intersect.hpo_similarity >= hpo_similarity])}
    #     if include_self is False and patient_id in intersections:
    #         del intersections[patient_id]
    #     return intersections

    # def make_patient_intersection_group(self, patient_id, length_similarity=0,
    #                                     loci_similarity=0, gene_similarity=0,
    #                                     hi_similarity=0, hpo_similarity=0,
    #                                     include_self=True, as_patient_database=True):
    #     intersections = self.filter_patient_intersections(
    #         patient_id, length_similarity, loci_similarity, gene_similarity,
    #         hi_similarity, hpo_similarity, include_self
    #         )
    #     patients = [intersection.get_other_patient(patient_id)
    #                 for intersection in intersections.values()]
    #     if as_patient_database:
    #         patients = PatientDatabase({patient.id: patient for patient in patients})
    #     return patients

    def compare_patient_pheno_prevalence(self, patient_id, phenotypes,
                                         length_similarity=0, loci_similarity=0,
                                         gene_similarity=0, hi_gene_similarity=0,
                                         dom_gene_match=True, hpo_similarity=0):
        group = self.filter_patient_intersects(
            patient_id, length_similarity, loci_similarity,
            gene_similarity, hi_gene_similarity, dom_gene_match, hpo_similarity,
            include_self=True, as_patient_database=True
            )
        size = group.size
        prevs = {
            phenotype: sum([patient.hpo[phenotype] == "T" for patient in group
                            if phenotype in patient.hpo])
            for phenotype in phenotypes
            }
        prevs = {phenotype: PhenotypePrevalence(patient_id, group, phenotype, size, prevalence)
                 for phenotype, prevalence in prevs.items()}
        prevs = PatientGroupPrevalences(patient_id, group, prevs)
        return prevs

    def compare_all_patient_pheno_prevalences(self, phenotypes, length_similarity=0,
                                              loci_similarity=0, gene_similarity=0,
                                              hi_gene_similarity=0, dom_gene_match=True,
                                              hpo_similarity=0):
        params = [phenotypes, length_similarity, loci_similarity, gene_similarity,
                  hi_gene_similarity, dom_gene_match, hpo_similarity]
        all_prevs = {patient: self.compare_patient_pheno_prevalence(patient, *params)
                     for patient in self.index.keys()}
        all_prevs = HomogeneityDatabase(all_prevs)
        return all_prevs

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
    def predict_phenotypes(comparison_group, show=10,
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

        all_hpo = {hpo: TraitPrediction(hpo, population, *counts.values())
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

    def test_phenotype_prediction(self, patient_id,
                                  length_similarity=0, loci_similarity=0,
                                  gene_similarity=0, hi_gene_similarity=0,
                                  dom_gene_match=True, hpo_similarity=0,
                                  skip_no_hpos=True):
        params = [patient_id, length_similarity, loci_similarity, gene_similarity,
                  hi_gene_similarity, dom_gene_match, hpo_similarity]

        index_hpos = {hpo for hpo, response in self.patient_db[patient_id].hpo.items()
                      if response == "T"}

        comparison_group = self.filter_patient_intersects(*params, include_self=False,
                                                          as_patient_database=True)
        comparison_group = list(comparison_group.patients.values())

        if skip_no_hpos:
            comparison_group = [patient for patient in comparison_group if patient.hpo]

        predictions = self.predict_phenotypes(comparison_group, show=0,
                                              additional_hpos=index_hpos,
                                              additional_groupname="index")
        predictions = PatientPredictions(
            patient=patient_id,
            patient_group=PatientDatabase({patient.id: patient for patient in comparison_group}),
            predictions=predictions
            )
        return predictions

    def test_all_phenotype_predictions(self, length_similarity=0, loci_similarity=0,
                                       gene_similarity=0, hi_gene_similarity=0,
                                       dom_gene_match=True, hpo_similarity=0,
                                       skip_no_hpos=True, filter_unknowns=True):
        params = [length_similarity, loci_similarity, gene_similarity, hi_gene_similarity,
                  dom_gene_match, hpo_similarity, skip_no_hpos]

        all_predictions = {patient_id: self.test_phenotype_prediction(patient_id, *params)
                           for patient_id in self.index}

        if filter_unknowns:
            all_predictions = {patient_id: results
                               for patient_id, results in all_predictions.items()
                               if len(results) > 0}

        all_predictions = PredictionDatabase(all_predictions)
        return all_predictions

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
                 dom_gene_compare, dom_gene_match,
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
