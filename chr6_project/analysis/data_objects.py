#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Data objects.

@author: T.D. Medina
"""

from collections import namedtuple
import csv
from datetime import datetime

from numpy import mean, median

from chr6_project.analysis.utilities import (
    merge_range_list,
    overlap,
    REFERENCE_CHR,
    )
from chr6_project.analysis.gene_set import is_haploinsufficient


SequenceContig = namedtuple("SequenceContig",
                            ["name", "length", "cumulative_start"])


class GenomeDict:
    """GATK-style genome dictionary containing contig names and lengths."""

    def __init__(self, file=None):
        self.sequences = {}
        self.index = {}
        self.total = 0

        if file is not None:
            self.sequences = self.read_genome_dict(file)
            self.index = {name: i for i, name in enumerate(self.sequences)}
            self.total = sum([seq.length for seq in self.sequences.values()])

    def __len__(self):
        return len(self.index)

    @staticmethod
    def read_genome_dict(file):
        """Read in GATK genome dictionary file."""
        with open(file) as infile:
            lines = infile.readlines()
        lines = [line for line in lines if line.startswith("@SQ")]
        for i, line in enumerate(lines):
            line = line.strip().split("\t")
            line = [line[1].replace("SN:", "", 1), int(line[2].replace("LN:", "", 1)), 0]
            if i != 0:
                line[2] = lines[i - 1][2] + line[1]
            lines[i] = SequenceContig(*line)
        lines = {seq.name: seq for seq in lines}
        return lines

    def abs_pos(self, chromosome, position):
        """Calculate absolute position from chromosome position.

        Uses genome contig order to establish relative positon within
        cumulative chromosome lengths.
        """
        return self[chromosome].cumulative_start + position

    def __getitem__(self, key):
        """Retrieve item using key."""
        return self.sequences[key]


class Cytoband:
    """Cytogenetic banding object."""

    def __init__(self, chromosome, start, stop, band, stain):
        self.chr = chromosome
        self.locus = range(start, stop)
        self.band = band
        self.stain = stain

    def __str__(self):
        """Make custom pretty string representation."""
        string = ("Cytoband("
                  f"locus='{self.chr}:{self.locus.start}-{self.locus.stop - 1}', "
                  f"band='{self.band}', "
                  f"stain='{self.stain}')")
        return string

    def __repr__(self):
        """Make custom technical string representation."""
        string = f"Cytoband({self.chr}:{min(self.locus)}-{max(self.locus)})"
        return string


class Cytomap:
    """Cytogenetic banding map of a genome."""

    def __init__(self, file="C:/Users/Ty/Documents/cytoBand.txt"):
        self.path = file
        self.cytobands = self.make_cytobands(self.path)

    @classmethod
    def make_cytobands(cls, filepath):
        """Set up cytogenetic map from data file."""
        cytobands = cls.read_cytoband_file(filepath)
        cytobands = cls.split_cytobands_by_chr(cytobands)
        return cytobands

    @staticmethod
    def read_cytoband_file(file):
        """Read cytogenetic banding from file."""
        with open(file) as infile:
            cytobands = infile.readlines()
        cytobands = [line.strip().split("\t") for line in cytobands]
        cytobands = [Cytoband(line[0].strip("chr"), int(line[1]), int(line[2]), line[3], line[4])
                     for line in cytobands]
        return cytobands

    @staticmethod
    def split_cytobands_by_chr(cytobands):
        """Organize cytogenetic bands by chromosome."""
        cytomap = {chrom: [] for chrom in {cytoband.chr for cytoband in cytobands}}
        for cytoband in cytobands:
            cytomap[cytoband.chr].append(cytoband)
        return cytomap

    @staticmethod
    def split_cytomap_by_arm(cytomap):
        """Split cytogenetic bands by chromosome arm."""
        arm_map = {key: {"p": [], "q": []} for key in cytomap.keys()}
        for chrom_list in cytomap.values():
            for cytoband in chrom_list:
                arm_map[cytoband.chr][cytoband.band[0]].append(cytoband)
        return arm_map

    def get_band(self, chromosome, coordinate):
        """Get corresponding cytogenetic band at a genomic locus."""
        for band in self.cytobands[chromosome]:
            if coordinate in band.locus:
                return band
        return None

    def get_bands(self, chromosome, range_: range):
        """Get corresponding cytogenetic bands in a genomic locus range."""
        bands = []
        for band in self.cytobands[chromosome]:
            if overlap(range_, band.locus):
                bands.append(band)
        return bands


class CNV:
    """Data object containing CNV information."""

    def __init__(self, chromosome, start, stop, change, ID=None, genes=None):
        self.id = ID
        self.chromosome = str(chromosome)
        self.range = range(start, stop + 1)
        self.change = change
        self.length = len(self.range)
        self.genes = genes

    def __len__(self):
        return self.length

    def __repr__(self):
        """Make string representation of object."""
        string = (
            f"CNV({self.change}: "
            f"{self.chromosome}:{self.range.start}-{self.range.stop})"
            )
        return string

    def __str__(self):
        """Make pretty string of object."""
        string = ("CNV:\n"
                  f"  Change = {self.change}\n"
                  f"  Locus = {self.chromosome}:"
                  f"{self.range.start}-{self.range.stop}\n"
                  f"  Length = {self.length:,} bp")
        return string


# TODO: Add a method to add a patient.
class PatientDatabase:
    """Database containing all Patient objects."""

    def __init__(self, patients):
        self.patients = patients
        self.index = self.make_index()
        self.cnvs = self.organize_cnvs()
        self.size = len(self.patients)

    def __len__(self):
        return self.size

    def __getitem__(self, key):
        """Get patient by patient ID."""
        return self.patients[key]

    def __iter__(self):
        """Initialize iterable."""
        self.__iteri__ = 0
        return self

    def __next__(self):
        """Iterate iterable."""
        if self.__iteri__ == self.size:
            raise StopIteration
        result = self.patients[self.index[self.__iteri__]]
        self.__iteri__ += 1
        return result

    def make_index(self):
        """Make database index for iteration purposes."""
        index = dict(enumerate(self.patients))
        return index

    def organize_cnvs(self):
        """Get and organize CNVs from all patients."""
        cnvs = sorted(
            [cnv for patient in self.patients.values() for cnv in patient.cnvs],
            key=lambda x: x.range.start
            )
        cnv_dict = {chromosome: [] for chromosome in {cnv.chromosome for cnv in cnvs}}
        for cnv in cnvs:
            cnv_dict[cnv.chromosome].append(cnv)
        return cnv_dict

    def list_ids(self):
        id_list = list(self.patients.keys())
        return id_list

    def get_median_cnv_position(self, chromosome):
        if chromosome not in self.cnvs:
            return 0
        cnv_median = median([(cnv.range.start + cnv.range.stop)/2 for cnv in self.cnvs[chromosome]])
        cnv_median = int(cnv_median)
        return cnv_median

    def get_mean_cnv_position(self, chromosome):
        cnv_mean = mean([(cnv.range.start + cnv.range.stop)/2 for cnv in self.cnvs[chromosome]])
        return cnv_mean

    def add_predictions(self, predictions):
        for patient in self:
            patient.predictions = predictions[patient.id]

    def summary(self):
        """Calculate summary counts."""
        len_cnv = len([cnv for chromosome in self.cnvs.values() for cnv in chromosome])
        summary_txt = (f"Patients: {len(self.patients)}\n"
                       f"CNVs: {len_cnv}")
        return summary_txt

    def make_phenotype_table(self):
        """Tabulate phenotype info for all patients."""
        all_hpos = sorted(list({hpo.id for patient in self for hpo in patient.hpo}))
        all_hpos = {name: pos for pos, name in enumerate(all_hpos)}
        hpo_count = len(all_hpos)

        header = ["patient"] + list(all_hpos.keys())
        rows = {}

        for patient in self:
            row = [0] * hpo_count
            for hpo in patient.hpo:
                row[all_hpos[hpo.id]] = 1
            row.insert(0, patient.id)
            row = dict(zip(header, row))
            rows[row["patient"]] = row
        return rows

    def make_gene_table(self):
        """Tabulate gene info for all patients."""
        all_genes = sorted(list({gene.gene_id for patient in self for gene in patient.all_genes()}))
        all_genes = {name: pos for pos, name in enumerate(all_genes)}
        gene_count = len(all_genes)

        header = ["patient"] + list(all_genes.keys())
        rows = {}

        for patient in self:
            row = [0] * gene_count
            for gene in patient.all_genes():
                row[all_genes[gene.gene_id]] = 1
            row.insert(0, patient.id)

            row = dict(zip(header, row))
            rows[row["patient"]] = row
        return rows

    def make_genotypes_table(self, absolute_pos=True, genome_dict=None):
        """Tabulate CNV information for all single-CNV patients."""
        if absolute_pos and genome_dict is None:
            raise ValueError("Asbolute positions require a reference genome dictionary.")
        header = ["patient", "start", "length"]

        rows = {}
        for patient in self:
            if len(patient.cnvs) != 1:
                continue
            cnv = patient.cnvs[0]
            if absolute_pos:
                row = dict(zip(header, [patient.id, genome_dict.abs_pos(cnv.chromosome, cnv.range.start), cnv.length]))
            else:
                row = dict(zip(header, [patient.id, cnv.range.start, cnv.length]))
            rows[row["patient"]] = row
        return rows

    def make_combination_table(self, genotypes=True, genes=True,
                               phenotypes=False, absolute_pos=True,
                               genome_dict=None):
        """Tabulate multiple types of patient data."""
        subtables = []
        if genotypes:
            subtables.append(self.make_genotypes_table(absolute_pos,
                                                       genome_dict))
        if genes:
            subtables.append(self.make_gene_table())
        if phenotypes:
            subtables.append(self.make_phenotype_table())

        names = set(subtables[0].keys())
        names.intersection_update(*[
            set(subtable.keys()) for subtable in subtables
            ])

        rows = {}
        for name in names:
            row = {}
            for subtable in subtables:
                row.update(subtable[name])
            rows[name] = row
        return rows

    @staticmethod
    def write_table(table, out):
        """Write CSV file of tabulated data."""
        table_values = list(table.values())
        with open(out, "w", newline="") as outfile:
            writer = csv.DictWriter(outfile, table_values[0].keys())
            writer.writeheader()
            writer.writerows(table_values)


class Patient:
    """Chromosome 6 Project patient data."""

    def __init__(self, ID):
        self.id = ID
        self.genotypes = []
        self.phenotypes = []
        self.cnvs = []
        # self.affected_genes = {}
        # self.affected_gene_ids = {}
        self.hpo = {}
        self.predictions = {}

        self.origin = None
        self.birthdate = None
        self.age = None

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id

    def __repr__(self):
        """Construct string representation."""
        return (f"ID: {self.id}\n"
                f"  Genotypes: {len(self.genotypes)}\n"
                f"  Phenotypes: {len(self.phenotypes)}")

    def extract_cnvs(self):
        """Pull CNV data from genotype information."""
        cnvs = []
        for genotype in self.genotypes:
            chrom = genotype["Chromosome"]
            start = genotype["Start positie in Hg19"]
            stop = genotype["Stop positie in Hg19"]
            change = genotype["imbalance"]

            if chrom not in REFERENCE_CHR or not start or not stop:
                continue
            cnvs.append(CNV(chrom, start, stop, change, self.id))
        cnvs = sorted(cnvs, key=lambda x: REFERENCE_CHR.index(x.chromosome))
        return cnvs

    def get_median_cnv_position(self, chromosome):
        cnvs = [cnv for cnv in self.cnvs if cnv.chromosome == chromosome]
        if not cnvs:
            return None
        cnv_median = median([(cnv.range.start + cnv.range.stop)/2 for cnv in cnvs])
        cnv_median = int(cnv_median)
        return cnv_median

    def get_mean_cnv_position(self, chromosome):
        cnvs = [cnv for cnv in self.cnvs if cnv.chromosome == chromosome]
        if not cnvs:
            return None
        cnv_mean = mean([(cnv.range.start + cnv.range.stop)/2 for cnv in cnvs])
        return cnv_mean

    def get_affected_ranges(self):
        """Get CNV ranges from all CNVs."""
        ranges = {cnv.chromosome: [] for cnv in self.cnvs}
        for cnv in self.cnvs:
            ranges[cnv.chromosome].append(cnv.range)
        ranges = {chromosome: merge_range_list(cnv_ranges)
                  for chromosome, cnv_ranges in ranges.items()}
        return ranges

    def identify_gene_overlaps(self, gene_set_obj):
        """Set genes affected per CNV."""
        for cnv in self.cnvs:
            cnv.genes = gene_set_obj.get_locus(cnv.chromosome, cnv.range.start,
                                               cnv.range.stop)

    def all_genes(self):
        """Get all genes affected by all CNVs."""
        all_genes = {gene for cnv in self.cnvs for gene in cnv.genes}
        return all_genes

    def all_HI_genes(self, pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=.86):
        """Get all haploinsufficient genes affected by all CNVs."""
        hi_genes = {gene for cnv in self.cnvs for gene in cnv.genes
                    if is_haploinsufficient(gene, pLI_threshold, HI_threshold, phaplo_threshold)}
        return hi_genes

    def all_dominant_genes(self):
        dom_genes = {gene for cnv in self.cnvs for gene in cnv.genes
                     if gene.dominant}
        return dom_genes

    def get_true_hpos(self):
        trues = {term for term, response in self.hpo.items()
                 if response == "T"}
        return trues

    def expand_hpo_terms(self):
        expanded = {parent for term, response in self.hpo.items()
                    for parent in term.superclasses()
                    if response == "T"}
        for hpo in expanded:
            self.hpo[hpo] = "T"

    def convert_birthday_to_datetime(self):
        if not self.phenotypes:
            return
        if not self.phenotypes[0]["birthdate"]:
            return
        self.birthdate = datetime.strptime(
            self.phenotypes[0]["birthdate"],
            "%Y-%m-%d")
        years = datetime.today().year - self.birthdate.year
        if (datetime.today().month <= self.birthdate.year
                and datetime.today().day < self.birthdate.day):
            years -= 1
        self.age = years
        return
