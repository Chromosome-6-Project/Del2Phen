"""The Chromosome 6 Project - Data objects.

@author: T.D. Medina
"""

from collections import defaultdict, namedtuple
from textwrap import dedent
from typing import Dict, List, Optional, Set, Union

import pandas as pd
from numpy import mean, median

from chr6_project.analysis.gene_set import Gene, GeneSet, is_haploinsufficient
from chr6_project.analysis.utilities import (
    merge_range_list,
    overlap,
    )

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

        Uses genome contig order to establish relative position within
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

    def __init__(self, chromosome, start, stop, copy_number, id=None, genes=None, **kwargs):
        self.id = id
        self.chromosome = str(chromosome)
        self.range = range(start, stop + 1)
        self.copy_number = copy_number
        self.genes = genes

    def __len__(self):
        return self.length

    def __repr__(self):
        string = (
            "CNV("
            f"chromosome='{self.chromosome}', "
            f"start={self.start}, "
            f"stop={self.stop}, "
            f"copy_number='{self.copy_number}', "
            f"id='{self.id}'"
            ")"
            )
        return string

    def __str__(self):
        string = ("CNV:\n"
                  f"  Copy number = {self.copy_number}\n"
                  f"  Locus = {self.chromosome}:"
                  f"{self.range.start}-{self.range.stop}\n"
                  f"  Length = {self.length:,} bp")
        return string

    @property
    def start(self):
        return self.range.start

    @property
    def stop(self):
        return self.range.stop

    @property
    def length(self):
        return len(self.range)

    def set_genes(self, geneset: GeneSet):
        """Assign overlapping genes to CNV based on provided GeneSet object."""
        self.genes = geneset.get_locus(self.chromosome, self.start, self.stop)


class Patient:
    """Patient CNV data object."""

    def __init__(self, patient_id):
        self.id = patient_id
        self.phenotypes = dict()
        self.cnvs = []

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id

    def __repr__(self):
        return f"Patient(id={self.id})"

    def filter_cnvs(self, chromosomes: Optional[Union[str, List[str]]] = None,
                    copy_numbers: Optional[Union[int, List[int]]] = None) -> List[CNV]:
        """Return a list of CNVs filtered by affected chromosome and copy number.

        Desired chromosome may be provided as a string, or list of strings for
        multiple chromosomes. Desired copy number may be provided as an integer, or
        list of integers for multiple copy numbers. Leaving either argument empty
        defaults to not removing any CNVs by that argument. Example: filter_cnvs(
        chromosomes=['6', '12'], copy_numbers=3) will return all CNVs for this patient
        that are on either chromosome 6 or 12 which are a single duplication, i.e.,
        copy number of 3.
        """
        if isinstance(chromosomes, str):
            chromosomes = {chromosomes}
        if isinstance(copy_numbers, int):
            copy_numbers = {copy_numbers}
        cnvs = [cnv for cnv in self.cnvs
                if (chromosomes is None or cnv.chromosome in chromosomes)
                and (copy_numbers is None or cnv.copy_number in copy_numbers)]
        return cnvs

    def get_median_cnv_position(
            self,
            chromosome: str,
            copy_numbers: Optional[Union[int, List[int]]] = None
            ) -> Optional[int]:
        """Calculate median position of all CNVs on a given chromosome.

        Median position is calculated as the median center point of all patient CNVs.
        CNVs can be filtered by specifying desired copy number.
        """
        cnvs = self.filter_cnvs(chromosome, copy_numbers)
        if not cnvs:
            return None
        cnv_median = int(median([(cnv.range.start + cnv.range.stop)/2 for cnv in cnvs]))
        return cnv_median

    def get_mean_cnv_position(
            self,
            chromosome: str,
            copy_numbers: Optional[Union[int, List[int]]] = None) -> Optional[float]:
        """Calculate mean position of all CNVs on a given chromosome.

        Mean position is calculated as the mean center point of all patient CNVs.
        CNVs can be filtered by specifying desired copy number.
        """
        cnvs = self.filter_cnvs(chromosome, copy_numbers)
        if not cnvs:
            return None
        cnv_mean = mean([(cnv.range.start + cnv.range.stop)/2 for cnv in cnvs])
        return cnv_mean

    def get_affected_ranges(
            self,
            chromosomes: Optional[Union[str, List[str]]] = None,
            copy_numbers: Optional[Union[int, List[int]]] = None) -> Dict[str, range]:
        """Return all ranges affected by a CNV.

        Overlapping CNVs are combined into a single range. CNVs can be filtered by
        specifying desired copy number and/or chromosomes.
        """
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        ranges = defaultdict(list)
        for cnv in cnvs:
            ranges[cnv.chromosome].append(cnv.range)
        ranges = {chromosome: merge_range_list(cnv_ranges)
                  for chromosome, cnv_ranges in ranges.items()}
        return ranges

    def assign_genes_to_cnvs(self, gene_set_obj: GeneSet) -> None:
        """Assign overlapping genes to each patient CNV based on provided GeneSet
        object."""
        for cnv in self.cnvs:
            cnv.genes = gene_set_obj.get_locus(cnv.chromosome, cnv.range.start,
                                               cnv.range.stop)

    def get_all_genes(self, chromosomes: Optional[Union[str, List[str]]] = None,
                      copy_numbers: Optional[Union[int, List[int]]] = None) -> Set[Gene]:
        """Return all genes overlapping selected CNVs.

        CNVs can be filtered by specifying desired copy number and/or chromosomes.
        """
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        all_genes = {gene for cnv in cnvs for gene in cnv.genes}
        return all_genes

    def get_all_HI_genes(self, chromosomes: Optional[Union[str, List[str]]] = None,
                         copy_numbers: Optional[Union[int, List[int]]] = None,
                         pLI_threshold: float = 0.9,
                         HI_threshold: float = 10,
                         phaplo_threshold: float = .86,
                         mode="confirm") -> Set[Gene]:
        """Return all haploinsufficient genes overlapping selected CNVs.

        The definition of a predicted haploinsufficient gene can be altered by
        adjusting the relevant score thresholds. CNVs can be filtered by specifying
        desired copy number and/or chromosomes.
        """
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        hi_genes = {gene for cnv in cnvs for gene in cnv.genes
                    if is_haploinsufficient(gene, pLI_threshold, HI_threshold,
                                            phaplo_threshold, mode)}
        return hi_genes

    def get_all_dominant_genes(
            self,
            chromosomes: Optional[Union[str, List[str]]] = None,
            copy_numbers: Optional[Union[int, List[int]]] = None) -> Set[Gene]:
        """Return all dominant-effect genes overlapping selected CNVs.

        Dominant-effect genes are labelled during GeneSet creation and are
        user-defined. CNVs can be filtered by specifying desired copy number and/or
        chromosomes.
        """
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        dom_genes = {gene for cnv in cnvs for gene in cnv.genes
                     if gene.dominant}
        return dom_genes

    # TODO: Needs to be updated to reflect Boolean responses, and to somehow allow
    #  searching for NA responses.
    # def filter_phenotypes_by_response(self, responses: Union[bool, List[str]]) -> Set[
    #     Term]:
    #     if isinstance(responses, bool):
    #         responses = {responses}
    #     filtered = {term for term, response in self.phenotypes.items()
    #                 if response in responses}
    #     return filtered


class PatientDatabase:
    """Database containing all Patient objects."""

    def __init__(self, patients):
        self.patients = patients
        self.index = self._make_index()
        self.cnvs = self._organize_cnvs()
        self.size = len(self.patients)
        self.phenotypes = {pheno for patient in self for pheno in patient.phenotypes}

    def __len__(self):
        return self.size

    def __getitem__(self, key):
        """Get patient by patient ID."""
        return self.patients[key]

    def __iter__(self):
        self.__iteri__ = 0
        return self

    def __next__(self):
        if self.__iteri__ == self.size:
            raise StopIteration
        result = self.patients[self.index[self.__iteri__]]
        self.__iteri__ += 1
        return result

    def _make_index(self):
        """Make database index for iteration purposes."""
        index = dict(enumerate(self.patients))
        return index

    def _organize_cnvs(self):
        """Get and organize CNVs from all patients."""
        cnvs = sorted(
            [cnv for patient in self.patients.values() for cnv in patient.cnvs],
            key=lambda x: x.start
            )
        cnv_dict = {chromosome: [] for chromosome in {cnv.chromosome for cnv in cnvs}}
        for cnv in cnvs:
            cnv_dict[cnv.chromosome].append(cnv)
        return cnv_dict

    def list_ids(self):
        """Return list of all patient IDs."""
        id_list = list(self.patients.keys())
        return id_list

    def add_patient(self, patient: Patient):
        """Add a Patient object to the patient database."""
        self.patients[patient.id] = patient
        self.index[max(self.index)+1] = patient.id
        self.size += 1
        self.phenotypes.update(set(patient.phenotypes))
        for cnv in patient.cnvs:
            self.cnvs[cnv.chromosome].append(cnv)

    def filter_cnvs(self,
                    chromosomes: Optional[Union[str, List[str]]] = None,
                    copy_numbers: Optional[Union[int, List[int]]] = None) -> List[CNV]:
        """Return a list of CNVs from all patients filtered by chromosome and copy number.

        Desired chromosome may be provided as a string, or list of strings for
        multiple chromosomes. Desired copy number may be provided as an integer, or
        list of integers for multiple copy numbers. Leaving either argument empty
        defaults to not removing any CNVs by that argument. Example: filter_cnvs(
        chromosomes=['6', '12'], copy_numbers=3) will return all CNVs for this patient
        that are on either chromosome 6 or 12 which are a single duplication, i.e.,
        copy number of 3.
        """
        if isinstance(chromosomes, str):
            chromosomes = {chromosomes}
        elif chromosomes is None:
            chromosomes = set(self.cnvs.keys())
        if isinstance(copy_numbers, int):
            copy_numbers = {copy_numbers}

        cnvs = [cnv for chromosome in chromosomes for cnv in self.cnvs[chromosome]
                if (copy_numbers is None or cnv.copy_number in copy_numbers)]
        return cnvs

    def get_median_cnv_position(
            self,
            chromosome: str,
            copy_numbers: Optional[Union[int, List[int]]] = None
            ) -> Optional[int]:
        """Calculate median position of CNVs from all patients on a given chromosome.

        Median position is calculated as the median center point of all patient CNVs.
        CNVs can be filtered by specifying desired copy number.
        """
        if chromosome not in self.cnvs:
            return None
        cnvs = self.filter_cnvs(chromosome, copy_numbers)
        cnv_median = median([(cnv.start + cnv.stop)/2
                             for cnv in cnvs])
        cnv_median = int(cnv_median)
        return cnv_median

    def get_mean_cnv_position(
            self,
            chromosome: str,
            copy_numbers: Optional[Union[int, List[int]]] = None
            ) -> Optional[int]:
        """Calculate mean position of CNVs from all patients on a given chromosome.

        Mean position is calculated as the mean center point of all patient CNVs.
        CNVs can be filtered by specifying desired copy number.
        """
        if chromosome not in self.cnvs:
            return None
        cnvs = self.filter_cnvs(chromosome, copy_numbers)
        cnv_mean = mean([(cnv.start + cnv.stop)/2
                         for cnv in cnvs])
        return cnv_mean

    def summarize_cnv_sizes(self,
                            chromosomes: Optional[Union[str, List[str]]] = None,
                            copy_numbers: Optional[Union[int, List[int]]] = None):
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        sizes = [cnv.length for cnv in cnvs]
        summary = pd.DataFrame(sizes, columns=["CNV Sizes"]).describe()
        return summary

    def summarize_hi_gene_counts(self,
                                 chromosomes: Optional[Union[str, List[str]]] = None,
                                 copy_numbers: Optional[Union[int, List[int]]] = None,
                                 pLI_threshold=0.9, HI_threshold=10,
                                 phaplo_threshold=0.86, mode="confirm"):
        params = dict(pLI_threshold=pLI_threshold, HI_threshold=HI_threshold,
                      phaplo_threshold=phaplo_threshold, mode=mode)
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        counts = [len([gene for gene in cnv.genes
                       if gene.is_haploinsufficient(**params)])
                  for cnv in cnvs]
        summary = pd.DataFrame(counts, columns=["HI Gene Count"]).describe()
        return summary

    def count_cnvs_with_dominant_genes(
            self,
            chromosomes: Optional[Union[str, List[str]]] = None,
            copy_numbers: Optional[Union[int, List[int]]] = None
            ):
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        de_count = sum(any([gene.dominant for gene in cnv.genes])
                       for cnv in cnvs)
        return de_count

    def summarize(self):
        """Print and return a summary of the patient cohort."""
        max_cnv = max((cnv for chrom in self.cnvs.values() for cnv in chrom),
                      key=lambda x: x.length)
        summary = f"""
        Patients: {self.size}
        CNVs: {len(self.cnvs)}
            Largest: {max_cnv.__repr__()}
        Phenotypes: {len(self.phenotypes)}
        """
        print(dedent(summary))
        summary_dict = dict(
            patients=self.size,
            cnvs=len(self.cnvs),
            largest_cnv=max_cnv,
            phenotypes=len(self.phenotypes)
            )
        return summary_dict

    def tabulate_phenotypes(self):
        """Tabulate phenotype info for all patients."""
        hpo_table = {patient.id: {hpo.id: response for hpo, response in
                                  patient.phenotypes.items()}
                     for patient in self}
        hpo_table = pd.DataFrame.from_dict(hpo_table, orient="index")
        hpo_table.index.name = "id"
        return hpo_table

    def tabulate_cnvs(self, chromosomes: Optional[Union[str, List[str]]] = None,
                      copy_numbers: Optional[Union[int, List[int]]] = None):
        """Tabulate CNV info for all patients."""
        cnvs = self.filter_cnvs(chromosomes, copy_numbers)
        cnv_table = [(cnv.id, cnv.chromosome, cnv.start, cnv.stop, cnv.copy_number)
                     for cnv in cnvs]
        cnv_table = pd.DataFrame(cnv_table, columns=["id", "chromosome", "start",
                                                     "stop", "copy_number"])
        return cnv_table
