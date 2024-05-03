#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Gene annotations.

@author: T.D. Medina
"""

import gzip
from collections import namedtuple
import importlib.resources as pkg_resources

import mygene
import pandas as pd
from venn import venn

from chr6_project.analysis.utilities import overlap, REFERENCE_CHR
from chr6_project import resources


# %% Annotation Types
class GeneAnnotation:
    """Record of a single gene annotation."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, gene_id, gene_name, transcript_id,
                 feature, seqname, start, end, strand,
                 score=".", frame=".", exon_id=None, exon_number=None,
                 source="."):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id

        self.feature = feature

        self.exon_id = exon_id
        self.exon_number = exon_number

        self.seqname = seqname
        self.start = start
        self.end = end

        self.strand = strand
        self.frame = frame
        self.score = score

        self.source = source

        self._hash = hash(str(self.__dict__))

    def __repr__(self):
        """Get official string representation."""
        attrs = [f"{x}={y}" for x, y in self.__dict__.items()
                 if y != "." and y is not None and not x.startswith("_")]
        string = f"{type(self).__name__}(" + ", ".join(attrs) + ")"
        return string

    def __hash__(self):
        """Hash str representation of object dictionary."""
        return self._hash

    def is_transcript(self):
        """Test if gene feature is transcript."""
        return self.feature == "transcript"


class Transcript(GeneAnnotation):
    """Gene sub-annotation of a single transcript.

    Designed to be nested inside a parent-level gene annotation.
    """
    def __init__(self, gene_id, gene_name, transcript_id,
                 seqname, start, end, strand,
                 score=".", source=".", annotations=None, **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "transcript", seqname, start, end, strand,
                         score,
                         source=source)
        self._annotations = annotations
        if self._annotations is None:
            self._annotations = []

        self._hash = self._make_hash()

    def _make_hash(self):
        return hash(self.__repr__())


class Exon(GeneAnnotation):
    """Gene sub-annotation of a single exon of a transcript.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "exon", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class UTR3(GeneAnnotation):
    """Gene sub-annotation of a transcript's UTR-3 region.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "3UTR", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class UTR5(GeneAnnotation):
    """Gene sub-annotation of a transcript's UTR-5 region.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "5UTR", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class StartCodon(GeneAnnotation):
    """Gene sub-annotation of a transcript's start codon.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, frame,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "start_codon", seqname, start, end, strand,
                         score, frame, exon_id, exon_number,
                         source)


class StopCodon(GeneAnnotation):
    """Gene sub-annotation of a transcript's stop codon.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, frame,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "stop_codon", seqname, start, end, strand,
                         score, frame, exon_id, exon_number,
                         source)


class CodingSequence(GeneAnnotation):
    """Gene sub-annotation of a transcript coding sequence.

    Designed to be nested inside a transcript annotation.
    """

    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, frame,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "CDS", seqname, start, end, strand,
                         score, frame, exon_id, exon_number,
                         source)


# %% GeneSet
class Gene:
    """Ensembl Gene object with unique gene ID."""

    def __init__(self, gene_id=None, gene_name=None, source=".",
                 seqname=None, start=None, end=None, strand=None,
                 transcripts=None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_symbol = None
        self.source = source

        self.seqname = seqname
        self.start = start
        self.end = end

        self.strand = strand

        self.transcripts = transcripts

        if self.transcripts is not None:
            self._set_attrs_from_transcripts()

        self.hi_score = None
        self.pli_score = None
        self.loeuf_score = None
        self.phaplo_score = None

        self.dominant = False

        self._hash = self._make_hash()

    def __repr__(self):
        """Official string representation."""
        locus = f"{self.seqname}:{self.start}-{self.end}"
        string = (f"Gene(gene_id={self.gene_id}, gene_name={self.gene_name}, "
                  f"source={self.source}, locus={locus}, "
                  f"strand={self.strand}, transcripts={len(self.transcripts)})")
        return string

    def __hash__(self):
        return self._hash

    def _set_attrs_from_transcripts(self):
        for attr in ["gene_id", "gene_name", "source", "seqname", "strand"]:
            attrs = {getattr(trans, attr) for trans in self.transcripts}
            if len(attrs) != 1:
                raise ValueError(f">1 {attr} detected from transcripts:"
                                 f"{self.gene_id}.{attr}: {attrs}")
            self.__setattr__(attr, list(attrs)[0])
        start = min([trans.start for trans in self.transcripts])
        end = max([trans.end for trans in self.transcripts])
        self.start = start
        self.end = end

    def _make_hash(self):
        string = self.__repr__()
        string += "\t".join([x.transcript_id for x in self.transcripts])
        return hash(string)

    def is_haploinsufficient(self, pLI_threshold=0.9, HI_threshold=10,
                             phaplo_threshold=0.86, mode="confirm"):
        """Check if gene has sufficient pLI, HI, or pHaplo score."""
        return is_haploinsufficient(self, pLI_threshold, HI_threshold,
                                    phaplo_threshold, mode)


class GeneSet:
    """Database of gene annotations."""

    def __init__(self, path=None):
        self.path = path
        self.genes = {}
        self.chromosomes = {}

        if self.path:
            self.genes = self.make_gene_set(path)
            self.organize_chromosomes()

        self.size = len(self.genes)

    def __len__(self):
        return self.size

    def make_gene_set(self, path):
        """Construct Genes and GeneSet objects from data file."""
        genes = self.read_annotations(path)
        genes = self.make_annotation_objects(genes)
        genes = self.make_genes(genes)
        return genes

    @staticmethod
    def read_annotations(file):
        """Read gene annotations from file."""
        with gzip.open(file) as infile:
            data = infile.readlines()

        data = [line.decode().rstrip(";\n").split("\t")
                for line in data]
        data = [line for line in data
                if line[0].lstrip("chr") in REFERENCE_CHR]

        for line in data:
            line[0] = line[0].lstrip("chr")
            line[3] = int(line[3])
            line[4] = int(line[4])
            ids = {}
            for field in line[-1].split(";"):
                field = field.strip().replace('"', "").split(" ")
                ids[field[0]] = field[1]
            line[-1] = ids

        return data

    @staticmethod
    def make_annotation_objects(data):
        """Construct gene sub-annotation objects."""
        fields = ["seqname", "source", "feature", "start", "end",
                  "score", "strand", "frame"]
        classes = {"transcript": Transcript, "exon": Exon,
                   "3UTR": UTR3, "5UTR": UTR5,
                   "start_codon": StartCodon, "stop_codon": StopCodon,
                   "CDS": CodingSequence}
        objects = set()
        while data:
            line = data.pop()
            line_dict = dict(zip(fields, line[:-1]))
            line_dict.update(line[-1])
            objects.add(classes[line_dict["feature"]](**line_dict))
        return objects

    @staticmethod
    def make_genes(data):
        """Make Gene objects from sub-annotations."""
        transcripts = {entry for entry in data
                       if isinstance(entry, Transcript)}
        data = data - transcripts
        transcripts = {trans.transcript_id: trans for trans in transcripts}

        while data:
            entry = data.pop()
            transcripts[entry.transcript_id]._annotations.append(entry)
        transcripts = list(transcripts.values())

        for transcript in transcripts:
            transcript._annotations.sort(key=lambda x: x.exon_number)

        genes = {trans.gene_id: [] for trans in transcripts}
        while transcripts:
            trans = transcripts.pop()
            genes[trans.gene_id].append(trans)
        genes = [Gene(transcripts=trans_list) for trans_list in genes.values()]
        for gene in genes:
            gene.transcripts.sort(key=lambda x: x.end)
            gene.transcripts.sort(key=lambda x: x.start)
        genes = {gene.gene_id: gene for gene in genes}
        return genes

    def organize_chromosomes(self):
        """Group genes by chromosome in a chromosome dictionary."""
        chrom_dict = {}
        # chrom_dict = {gene.seqname: {} for gene in self.genes}
        for gene in self.genes.values():
            if gene.seqname not in chrom_dict:
                chrom_dict[gene.seqname] = {}
            chrom_dict[gene.seqname][gene.gene_id] = gene
        self.chromosomes = chrom_dict

    def get_locus(self, seqname, start, stop=None):
        """Return all genes that intersect a base or range."""
        if seqname not in self.chromosomes:
            return []
        if stop is None:
            stop = start
        query = range(start, stop + 1)
        results = []
        for gene in self.chromosomes[seqname].values():
            if overlap(query, range(gene.start, gene.end+1)):
                results.append(gene)
        return results

    def add_pLI_scores(self, gnomad_info):
        """Add pLI and LOEUF scores to genes from gnomad data."""
        for gene, info in gnomad_info.items():
            if gene in self.genes:
                self.genes[gene].pli_score = info.pLI
                self.genes[gene].loeuf_score = info.oe_lof_upper

    def add_HI_scores(self, HI_info):
        """Add HI scores to genes."""
        for gene, info in HI_info.items():
            if gene in self.genes:
                self.genes[gene].hi_score = info.HI_score

    def add_phaplo_scores(self, phaplo_data):
        """Add pHaplo scores to genes."""
        for gene, score in phaplo_data.items():
            if gene in self.genes:
                self.genes[gene].phaplo_score = score

    def annotate_dominant_genes(self, gene_list):
        for gene in gene_list:
            if gene in self.genes:
                self.genes[gene].dominant = True

    def get_all_gene_ids(self):
        """Get gene IDs from all genes."""
        return set(self.genes.keys())

    def plot_haploinsufficient_venn(self, chromosome, pLI_threshold=0.9,
                                    HI_threshold=10, loeuf_threshold=.25,
                                    phaplo_threshold=0.86):
        venn_dict = {
            "HI": {gene.gene_id for gene in self.chromosomes[chromosome].values()
                   if gene.hi_score is not None and gene.hi_score <= HI_threshold},
            "pLI": {gene.gene_id for gene in self.chromosomes[chromosome].values()
                    if gene.pli_score is not None and gene.pli_score >= pLI_threshold},
            # "LOEUF": {gene.gene_id for gene in self.chromosomes[chromosome].values()
            #           if gene.loeuf_score is not None and gene.loeuf_score <= loeuf_threshold},
            "pHaplo": {gene.gene_id for gene in self.chromosomes[chromosome].values()
                       if gene.phaplo_score is not None and gene.phaplo_score >= phaplo_threshold}
            }
        venn(venn_dict)


# %% Haploinsufficients
pLI_info = namedtuple("pLI_info", [
    "gene", "transcript", "obs_mis", "exp_mis", "oe_mis", "mu_mis",
    "possible_mis", "obs_mis_pphen", "exp_mis_pphen", "oe_mis_pphen",
    "possible_mis_pphen", "obs_syn", "exp_syn", "oe_syn", "mu_syn",
    "possible_syn", "obs_lof", "mu_lof", "possible_lof", "exp_lof",
    "pLI", "pNull", "pRec", "oe_lof", "oe_syn_lower", "oe_syn_upper",
    "oe_mis_lower", "oe_mis_upper", "oe_lof_lower", "oe_lof_upper",
    "constraint_flag", "syn_z", "mis_z", "lof_z", "oe_lof_upper_rank",
    "oe_lof_upper_bin", "oe_lof_upper_bin_6", "n_sites", "classic_caf",
    "max_af", "no_lofs", "obs_het_lof", "obs_hom_lof", "defined", "p",
    "exp_hom_lof", "classic_caf_afr", "classic_caf_amr", "classic_caf_asj",
    "classic_caf_eas", "classic_caf_fin", "classic_caf_nfe",
    "classic_caf_oth", "classic_caf_sas", "p_afr", "p_amr", "p_asj", "p_eas",
    "p_fin", "p_nfe", "p_oth", "p_sas", "transcript_type", "gene_id",
    "transcript_level", "cds_length", "num_coding_exons", "gene_type",
    "gene_length", "exac_pLI", "exac_obs_lof", "exac_exp_lof",
    "exac_oe_lof", "brain_expression", "chromosome",
    "start_position", "end_position"
    ])

HI_info = namedtuple("HI_info", ["gene_symbol", "chromosome", "start",
                                 "stop", "HI_score"])


def read_gnomad_pli_data(file, pLI_max=1):
    data = pd.read_csv(file, sep="\t")
    data_objs = {}
    for i, row in data.iterrows():
        if row["pLI"] > pLI_max:
            continue
        if pd.isna(row["oe_lof_upper"]):
            row["oe_lof_upper"] = None
        pli_obj = pLI_info(*row)
        data_objs[pli_obj.gene_id] = pli_obj
    return data_objs


def read_HI_gene_data(file, HI_max=100):
    """Read HI gene information file."""
    with open(file) as infile:
        infile.readline()
        data = infile.readlines()
    data = [x.lstrip("chr").rstrip("\n").split("\t") for x in data]
    data = [[x[3].split("|")[0], x[0], int(x[1]), int(x[2]),
            float(x[3].split("|")[-1].rstrip("%"))]
            for x in data]
    data = {x[0]: HI_info(*x) for x in data}
    data = {x: y for x, y in data.items() if y.HI_score < HI_max}
    return data


def read_phaplo_gene_data(file):
    """Read pHaplo gene information file."""
    with open(file) as infile:
        infile.readline()
        data = infile.readlines()
    data = [x.strip().split("\t") for x in data]
    data = {x[0]: float(x[1]) for x in data}
    return data


def read_dominant_gene_list(file):
    with open(file) as infile:
        gene_list = infile.readlines()
    gene_list = [gene.strip() for gene in gene_list]
    return gene_list


def lookup_gene_symbols(symbol_list):
    mg = mygene.MyGeneInfo()
    results = mg.querymany(symbol_list, scopes="symbol",
                           species="human", fields="ensembl.gene")
    results_final = {}
    results_bad = {}
    results_none = {}
    results_multi = {}
    for result in results:
        query = result["query"]
        if "notfound" in result and result["notfound"]:
            results_none[query] = result
            continue
        if "ensembl" not in result:
            results_bad[query] = result
            continue
        if isinstance(result["ensembl"], list):
            results_multi[query] = result
            continue
        if query not in results_final:
            results_final[query] = result
            continue
        if result["_score"] > results_final[query]["_score"]:
            results_final[query] = result
    results = {"final": results_final, "multi": results_multi,
               "none": results_none, "bad": results_bad}
    return results


def replace_HI_gene_ids(hi_genes, mygene_lookup):
    new_hi_genes = {}
    for gene_symbol, gene in hi_genes.items():
        results = find_symbol_in_results(gene_symbol, mygene_lookup)
        if results is None:
            continue
        for result in results:
            new_hi_genes[result] = gene
    return new_hi_genes


def make_hi_genes(file):
    data = read_HI_gene_data(file)
    lookup = lookup_gene_symbols(data.keys())
    data = replace_HI_gene_ids(data, lookup)
    return data


# %% MyGene
def find_symbol_in_results(symbol, results):
    if symbol in results["final"]:
        return [results["final"][symbol]["ensembl"]["gene"]]
    if symbol in results["multi"]:
        return [gene["gene"] for gene in results["multi"][symbol]["ensembl"]]
    if symbol in results["bad"] or results["none"]:
        return None


def symbol_lookup_multi(mygene_instance, gene_symbols):
    results = mygene_instance.querymany(gene_symbols, scopes="symbol",
                                        species="human", fields="ensembl.gene")
    results_final = {}
    results_bad = {}
    results_none = {}
    results_multi = {}
    for result in results:
        query = result["query"]
        if "notfound" in result and result["notfound"]:
            results_none[query] = result
            continue
        if "ensembl" not in result:
            results_bad[query] = result
            continue
        if isinstance(result["ensembl"], list):
            results_multi[query] = result
            continue
        if query not in results_final:
            results_final[query] = result
            continue
        if result["_score"] > results_final[query]["_score"]:
            results_final[query] = result
    results = {"final": results_final, "multi": results_multi,
               "none": results_none, "bad": results_bad}
    return results


# TODO: Finish this to fix missing genes.
def symbol_relookup_missing(mygene_lookup, mygene_instance):
    missing = mygene_lookup["none"]
    results = {}
    for symbol in missing:
        result = mygene_instance.query(symbol, scope="symbol",
                                       species="human", fields="ensembl.gene")
        # ...
    return results


# %% Helper Functions
def is_haploinsufficient(gene, pLI_threshold=0.9, HI_threshold=10,
                         phaplo_threshold=0.86, mode="confirm"):
    """Check if gene has sufficiently low pLI or HI score."""
    pli_pass = gene.pli_score is not None and gene.pli_score >= pLI_threshold
    hi_pass = gene.hi_score is not None and gene.hi_score <= HI_threshold
    # loeuf_pass = gene.loeuf_score is not None and gene.loeuf_score <= loeuf_threshold
    phaplo_pass = gene.phaplo_score is not None and gene.phaplo_score >= phaplo_threshold
    if mode == "all":
        return all([hi_pass, pli_pass, phaplo_pass])
    elif mode == "any":
        return any([hi_pass, pli_pass, phaplo_pass])
    elif mode == "2":
        return (hi_pass + pli_pass + phaplo_pass) >= 2
    elif mode == "confirm":
        none_count = [gene.pli_score, gene.hi_score, gene.phaplo_score].count(None)
        if none_count == 0:
            return (hi_pass + pli_pass + phaplo_pass) >= 2
        elif none_count == 1:
            return (hi_pass + pli_pass + phaplo_pass) >= 1
        else:
            return False


# %% Main
def _read_defaults():
    geneset = list(pkg_resources.path(resources, "hg19.ensGene.chr6.gtf.gz").gen)[0]
    pli = list(pkg_resources.path(resources, "gnomad.v2.1.1.lof_metrics.by_gene.6.tsv").gen)[0]
    hi = list(pkg_resources.path(resources, "HI_Predictions.v3.chr6.bed").gen)[0]
    phaplo = list(pkg_resources.path(resources, "phaplo_scores.chr6.tsv").gen)[0]
    dom = list(pkg_resources.path(resources, "dominant_genes.txt").gen)[0]
    pli_genes = read_gnomad_pli_data(pli)
    phaplo_genes = read_phaplo_gene_data(phaplo)
    hi_genes = make_hi_genes(hi)
    dominant_genes = read_dominant_gene_list(dom)
    return geneset, pli_genes, hi_genes, phaplo_genes, dominant_genes


def read_geneset_gtf(geneset_gtf):
    """Load geneset."""
    _, pli_genes, hi_genes, phaplo_genes, dominant_genes = _read_defaults()
    geneset = GeneSet(geneset_gtf)
    geneset.add_pLI_scores(pli_genes)
    geneset.add_HI_scores(hi_genes)
    geneset.add_phaplo_scores(phaplo_genes)
    geneset.annotate_dominant_genes(dominant_genes)
    return geneset


def make_geneset():
    """Load geneset."""
    geneset, pli_genes, hi_genes, phaplo_genes, dominant_genes = _read_defaults()
    geneset = GeneSet(geneset)
    geneset.add_pLI_scores(pli_genes)
    geneset.add_HI_scores(hi_genes)
    geneset.add_phaplo_scores(phaplo_genes)
    geneset.annotate_dominant_genes(dominant_genes)
    return geneset
