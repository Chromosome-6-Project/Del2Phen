#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chromosome 6 Project - Gene annotations.

@author: T.D. Medina
"""

import gzip
from utilities import overlap, REFERENCE_CHR


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
    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "exon", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class UTR3(GeneAnnotation):
    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "3UTR", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class UTR5(GeneAnnotation):
    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "5UTR", seqname, start, end, strand,
                         score, ".", exon_id, exon_number,
                         source)


class StartCodon(GeneAnnotation):
    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, frame,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "start_codon", seqname, start, end, strand,
                         score, frame, exon_id, exon_number,
                         source)


class StopCodon(GeneAnnotation):
    def __init__(self, gene_id, gene_name, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, frame,
                 score=".", source=".", **kwargs):
        super().__init__(gene_id, gene_name, transcript_id,
                         "stop_codon", seqname, start, end, strand,
                         score, frame, exon_id, exon_number,
                         source)


class CodingSequence(GeneAnnotation):
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
        self.source = source

        self.seqname = seqname
        self.start = start
        self.end = end

        self.strand = strand

        self.transcripts = transcripts

        if self.transcripts is not None:
            self._set_attrs_from_transcripts()

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
                raise ValueError(f">1 {attr} detected from transcripts.")
            self.__setattr__(attr, list(attrs)[0])
        start = min([trans.start for trans in self.transcripts])
        end = min([trans.end for trans in self.transcripts])
        self.start = start
        self.end = end

    def _make_hash(self):
        string = self.__repr__()
        string += "\t".join([x.transcript_id for x in self.transcripts])
        return hash(string)


class GeneSet:
    """Database of gene annotations."""

    def __init__(self, path=None):
        self.path = path
        self.genes = []

        if self.path:
            self.genes = self.make_gene_set(path)

    def make_gene_set(self, path):
        flat_annotes = self.read_annotations(path)
        annote_objects = self.make_annotation_objects(flat_annotes)
        genes = self.make_genes(annote_objects)
        organized_genes = self.organize_genes_by_seqname(genes)
        return organized_genes

    @staticmethod
    def read_annotations(file):
        """Read gene annotations from file.

        Example line:
            ensGene	exon	95124	95454	.	+	.	gene_id "ENSG00000271530"; transcript_id "ENST00000604449"; exon_number "1"; exon_id "ENST00000604449.1"; gene_name "ENSG00000271530";
        """
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
        transcripts = {entry for entry in data if isinstance(entry, Transcript)}
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
        return genes

    @staticmethod
    def organize_genes_by_seqname(genes):
        chrom_dict = {gene.seqname: [] for gene in genes}
        for gene in genes:
            chrom_dict[gene.seqname].append(gene)
        return chrom_dict

    def get_locus(self, seqname, start, stop=None):
        """Return all genes that intersect a base or range."""
        if seqname not in self.genes:
            return []
        if stop is None:
            stop = start
        query = range(start, stop + 1)
        results = []
        for gene in self.genes[seqname]:
            if overlap(query, range(gene.start, gene.end+1)):
                results.append(gene)
        return results


# def main():
#     """Load geneset."""
#     geneset = GeneSet("C:/Users/tyler/Documents/Chr6/hg19.ensGene.gtf.gz")
#     return geneset


# if __name__ == "__main__":
#     geneset = main()
