"""The Chromosome 6 Project Data Management."""
import csv
import gzip
from collections import namedtuple
from math import log
import matplotlib.pyplot as plt
import numpy as np
from pronto import Ontology
import mygene

REFERENCE_CHR = [str(i) for i in range(1, 23)] + ["X", "Y"]


# %% Utilities
def jaccard(set1, set2):
    """Calculate Jaccard Index of two sets."""
    if isinstance(set1, (int, float)) and isinstance(set2, (int, float)):
        if set1 == 0 or set2 == 0:
            return 0
        return min([set1, set2]) / max([set1, set2])
    set1 = set(set1)
    set2 = set(set2)
    union = set1 | set2
    if len(union) == 0:
        return 0, set()
    intersect = set1 & set2
    jaccard_index = len(intersect) / len(union)
    return jaccard_index, intersect


def cnv_overlap(cnv1, cnv2):
    """Check if two CNVs overlap each other."""
    if cnv1.chromosome != cnv2.chromosome:
        return False
    return overlap(cnv1.range, cnv2.range)


def overlap(range1: range, range2: range):
    """Test if two ranges intersect."""
    # Needs to be <, not <=, because of range half-open notation.
    # e.g. range(1,3) should not overlap range(3, 5) because
    # range(1,3) is [1,2] and range(3,5) is [3,4].
    if range1.start < range2.stop and range2.start < range1.stop:
        return True
    return False


def get_range_intersect(range1: range, range2: range):
    """Make range from the intersect of two ranges."""
    if not overlap(range1, range2):
        return range(0)
    start = max(range1.start, range2.start)
    stop = min(range1.stop, range2.stop)
    return range(start, stop)


def merge_ranges(range1: range, range2: range):
    """Make range from the union of two ranges."""
    if not overlap(range1, range2):
        raise ValueError("Ranges do not overlap.")
    merged = range(min(range1.start, range1.start),
                   max(range1.stop, range2.stop))
    return merged


def merge_range_list(ranges):
    """Merge all overlapping ranges in a list of ranges."""
    ranges_copy = sorted(ranges.copy(), key=lambda x: x.stop)
    ranges_copy = sorted(ranges_copy, key=lambda x: x.start)
    merged = []

    while ranges_copy:
        range1 = ranges_copy[0]
        del ranges_copy[0]
        merges = []
        for i, range2 in enumerate(ranges_copy):
            if overlap(range1, range2):
                range1 = merge_ranges(range1, range2)
                merges.append(i)
        merged.append(range1)
        for i in reversed(merges):
            del ranges_copy[i]
    merged = sorted(merged, key=lambda x: x.start)

    return merged


def length_of_range_intersects(ranges):
    """Calculate total length of all range intersects in list of ranges."""
    ranges_copy = sorted(ranges.copy(), key=lambda x: x.stop)
    ranges_copy = sorted(ranges_copy, key=lambda x: x.start)
    intersects = []

    while ranges_copy:
        range1 = ranges_copy[0]
        del ranges_copy[0]
        for range2 in ranges_copy:
            if not overlap(range1, range2):
                continue
            intersects.append(get_range_intersect(range1, range2))
    intersects = merge_range_list(intersects)
    total = sum([len(x) for x in intersects])
    return total


def is_patient(patient):
    """Check if object is a Patient object.

    HI_Gene objects will still return False.
    """
    if not isinstance(patient, HI_Gene) and isinstance(patient, Patient):
        return True
    return False


def is_gene(gene):
    """Check if object is a HI_Gene object."""
    if isinstance(gene, HI_Gene):
        return True
    return False


# %% Data Classes

sequence_contig = namedtuple("Sequence_Contig",
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
            lines[i] = sequence_contig(*line)
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


class PatientIntersect:
    """Record of similarity between two patients."""

    def __init__(self, patient1, patient2,
                 genes, gene_count, gene_similarity,
                 hpos, hpo_count, hpo_similarity):
        self.patients = [patient1, patient2]
        self.genes = genes
        self.gene_count = gene_count
        self.gene_similarity = gene_similarity
        self.hpos = hpos
        self.hpo_count = hpo_count
        self.hpo_similarity = hpo_similarity

    def __repr__(self):
        """Get official string representation."""
        string = (f"PatientIntersect(patients={self.patients}, "
                  f"gene_count={self.gene_count}, hpo_count={self.hpo_count})")
        return string

    def __str__(self):
        """Get pretty-printing string representation."""
        string = (f"{self.patients[0]} vs. {self.patients[1]}:\n"
                  f"  Genes: {self.gene_count}\n"
                  f"  HPO terms: {self.hpo_count}")
        return string

    def gene_info(self):
        """Get gene portion of intersect."""
        return self.genes, self.gene_count, self.gene_similarity

    def hpo_info(self):
        """Get HPO portion of intersect."""
        return self.hpos, self.hpo_count, self.hpo_similarity


class PatientIntersect2:
    """Record of similarity between two patients."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, patient_1, patient_2,
                 length_compare, loci_compare,
                 gene_compare, hpo_compare):
        self.patients = [patient_1, patient_2]

        self.length_similarity = length_compare

        self.loci_similarity = loci_compare[0]
        self.loci_shared_size = loci_compare[1]

        self.gene_similarity = gene_compare[0]
        self.genes = gene_compare[1]
        self.gene_count = len(gene_compare[1])

        self.hpo_similarity = hpo_compare[0]
        self.hpos = hpo_compare[1]
        self.hpo_count = len(hpo_compare[1])

    def __repr__(self):
        """Get official string representation."""
        string = (f"PatientIntersect(patients={self.patients}, "
                  f"length_similarity={self.length_similarity}, "
                  f"loci_similarity={self.loci_similarity}, "
                  f"gene_similarity={self.gene_similarity}, "
                  f"hpo_similarity={self.hpo_similarity})")
        return string

    def __str__(self):
        """Get pretty-printing string representation."""
        string = (f"Similarities of {self.patients[0]} vs. {self.patients[1]}:\n"
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


# TODO: Add a method to add a patient.
class ComparisonTable:
    """Data object holding all patient vs. patient comparisons."""

    def __init__(self, patient_db):
        self.patient_db = patient_db
        self.raw = self.compare_patients()
        self.index = self.make_index()
        self.array = self.make_array2()
        self.size = len(self.index)

        self.__iteri__ = 0
        self.__iterj__ = 0

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

    @classmethod
    def compare_all(cls, patient_1, patient_2):
        """Compare all metrics between two patients."""
        length_compare = cls.compare_length(patient_1, patient_2)
        loci_compare = cls.compare_loci(patient_1, patient_2)
        gene_compare = cls.compare_genes(patient_1, patient_2)
        hpo_compare = cls.compare_hpos(patient_1, patient_2)
        comparison = PatientIntersect2(
            patient_1.id, patient_2.id,
            length_compare, loci_compare, gene_compare, hpo_compare
            )
        return comparison

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
        if (is_gene(patient_1) or is_gene(patient_2)) and jaccard_index > 0:
            jaccard_index = 1.0
        return jaccard_index, intersect

    @staticmethod
    def compare_hpos(patient_1, patient_2):
        """Compare HPO terms between two patients."""
        jaccard_index, intersect = jaccard(patient_1.hpo, patient_2.hpo)
        return jaccard_index, intersect

    def make_array(self, pack=True):
        """Structure patient comparisons as a numpy array."""
        array = []
        for pid in self.index:
            values = []
            for p2 in self.index:
                if p2 in self.raw[pid]:
                    lookup = self.raw[pid][p2]
                else:
                    lookup = self.raw[p2][pid]
                if pack:
                    lookup = PatientIntersect(pid, p2, *lookup)
                values.append(lookup)
            array.append(values)
        if pack:
            array = np.array(array)
        else:
            array = np.array(
                array,
                dtype=[
                    ("genes", object), ("genecount", int), ("genesim", float),
                    ("hpos", object), ("hpocount", int), ("hposim", float)
                    ]
                )
        return array

    def make_array2(self):
        array = []
        for pid in self.index:
            values = []
            for p2 in self.index:
                if p2 in self.raw[pid]:
                    values.append(self.raw[pid][p2])
                else:
                    values.append(self.raw[p2][pid])
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
        if pid1 not in self.index:
            raise KeyError("ID not found.")
        return self.array[self.index[pid1]][self.index[pid2]]

    def write_all_comparisons(self, outfile, patients_only=True):
        properties = ["length_similarity",
                      "loci_similarity", "loci_shared_size",
                      "gene_similarity", "gene_count",
                      "hpo_similarity", "hpo_count"]
        write_me = ["\t".join(["patient1", "patient2"] + properties) + "\n"]
        for intersect in self:
            p1, p2 = [self.patient_db[p] for p in intersect.patients]
            if p1 == p2:
                continue
            if (is_gene(p1) or is_gene(p2)) and patients_only:
                continue
            this_intersect = "\t".join(
                [f"{intersect.patients[0]}", f"{intersect.patients[1]}"]
                + [f"{intersect.__getattribute__(prop)}" for prop in properties]
                ) + "\n"
            write_me.append(this_intersect)
        with open(outfile, "w") as out:
            out.writelines(write_me)


# XXX: This probably doesn't work anymore.
def write_comparison_table(table, patients, out, self_match="size"):
    """Write comparison table to file. Deprecated."""
    header = "PatientID"
    writer = []
    for pid in table:
        header += f",{pid}"
        string = pid
        for p2 in table:
            if pid == p2:
                if self_match == "size":
                    string += f",{len(patients[pid].all_genes())}"
                else:
                    string += f",{self_match}"
            elif p2 in table[pid]:
                string += f",{table[pid][p2]}"
            else:
                string += f",{table[p2][pid]}"
        string += "\n"
        writer.append(string)
    writer = [header + "\n"] + writer
    with open(out, "w") as outfile:
        outfile.writelines(writer)


class GeneAnnotation:
    """Record of a single gene annotation."""
    # pylint: disable=too-many-instance-attributes

    def __init__(self, chromosome, source, feature, start, end, score, strand,
                 frame, gene_id=None, gene_name=None, transcript_id=None,
                 exon_id=None, exon_number=None):
        self.chromosome = chromosome
        self.source = source
        self.feature = feature
        self.range = range(start, end + 1)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.exon_id = exon_id
        self.exon_number = exon_number
        self._hash = hash(str(self.__dict__))

    def __repr__(self):
        """Get official string representation."""
        string = ("GeneAnnotation("
                  f"gene_id={self.gene_id}, "
                  f"locus={self.chromosome}:{self.range.start}-{self.range.stop - 1}"
                  f"")
        return string

    def __hash__(self):
        """Hash str represenation of object dictionary."""
        return self._hash

    def is_transcript(self):
        """Test if gene feature is transcript."""
        return self.feature == "transcript"


class GeneSet:
    """Database of gene annotations."""

    def __init__(self, file=None, genes_only=True):
        self.path = file
        self.genes = []
        self.gene_ids = []

        if self.path:
            self.genes = self.read_genes(self.path)

    @staticmethod
    def read_genes(file):
        """Read gene annotations from file."""
        with gzip.open(file) as infile:
            data = infile.readlines()
            # i = 0
            # while i < 50:
            #     data.append(infile.readline())
            #     i += 1
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
        data = [GeneAnnotation(*line[:-1], **line[-1]) for line in data]
        data_map = {chrom: [] for chrom in {gene.chromosome for gene in data}}
        for gene in data:
            data_map[gene.chromosome].append(gene)
        return data_map

    def get_locus(self, chromosome, start, stop=None):
        """Return all gene annotations that interect a base or range."""
        if stop is None:
            stop = start
        query = range(start, stop + 1)
        results = []
        for gene in self.genes[chromosome]:
            if overlap(query, gene.range):
                results.append(gene)
        return results


class DataManager:
    """Collection of Chromosome 6 Project data handling methods."""

    def __init__(self):
        pass

    @classmethod
    def read_data(cls, data_file):
        """Import general CSV data."""
        data = []
        with open(data_file) as infile:
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

    @classmethod
    def fix_genotype_data(cls, data):
        """Apply all data fixes."""
        fixed = cls.trim_chromosome_names(data)
        fixed = cls.fix_hg19_positions(fixed)
        return fixed

    @staticmethod
    def make_patients(genotypes, phenotypes,
                      geneset=None, hpos=None, ontology=None):
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
                patient_hpos = []
                for hpo_id, value in hpos[patient.id].items():
                    if hpo_id == "label":
                        continue
                    if value == "T":
                        patient_hpos.append(ontology[hpo_id])
                patient.hpo = patient_hpos
        for patient in patients.values():
            if patient.genotypes:
                patient.origin = patient.genotypes[0]["origin info"]
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
    def make_UCSC_browser_tracks(patients, out, filter_chrom=None):
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
    def read_HI_genes(file):
        """Read HI gene information file."""
        with open(file) as infile:
            infile.readline()
            data = infile.readlines()
        data = [x.lstrip("chr").rstrip("\n").split("\t") for x in data]
        data = {x[3].split("|")[0]: [x[0], int(x[1]), int(x[2]),
                                     float(x[3].split("|")[-1].rstrip("%"))]
                for x in data}
        return data

# !!!
    @classmethod
    def make_HI_objects(cls, hi_genes, geneset, gene_lookup):
        """Construct HI_Gene objects from HI gene info."""
        hi_gene_objs = {}
        for geneID, hi_gene in hi_genes.items():
            this_hi = HI_Gene(geneID)
            this_hi.genotypes.append(dict(zip(
                ["Chromosome", "Start positie in Hg19", "Stop positie in Hg19", "imbalance"],
                hi_gene[:3] + ["HI Gene"]
                )))
            this_hi.cnvs = this_hi.extract_cnvs()
            this_hi.identify_gene_overlaps(geneset, True)

            symbol_results = cls.find_symbol_in_results(geneID, gene_lookup)
            if symbol_results is not None:
                refine = [gene for gene in this_hi.cnvs[0].genes
                          if gene.gene_id in symbol_results]
                if not refine:
                    print(f"Warning: No refined matches for {geneID}.")
                else:
                    this_hi.cnvs[0].genes = refine
                    this_hi.refined = True

            this_hi.score = hi_gene[3]
            this_hi.origin = "HI Gene"
            hi_gene_objs[geneID] = this_hi
        return hi_gene_objs

    @staticmethod
    def find_symbol_in_results(symbol, results):
        if symbol in results["final"]:
            return [results["final"][symbol]["ensembl"]["gene"]]
        if symbol in results["multi"]:
            return [gene["gene"] for gene in results["multi"][symbol]["ensembl"]]
        if symbol in results["bad"] or results["none"]:
            return None

    @staticmethod
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

    # @staticmethod
    # def read_gnomad_pli_data(file):
    #     with open(file) as infile:
    #         reader = csv.DictReader(infile, delimiter="\t")
    #         data = [row for row in reader]
    #         # header = infile.readline()
    #         # data = infile.readlines()
    #     # header = header.rstrip("\n").split("\t")
    #     # data = [line.rstrip("\n").split("\t") for line in data]
    #     data = {x["gene"]: x for x in data}
    #     return data


    # @staticmethod
    # def symbol_lookup(mygene_instance, gene_symbol):
    #     results = mygene_instance.query(gene_symbol, fields="ensembl.gene", species="human")
    #     if len(results["hits"]) > 1:
    #         print(f"Warning: Multiple hits for {gene_symbol}")
    #         return None
    #     if len(results["hits"]) == 0:
    #         print(f"Warning: No results found for {gene_symbol}")
    #     return results["hits"][0]["ensembl"]["gene"]


class Cytoband:
    def __init__(self, chromosome, start, stop, band, stain):
        self.chr = chromosome
        self.locus = range(start, stop)
        self.band = band
        self.stain = stain

    def __str__(self):
        string = ("Cytoband("
                  f"locus='{self.chr}:{self.locus.start}-{self.locus.stop - 1}', "
                  f"band='{self.band}', "
                  f"stain='{self.stain}')")
        return string

    def __repr__(self):
        string = f"Cytoband({self.chr}:{min(self.locus)}-{max(self.locus)})"
        return string


class Cytomap:
    def __init__(self, file="C:/Users/tyler/Documents/cytoBand.txt"):
        self.path = file
        self.cytobands = self.make_cytobands(self.path)

    @classmethod
    def make_cytobands(cls, filepath):
        cytobands = cls.read_cytoband_file(filepath)
        cytobands = cls.split_cytobands_by_chr(cytobands)
        return cytobands

    @staticmethod
    def read_cytoband_file(file):
        with open(file) as infile:
            cytobands = infile.readlines()
        cytobands = [line.strip().split("\t") for line in cytobands]
        cytobands = [Cytoband(line[0].strip("chr"), int(line[1]), int(line[2]), line[3], line[4])
                     for line in cytobands]
        return cytobands

    @staticmethod
    def split_cytobands_by_chr(cytobands):
        cytomap = {chrom: [] for chrom in {cytoband.chr for cytoband in cytobands}}
        for cytoband in cytobands:
            cytomap[cytoband.chr].append(cytoband)
        return cytomap

    @staticmethod
    def split_cytomap_by_arm(cytomap):
        arm_map = {key: {"p": [], "q": []} for key in cytomap.keys()}
        for chrom_list in cytomap.values():
            for cytoband in chrom_list:
                arm_map[cytoband.chr][cytoband.band[0]].append(cytoband)
        return arm_map

    def get_band(self, chromosome, coordinate):
        for band in self.cytobands[chromosome]:
            if coordinate in band.locus:
                return band
        return None

    def get_bands(self, chromosome, range_: range):
        bands = []
        for band in self.cytobands[chromosome]:
            if overlap(range_, band.locus):
                bands.append(band)
        return bands


class CNV:
    def __init__(self, chromosome, start, stop, change, ID=None, genes=None):
        self.id = ID
        self.chromosome = str(chromosome)
        self.range = range(start, stop + 1)
        self.change = change
        self.length = len(self.range)
        self.genes = genes

    def __repr__(self):
        string = (
            f"CNV({self.change}:"
            f"{self.chromosome}:{self.range.start}-{self.range.stop})"
            )
        return string

    def __str__(self):
        string = ("CNV:\n"
                  f"  Change = {self.change}\n"
                  f"  Locus = {self.chromosome}:"
                  f"{self.range.start}-{self.range.stop}\n"
                  f"  Length = {self.length:,} bp")
        return string


# TODO: Add a method to add a patient.
class PatientDatabase:
    def __init__(self, patients):
        self.patients = patients
        self.index = self.make_index()
        self.cnvs = self.organize_cnvs()
        self.size = len(self.patients)
        print(self.summary())

    def __getitem__(self, key):
        return self.patients[key]

    def __iter__(self):
        """Initialize iterable."""
        self.__iteri__ = 0
        return self

    def __next__(self):
        """Iterate 2-axis iterable."""
        if self.__iteri__ == self.size:
            raise StopIteration
        result = self.patients[self.index[self.__iteri__]]
        self.__iteri__ += 1
        return result

    def make_index(self):
        index = dict(enumerate(self.patients))
        return index

    def organize_cnvs(self):
        cnvs = sorted([cnv for patient in self.patients.values() for cnv in patient.cnvs], key=lambda x: x.range.start)
        cnv_dict = {chromosome: [] for chromosome in {cnv.chromosome for cnv in cnvs}}
        for cnv in cnvs:
            cnv_dict[cnv.chromosome].append(cnv)
        return cnv_dict

    def summary(self):
        len_cnv = len([cnv for chromosome in self.cnvs.values() for cnv in chromosome])
        summary_txt = (f"Patients: {len(self.patients)}\n"
                       f"CNVs: {len_cnv}")
        return summary_txt

    def make_phenotype_table(self):
        all_hpos = sorted(list({hpo.id for patient in self for hpo in patient.hpo}))
        all_hpos = {name: pos for pos, name in enumerate(all_hpos)}
        hpo_count = len(all_hpos)

        header = ["patient"] + list(all_hpos.keys())
        rows = {}

        for patient in self:
            if is_gene(patient):
                continue
            row = [0] * hpo_count
            for hpo in patient.hpo:
                row[all_hpos[hpo.id]] = 1
            row.insert(0, patient.id)
            row = dict(zip(header, row))
            rows[row["patient"]] = row
        return rows

    def make_gene_table(self):
        all_genes = sorted(list({gene.gene_id for patient in self for gene in patient.all_genes()}))
        all_genes = {name: pos for pos, name in enumerate(all_genes)}
        gene_count = len(all_genes)

        header = ["patient"] + list(all_genes.keys())
        rows = {}

        for patient in self:
            if is_gene(patient):
                continue

            row = [0] * gene_count
            for gene in patient.all_genes():
                row[all_genes[gene.gene_id]] = 1
            row.insert(0, patient.id)

            row = dict(zip(header, row))
            rows[row["patient"]] = row
        return rows

    def make_genotypes_table(self, absolute_pos=True, genome_dict=None):
        if absolute_pos and genome_dict is None:
            raise ValueError("Asbolute positions require a reference genome dictionary.")
        header  = ["patient", "start", "length"]

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

    def make_combination_table(self, genotypes=True, genes=True, phenotypes=False,
                               absolute_pos=True, genome_dict=None):

        # geno_table = {}
        # gene_table = {}
        # pheno_table = {}

        # if genotypes:
        #     geno_table = self.make_genotypes_table(absolute_pos, genome_dict)
        # if genes:
        #     gene_table = self.make_gene_table()
        # if phenotypes:
        #     pheno_table = self.make_phenotype_table()

        subtables = []
        if genotypes:
            subtables.append(self.make_genotypes_table(absolute_pos, genome_dict))
        if genes:
            subtables.append(self.make_gene_table())
        if phenotypes:
            subtables.append(self.make_phenotype_table())

        names = set(subtables[0].keys())
        names.intersection_update(*[set(subtable.keys()) for subtable in subtables])
        # names = {name for subtable in subtables for name in subtable}
        # names = (set(geno_table.keys())
        #          | set(gene_table.keys())
        #          | set(pheno_table.keys()))

        rows = {}
        for name in names:
            row = {}
            for subtable in subtables:
                row.update(subtable[name])
            # if any([genotypes and name not in geno_table, genes and name not in gene_table, phenotypes and name not in pheno_table]):
            #     continue
            # row = {}
            # row.update(geno_table[name])
            # row.update(gene_table[name])
            # row.update(pheno_table[name])
            rows[name] = row
        return rows

    @staticmethod
    def write_table(table, out):
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
        self.hpo = []
        self.origin = None

    def __repr__(self):
        """Construct string representation."""
        return (f"ID: {self.id}\n"
                f"  Genotypes: {len(self.genotypes)}\n"
                f"  Phenotypes: {len(self.phenotypes)}")

    def extract_cnvs(self):
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

    def get_affected_ranges(self):
        ranges = {cnv.chromosome: [] for cnv in self.cnvs}
        for cnv in self.cnvs:
            ranges[cnv.chromosome].append(cnv.range)
        ranges = {chromosome: merge_range_list(cnv_ranges)
                  for chromosome, cnv_ranges in ranges.items()}
        return ranges

    def identify_gene_overlaps(self, gene_set, transcripts_only=True):
        for cnv in self.cnvs:
            results = gene_set.get_locus(cnv.chromosome,
                                         cnv.range.start,
                                         cnv.range.stop)

            cnv.genes = [x for x in results if not transcripts_only or x.feature == "transcript"]
            # self.affected_genes[cnv.__repr__()] = results
            # self.affected_gene_ids[cnv.__repr__()] = set([x.gene_id for x in results])

    def all_genes(self, transcripts_only=True):
        all_genes = {gene for cnv in self.cnvs for gene in cnv.genes
                     if gene.is_transcript() or not transcripts_only}
        return all_genes


class HI_Gene(Patient):
    def __init__(self, ID):
        super().__init__(ID)
        self.score = None
        self.refined = False

# %% Graph Stuff

Node = namedtuple("Node", ["id", "label", "group", "color",
                           "value", "title", "HI", "ranges"])
Edge = namedtuple("Edge", ["edge_id", "id1", "id2", "width",
                           "color", "title", "sim"])


def build_network_nodes2(comparison_table):
    """Build patient node objects from comparison table."""
    colors = {"Literature case report": "blue",
              "Parental uploaded array report": "pink",
              None: "black",
              "HI Gene": "yellow"}

    nodes = {}
    for patient in comparison_table.index:
        lookup = comparison_table.lookup(patient)
        group = comparison_table.patient_db[patient].origin
        ranges = []
        for cnv in comparison_table.patient_db[patient].cnvs:
            ranges.append(f"{cnv.chromosome}:{cnv.range.start}:{cnv.range.stop}")
        ranges = ";".join(ranges)
        color = colors[group]
        if group == "HI Gene":
            hi = comparison_table.patient_db[patient].score
            value = 2
            title = (f"<p>{patient}<br>"
                     f"Group: {group}<br>"
                     f"HI score: {hi}<br>"
                     f"Transcripts: {value}<br></p>")
        else:
            hi = 0
            for intersect in comparison_table.array[comparison_table.index[patient]]:
                patient2 = comparison_table.patient_db[intersect.patients[1]]
                if (is_gene(patient2)
                        and intersect.gene_count > 0
                        and patient2.score <= 2):
                    hi += 1
            value = lookup.gene_count
            title = (f"<p>{patient}<br>"
                     f"Group: {group}<br>"
                     f"Affected transcripts: {value}<br>"
                     f"HPO terms: {lookup.hpo_count}</p>")
        node = Node(patient, patient, group, color, value, title, hi, ranges)
        nodes[patient] = node
    return nodes


def write_network_nodes2(nodes, out, normalize=True):
    """Write patient nodes to CSV file."""
    writer = ["id,label,group,color,value,title,hi,ranges\n"]
    for node in sorted(list(nodes.values()), key=lambda x: x[0]):
        if normalize:
            node = list(node)
            node[4] = log(node[4] + 1, 2)
        writer.append(",".join([str(x) for x in node]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)


def build_network_edges2(comparison_table):
    """Build patient-patient edge objects from comparison table."""
    edges = []
    for intersect in comparison_table:
        if not intersect.gene_count or len(set(intersect.patients)) == 1:
            continue
        id1 = intersect.patients[0]
        id2 = intersect.patients[1]
        edge_id = f"{id1}_{id2}"
        gene_count = intersect.gene_count
        gene_sim = intersect.gene_similarity
        hpo_count = intersect.hpo_count
        hpo_sim = intersect.hpo_similarity

        color = "gray"
        title = (f"<p>{id1}---{id2}:<br>"
                 f"Shared genes: {gene_count} ({gene_sim:.2%})<br>"
                 f"Shared HPO terms: {hpo_count} ({hpo_sim:.2%})<br></p>")
        edges.append(Edge(edge_id, id1, id2, gene_count,
                          color, title, gene_sim))
    return edges


def write_network_edges2(edges, out, normalize=True):
    """Write patient-patient edges to CSV file."""
    writer = ["id,from,to,width,color,title,gene_sim\n"]
    for edge in edges:
        if normalize:
            edge = list(edge)
            edge[3] = log(edge[3], 2)
        writer.append(",".join([str(x) for x in edge]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)


# XXX: Probably deprecated.
def make_array(table):
    """Make array from comparisons."""
    array = []
    for pid in table:
        values = []
        for p2 in table:
            if pid == p2:
                values.append(1)
            elif p2 in table[pid]:
                values.append(table[pid][p2])
            else:
                values.append(table[p2][pid])
        array.append(values)
    array = np.array(array)
    return array


# XXX: Probably deprecated.
def gene_comparison_heatmap(table):
    """Plot heatmap of gene comparisons."""
    data_array = make_array(table)

    fig, ax = plt.subplots()
    im = ax.imshow(data_array)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(table)))
    ax.set_yticks(np.arange(len(table)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(list(table.keys()))
    ax.set_yticklabels(list(table.keys()))

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    # for i in range(len(vegetables)):
    #     for j in range(len(farmers)):
    #         text = ax.text(j, i, harvest[i, j],
    #                        ha="center", va="center", color="w")

    ax.set_title("Test Heat Map")
    fig.tight_layout()
    plt.show()


def similarity_vs_hpo_scatter(patient_comparison):
    points = []
    for intersect in patient_comparison:
        p1 = patient_comparison.patient_db[intersect.patients[0]]
        p2 = patient_comparison.patient_db[intersect.patients[1]]
        if is_gene(p1) or is_gene(p2) or p1.id == p2.id:
            continue
        points.append((intersect.gene_similarity, intersect.hpo_count))
    points = list(zip(*points))
    plt.scatter(*points)
    plt.ylabel("Shared HPO terms", fontsize=24)
    plt.xlabel("Shared genes (percent similarity)", fontsize=24)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)


def plot_individual_factors(patient_comparison):
    plotters = [x for x in patient_comparison if not is_gene(x)]
    plotters = [x for x in plotters if x.patients[0] != x.patients[1]]
    hpos = [x.hpo_count for x in plotters]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.scatter([x.length_similarity for x in plotters], hpos)
    ax2.scatter([x.loci_similarity for x in plotters], hpos)
    ax3.scatter([x.gene_similarity for x in plotters], hpos)
    fig.suptitle("Similarity vs. Shared HPO Terms.", fontsize=30)
    ax1.set_ylabel("Shared HPO terms", fontsize=24)
    ax1.set_xlabel("Length Similarity", fontsize=24)
    ax2.set_xlabel("Loci Similarity", fontsize=24)
    ax3.set_xlabel("Gene Similarity", fontsize=24)

def make_cnv_histogram_info(patient_db, chromosome, genome_dict):
    cnvs = patient_db.cnvs[chromosome]

    bin_starts = range(1, genome_dict["6"].length, 1000000)
    bin_counts = {range(bin_start, bin_start + 1000000): 0 for bin_start in bin_starts[:-1]}
    bin_counts[range(bin_starts[-1], genome_dict["6"].length + 1)] = 0

    for cnv in cnvs:
        for bin_range in bin_counts:
            if overlap(cnv.range, bin_range):
                bin_counts[bin_range] += 1
    return bin_counts

def plot_cnv_histogram(hist_info):
    bins = [(x.start - 1)/1000000 for x in hist_info]
    heights = list(hist_info.values())

    # plt.bar(y_pos, performance, align='center', alpha=0.5)
    plt.bar(bins, heights, align="center", alpha=0.5, color="seagreen")
    plt.xticks(fontsize=30, fontname="Tahoma")
    plt.yticks(fontsize=30, fontname="Tahoma")
    plt.xlabel("Megabase", fontsize=36, fontname="Tahoma")
    plt.ylabel("Count", fontsize=36, fontname="Tahoma")
    plt.title("CNV coverage per chromosome 6 megabase", fontsize=36, fontname="Tahoma")

    plt.show()


# %% Main
def main(geno_file, pheno_file):
    """Run main."""
    genotypes = DataManager.read_data(geno_file)
    phenotypes = DataManager.read_data(pheno_file)
    patients = DataManager.make_patients(genotypes, phenotypes)
    DataManager.print_summary_counts(patients)


def test():
    """Test run."""
    # Read genome dictionary.
    print("Reading reference genome dicionary.")
    genomedict = GenomeDict("C:/Users/tyler/Documents/Chr6/human_g1k_v37_phiX.dict")

    # Read geneset.
    print("Reading gene set...")
    geneset = GeneSet("C:/Users/tyler/Documents/Chr6/hg19.ensGene.gtf.gz")

    # Read HPO ontology.
    print("Loading Human Phenotype Ontology...")
    ontology = Ontology("C:/Users/tyler/Documents/Chr6/HPO/hp2.obo")

    # Read patient genotypes.
    print("Reading patient genotype data...")
    genotypes = DataManager.read_data("C:/Users/tyler/Documents/Chr6/genotypes.csv")
    genotypes = DataManager.fix_genotype_data(genotypes)
    # genotypes = trim_chromosome_names(genotypes)

    # Read patient phenotypes.
    print("Reading patient phenotype data...")
    phenotypes = DataManager.read_data("C:/Users/tyler/Documents/Chr6/phenotypes.csv")

    # Read patient HPO terms.
    print("Reading patient HPO data...")
    hpos = DataManager.read_data("C:/Users/tyler/Documents/Chr6/c6_research_patients_2020-10-28_11_27_04.csv")
    hpos = DataManager.fix_patient_hpos(hpos)

    # Make HI Gene objects.
    print("Reading HI gene information...")
    mg = mygene.MyGeneInfo()
    hi_genes = DataManager.read_HI_genes("C:/Users/tyler/Documents/Chr6/HI_chr6.bed")
    symbol_lookup = DataManager.symbol_lookup_multi(mg, list(hi_genes.keys()))
    hi_genes = DataManager.make_HI_objects(hi_genes, geneset, symbol_lookup)
    hi_genes = {x: y for x, y in hi_genes.items() if y.refined}

    # Build patient objects.
    print("Building patient objects...")
    patients = DataManager.make_patients(genotypes, phenotypes, geneset,
                                         hpos, ontology)
    patients.update(hi_genes)
    patients = PatientDatabase(patients)
    # patients = PatientDatabase({patient.id: patient for patient in list(patients.values())[:10]})
    # DataManager.print_summary_counts(patients)
    comparison = ComparisonTable(patients)
    print("Done.")
    return (
        genotypes,
        phenotypes,
        patients,
        comparison,
        geneset,
        ontology,
        hi_genes,
        genomedict
        )


if __name__ == "__main__":
    # import sys
    # main(sys.argv[1], sys.argv[2])
    # genotypes, phenotypes, patients, geneset = test()
    my_genotypes, my_phenotypes, my_patients, my_comparison, my_geneset, my_ontology, my_hi_genes, my_genomedict = test()
