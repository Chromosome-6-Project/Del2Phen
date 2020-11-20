
import csv
import gzip
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
from pronto import Ontology
import mygene

REFERENCE_CHR = [str(i) for i in range(1, 23)] + ["X", "Y"]


def overlap(range1: range, range2: range):
    if range1.start <= range2.stop and range2.start <= range1.stop:
        return True
    return False


class PatientIntersect:
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
        string = (f"PatientIntersect(patients={self.patients}, "
                  f"gene_count={self.gene_count}, hpo_count={self.hpo_count})")
        return string

    def __str__(self):
        string = (f"{self.patients[0]} vs. {self.patients[1]}:\n"
                  f"  Genes: {self.gene_count}\n"
                  f"  HPO terms: {self.hpo_count}")
        return string

    def gene_info(self):
        return self.genes, self.gene_count, self.gene_similarity

    def hpo_info(self):
        return self.hpos, self.hpo_count, self.hpo_similarity

class ComparisonTable:
    def __init__(self, raw, patients=None, pack=True):
        self.raw = raw
        self.index = self.make_index()
        self.array = self.make_array(pack)
        self.patients = patients
        self.size = len(self.index)

        self.__iteri__ = 0
        self.__iterj__ = 0

    def __iter__(self):
        self.__iteri__ = 0
        self.__iterj__ = 0
        return self

    def __next__(self):
        if self.__iteri__ == self.size:
            raise StopIteration
        result = self.array[self.__iteri__][self.__iterj__]
        self.__iterj__ +=1
        if self.__iterj__ == self.size:
            self.__iteri__ += 1
            self.__iterj__ = self.__iteri__
        return result

    def tabulate(self):
        table = {}
        # self.

    def make_array(self, pack=True):
        array = []
        for pid in self.index:
            values = []
            for p2 in self.index:
                # if pid == p2:
                    # values.append(({}, 1, 1, {}, 1, 1))
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
            array = np.array(array, dtype=[("genes", object), ("genecount", int), ("genesim", float), ("hpos", object), ("hpocount", int), ("hposim", float)])
        return array

    def make_index(self):
        return {j: i for i, j in enumerate(self.raw)}

    def lookup(self, pid1, pid2=None):
        if not pid2:
            return self.array[self.index[pid1]][self.index[pid1]]
        if pid1 not in self.index:
            print("ID not found.")
            return
        return self.array[self.index[pid1]][self.index[pid2]]


# def compare_genes(patient1, patient2):
#     patient1_set = set([gene for genelist in patient1.affected_genes.values()
#                         for gene in genelist])
#     patient2_set = set([gene for genelist in patient2.affected_genes.values()
#                         for gene in genelist])
#     if not patient1_set or not patient2_set:
#         return "NA"
#     total_set = patient1_set | patient2_set
#     return len(patient1_set & patient2_set) / len(total_set)


# def make_patient_geneset(patient, ids_only=False):
#     if ids_only:
#         patient_genes = set([gene for genelist in patient.affected_gene_ids.values()
#                             for gene in genelist])
#     else:
#         patient_genes = set([gene for genelist in patient.affected_genes.values()
#                             for gene in genelist])
#     return patient_genes


# XXX: Just use compound compare instead. This still works though, in general.
# def build_gene_comparison_table(patients, skip_empty=True,
#                                 transcripts_only=True, mode="full"):
#     gene_sets = {patient.id: patient.all_genes(transcripts_only)
#                  for patient in patients.values()}
#     if skip_empty:
#         gene_sets = {k: v for k, v in gene_sets.items() if v}
#     patient_ids = sorted(list(gene_sets.keys()))
#     comparisons = {}
#     while patient_ids:
#         patient1 = patient_ids.pop()
#         patient1_set = gene_sets[patient1]
#         patient_compare = {}
#         for patient2 in patient_ids:
#             patient2_set = gene_sets[patient2]
#             if not patient1_set or not patient2_set:
#                 patient_compare[patient2] = None
#                 continue
#             if mode == "count":
#                 patient_compare[patient2] = len(patient1_set & patient2_set)
#             elif mode == "similarity":
#                 patient_compare[patient2] = (len(patient1_set & patient2_set)
#                                              / len(patient1_set | patient2_set))
#             elif mode == "full":
#                 genes = patient1_set & patient2_set
#                 count = len(genes)
#                 similarity = count / len(patient1_set | patient2_set)
#                 patient_compare[patient2] = (genes, count, similarity)
#         comparisons[patient1] = patient_compare
#     return comparisons


# XXX: As above, just use compound compare.
# def build_hpo_comparison_table(patients, skip_empty=True, mode="full"):
#     hpo_sets = {patient.id: set(patient.hpo) for patient in patients.values()}
#     if skip_empty:
#         hpo_sets = {k: v for k, v in hpo_sets.items() if v}
#     patient_ids = sorted(list(hpo_sets.keys()))
#     comparisons = {}
#     while patient_ids:
#         patient1 = patient_ids.pop()
#         patient1_set = hpo_sets[patient1]
#         patient_compare = {}
#         for patient2 in patient_ids:
#             patient2_set = hpo_sets[patient2]
#             if not patient1_set or not patient2_set:
#                 patient_compare[patient2] = None
#                 continue
#             if mode == "count":
#                 patient_compare[patient2] = len(patient1_set & patient2_set)
#             elif mode == "similarity":
#                 patient_compare[patient2] = (len(patient1_set & patient2_set)
#                                              / len(patient1_set | patient2_set))
#             elif mode == "full":
#                 hpos = patient1_set & patient2_set
#                 count = len(hpos)
#                 similarity = count / len(patient1_set | patient2_set)
#                 patient_compare[patient2] = (hpos, count, similarity)
#         comparisons[patient1] = patient_compare
#     return comparisons


def compound_comparison(patients): #, gene_comp, hpo_comp):
    patient_ids = sorted([patient.id for patient in patients.values()])
    genes = {patient.id: patient.all_genes() for patient in patients.values()}
    hpos = {patient.id: set(patient.hpo) for patient in patients.values()}

    comparisons = {}
    while patient_ids:
        patient1 = patient_ids.pop()
        genes1 = genes[patient1]
        hpo1 = hpos[patient1]

        patient_comparison = {}
        patient_comparison[patient1] = (genes1, len(genes1), 1,
                                        hpo1, len(hpo1), 1)

        for patient2 in patient_ids:
            genes2 = genes[patient2]
            hpo2 = hpos[patient2]

            gene_union = genes1 & genes2
            gene_len = len(gene_union)
            if gene_len == 0:
                gene_similarity = 0.0
            else:
                gene_similarity = gene_len / len(genes1 | genes2)

            hpo_union = hpo1 & hpo2
            hpo_len = len(hpo_union)
            if hpo_len == 0:
                hpo_similarity = 0.0
            else:
                hpo_similarity = hpo_len / len(hpo1 | hpo2)

            comparison = (gene_union, gene_len, gene_similarity,
                          hpo_union, hpo_len, hpo_similarity)
            patient_comparison[patient2] = comparison

        comparisons[patient1] = patient_comparison
    return comparisons


def write_comparison_table(table, patients, out, self_match="size"):
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
    return


# %% Network Stuff

Node = namedtuple("Node", ["id", "label", "group", "color", "value", "title"])
Edge = namedtuple("Edge", ["edge_id", "id1", "id2", "width", "color", "title", "sim"])

def build_network_nodes2(comparison_table, patients):
    colors = {"Literature case report": "blue",
              "Parental uploaded array report": "pink",
              None: "black",
              "HI Gene": "yellow"}

    nodes = {}
    for patient in comparison_table.index:
        lookup = comparison_table.lookup(patient)
        group = patients[patient].origin
        color = colors[group]
        value = lookup.gene_count
        title = (f"<p>{patient}<br>"
                 f"Group: {group}<br>"
                 f"Affected transcripts: {value}<br>"
                 f"HPO terms: {lookup.hpo_count}</p>")
        node = Node(patient, patient, group, color, value, title)
        nodes[patient] = node
    return nodes


def write_network_nodes2(nodes, out):
    writer = ["id,label,group,color,value,title\n"]
    for node in sorted(list(nodes.values()), key=lambda x: x[0]):
        writer.append(",".join([str(x) for x in node]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)
    return


def build_network_edges2(comparison_table):
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


def write_network_edges2(edges, out):
    writer = ["id,from,to,width,color,title,gene_sim\n"]
    for edge in edges:
        writer.append(",".join([str(x) for x in edge]) + "\n")
    with open(out, "w") as outfile:
        outfile.writelines(writer)
    return


# def build_network_nodes(table, patients):
#     colors = {"Literature case report": "blue",
#               "Parental uploaded array report": "pink",
#               None: "black"}

#     node_ids = set(table.keys())
#     nodes = {node_id: {} for node_id in node_ids}

#     for node in nodes:
#         nodes[node]["color"] = colors[patients[node].genotypes[0]["origin info"]]
#         nodes[node]["value"] = len({geneID for geneset in patients[node].affected_gene_ids.values()
#                                     for geneID in geneset})
#         # nodes[node]["x"] = (patients[node].cnvs[0].range.stop
#         #                     + patients[node].cnvs[0].range.start) / 2
#         # nodes[node]["y"] = 0
#     return nodes

# def build_network_edges(table, nodes):
#     edges = []
#     for id1, comparisons in table.items():
#         for id2, value in comparisons.items():
#             if not value:
#                 continue
#             edge_dict = {}
#             edge_dict["from"] = id1
#             edge_dict["to"] = id2
#             edge_dict["width"] = value
#             if nodes[id1]["color"] == nodes[id2]["color"]:
#                 edge_dict["color"] = nodes[id1]["color"]
#             else:
#                 edge_dict["color"] = "black"
#             edges.append(edge_dict)
#     return edges


# def write_network_nodes(nodes, out):
#     writer = ["id,color,value,x,y\n"]
#     for node_id, node in nodes.items():
#         node_line = (f"{node_id},{node['color']},{node['value']},"
#                      f"{node['x']},{node['y']}\n")
#         writer.append(node_line)
#     with open(out, "w") as outfile:
#         outfile.writelines(writer)
#     return


# def write_network_edges(edges, out):
#     writer = ["from,to,width,color\n"]
#     for edge in edges:
#         edge_line = f"{edge['from']},{edge['to']},{edge['width']},{edge['color']}\n"
#         writer.append(edge_line)
#     with open(out, "w") as outfile:
#         outfile.writelines(writer)


# def write_edge_list(table, out):
#     edge_writer = [["from", "to", "width", "color"]]
#     for id1, comps in table.items():
#         for id2, value in comps.items():
#             if not value:
#                 continue
#             edge_writer.append([id1, id2, str(value)])
#     edge_writer = [",".join(line) + "\n" for line in edge_writer]
#     with open(out, "w") as outfile:
#         outfile.writelines(edge_writer)
#     return


# def write_node_list(table, out):
#     colors = {"Literature case report": "blue",
#               "Parental uploaded array report": "pink",
#               None: "gray"}
#     node_writer = [["id", "value", "color"]]
#     nodes = set(table.keys())
#     for id1 in nodes:
#         color = colors[patients[id1].genotypes[0]["origin info"]]
#         node_writer.append(
#             [id1,
#              str(len({y for x in patients[id1].affected_gene_ids.values()
#                       for y in x})),
#              color]
#              )
#     node_writer = [",".join(line) + "\n" for line in node_writer]
#     with open(out, "w") as outfile:
#         outfile.writelines(node_writer)
#     return


def make_array(table):
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


def gene_comparison_heatmap(table):
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


# def plot_patient(patient):
#     plt.figure()
#     rect = mpl.patches.Rectangle((i-(w/2), band.start), w, band.stop - band.start, fc=f"{stain_colors[band.stain]}")
#     plt.gca().add_patch(rect)

# %% Main

class GeneAnnotation:
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

    def __repr__(self):
        string = ("GeneAnnotation("
                  f"gene_id={self.gene_id}, "
                  f"locus={self.chromosome}:{self.range.start}-{self.range.stop - 1}"
                  f"")
        return string

    def __hash__(self):
        return hash(str(self.__dict__))

    def is_transcript(self):
        return self.feature == "transcript"

class GeneSet:
    def __init__(self, file=None):
        self.path = file
        self.genes = []
        self.gene_ids = []

        if self.path:
            self.genes = self.read_genes(self.path)

    @staticmethod
    def read_genes(file):
        # data = [
        with gzip.open(file) as infile:
            data = infile.readlines()
            # i = 0
            # while i < 50:
            #     data.append(infile.readline())
            #     i += 1
        data = [line.decode().rstrip(";\n").split("\t")
                for line in data]
        # !!! Skipping non-autosomes.
        data = [line for line in data if line[0].lstrip("chr") in REFERENCE_CHR]
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
        if stop is None:
            stop = start
        query = range(start, stop + 1)
        results = []
        for gene in self.genes[chromosome]:
            if overlap(query, gene.range):
                results.append(gene)
        return results


class DataManager:
    def __init__(self):
        pass

    @classmethod
    def read_data(cls, data_file):
        """Import CSV data."""
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
        return

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
    def make_patients(genotypes, phenotypes, geneset=None, hpos=None, ontology=None):
        """Construct Patient objects from data."""
        subjects = (set([record["subject"] for record in genotypes])
                    | set([record["owner"] for record in phenotypes]))
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
        return

    @staticmethod
    def make_UCSC_browser_tracks(patients, out, filter_chrom=None):
        writer = ["browser hide all position chr6\n"]
        patient_list = sorted([patient for patient in patients.values()
                               if patient.cnvs], key=lambda x: x.cnvs[0].range.start)
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
        return

    @staticmethod
    def read_HI_genes(file):
        with open(file) as infile:
            infile.readline()
            data = infile.readlines()
        data = [x.lstrip("chr").rstrip("\n").split("\t") for x in data]
        data = {x[3].split("|")[0]: [x[0], int(x[1]), int(x[2]), float(x[4])] for x in data}
        return data

# !!!
    @classmethod
    def make_HI_objects(cls, hi_genes, geneset, gene_lookup):
        mg = mygene.MyGeneInfo()
        hi_gene_objs = {}
        for geneID, hi_gene in hi_genes.items():
            this_hi = HI_Gene(geneID)
            this_hi.genotypes.append(dict(zip(["Chromosome", "Start positie in Hg19", "Stop positie in Hg19", "imbalance"],
                                              hi_gene[:3] + ["HI Gene"])))
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

    @staticmethod
    def symbol_lookup(mygene_instance, gene_symbol):
        results = mygene_instance.query(gene_symbol, fields="ensembl.gene", species="human")
        if len(results["hits"]) > 1:
            print(f"Warning: Multiple hits for {gene_symbol}")
            return None
        if len(results["hits"]) == 0:
            print(f"Warning: No results found for {gene_symbol}")
        return results["hits"][0]["ensembl"]["gene"]


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


class Patient:
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
            cnvs.append(CNV(chrom, start, stop, change))
        cnvs = sorted(cnvs, key=lambda x: REFERENCE_CHR.index(x.chromosome))
        return cnvs

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


def main(geno_file, pheno_file):
    """Run main."""
    genotypes = DataManager.read_data(geno_file)
    phenotypes = DataManager.read_data(pheno_file)
    patients = DataManager.make_patients(genotypes, phenotypes)
    DataManager.print_summary_counts(patients)


def test():
    """Test run."""
    # Read geneset.
    geneset = GeneSet("C:/Users/tyler/Documents/Chr6/hg19.ensGene.gtf.gz")

    # Read HPO ontology.
    ontology = Ontology("C:/Users/tyler/Documents/Chr6/HPO/hp2.obo")

    # Read patient genotypes.
    genotypes = DataManager.read_data("C:/Users/tyler/Documents/Chr6/genotypes.csv")
    genotypes = DataManager.fix_genotype_data(genotypes)
    # genotypes = trim_chromosome_names(genotypes)

    # Read patient phenotypes.
    phenotypes = DataManager.read_data("C:/Users/tyler/Documents/Chr6/phenotypes.csv")

    # Read patient HPO terms.
    hpos = DataManager.read_data("C:/Users/tyler/Documents/Chr6/c6_research_patients_2020-10-28_11_27_04.csv")
    hpos = DataManager.fix_patient_hpos(hpos)

    # Build patient objects.
    patients = DataManager.make_patients(genotypes, phenotypes, geneset,
                                         hpos, ontology )
    DataManager.print_summary_counts(patients)
    return genotypes, phenotypes, patients, geneset, ontology


if __name__ == "__main__":
    # import sys
    # main(sys.argv[1], sys.argv[2])
    # genotypes, phenotypes, patients, geneset = test()
    my_genotypes, my_phenotypes, my_patients, my_geneset, my_ontology = test()
