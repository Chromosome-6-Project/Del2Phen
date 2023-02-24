"""The Chromosome 6 Project - Main Analysis Methods.

@author: T.D. Medina
"""

import csv

from chr6_project.analysis.patient_comparison import ComparisonTable
from chr6_project.analysis.data_objects import PatientDatabase, Patient
from chr6_project.analysis.gene_set import make_geneset, read_geneset_gtf
from chr6_project.molgenis_import.import_data import import_chr6_data
import chr6_project.analysis.hpo as c6_hpo


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
        start_report = "Start_positie_in_report"
        start_convert = "Start_positie_in_Hg19"
        stop_report = "Stop_positie_in_report"
        stop_convert = "Stop_positie_in_Hg19"

        for entry in fixed:
            build = entry["Human_genome_build_version"]
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
                    print(f"Error while parsing:\n{err}")
        return fixed

    @staticmethod
    def trim_chromosome_names(data):
        """Remove 'chr' from contig names."""
        trimmed = data
        for entry in trimmed:
            if not entry["Chromosoom"]:
                continue
            if entry["Chromosoom"].lower().startswith("chr"):
                entry["Chromosoom"] = entry["Chromosoom"].lstrip("chrCHR")
        return trimmed

    @staticmethod
    def fix_patient_hpos(data):
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
                      ontology=None):
        """Construct Patient objects from data."""
        if isinstance(phenotypes, dict):
            phenotypes = list(phenotypes.values())
        subjects = ({record["owner"] for record in genotypes}
                    | {record["id"] for record in phenotypes})
        patients = {subject: Patient(subject) for subject in subjects}
        for genotype in genotypes:
            patients[genotype["owner"]].genotypes.append(genotype)
        for phenotype in phenotypes:
            patients[phenotype["id"]].phenotypes.append(phenotype)
        if geneset:
            for patient in patients.values():
                patient.cnvs = patient.extract_cnvs()
                patient.identify_gene_overlaps(geneset)
        if ontology and hpos:
            for patient in patients.values():
                if patient.id not in hpos:
                    continue
                patient.hpo = hpos[patient.id]
                # patient.hpo = {ontology[hpo_id]: response
                #                for hpo_id, response in hpos[patient.id].items()}

                # patient_hpos = []
                # for hpo_id, value in hpos[patient.id].items():
                #     if hpo_id == "label":
                #         continue
                #     if value == "T":
                #         patient_hpos.append(ontology[hpo_id])
                # patient.hpo = patient_hpos

                # if expand_hpos:
                #     patient.expand_hpo_terms(ontology)
        for patient in patients.values():
            if patient.genotypes:
                if "origin info" in patient.genotypes[0]:
                    patient.origin = patient.genotypes[0]["origin info"]
                else:
                    patient.origin = None
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
    def filter_patients(patient_dict, drop_list=None, remove_ungenotyped=True,
                        remove_unphenotyped=True, remove_unsubmitted=True):
        if drop_list is None:
            drop_list = set()
        patients_filtered = {
            name: patient for name, patient in patient_dict.items()
            if all([
                name not in drop_list,
                not remove_ungenotyped or patient.genotypes,
                not remove_unphenotyped or patient.phenotypes,
                not remove_unsubmitted or patient.check_if_submitted()
                ])
            }
        return patients_filtered


def analyze(genotypes, phenotypes, patient_hpo, geneset_gtf, drop_list_file,
            hpo_termset_yaml=None, expand_hpos=False, remove_ungenotyped=True,
            remove_unphenotyped=True, remove_unsubmitted=True):
    # Read geneset.
    print("Loading gene set...")
    # Build GeneSet from source GTF file (gzipped) (slow, large file):
    # geneset = GeneSet("C:/Users/Ty/Documents/Chr6/hg19.ensGene.gtf.gz")

    # Build chr6 GeneSet only from source GTF file (gzipped) (slow, large file)
    # and add pLI and HI info automatically from default sources:
    geneset = read_geneset_gtf(geneset_gtf)

    # Or, load pre-made GeneSet from pickle (faster, large file):
    # with open("GeneSets/hg19.ensGene.pkl", "rb") as infile:
    #     geneset = pickle.load(infile)

    # Or, load pre-made GeneSet from bz2 pickle (less fast, small file):
    # with bz2.BZ2File("GeneSets/hg19.ensGene.pkl.bz2", "rb") as infile:
    #     geneset = cPickle.load(infile)

    # Read HPO ontology.
    # print("Loading Human Phenotype Ontology...")
    # ontology = c6_hpo.make_c6_hpo()

    # Read patient genotypes.
    print("Reading patient genotype resources...")
    genotypes = DataManager.read_data(genotypes)
    genotypes = DataManager.fix_genotype_data(genotypes)
    # genotypes = trim_chromosome_names(genotypes)

    # Read patient phenotypes.
    print("Reading patient phenotype resources...")
    phenotypes = DataManager.read_data(phenotypes)

    # Read patient HPO terms.
    print("Loading patient HPO resources...")
    ontology, hpos, termset = c6_hpo.make_c6_hpo(patient_hpo, hpo_termset_yaml,
                                                 expand_hpos)
    # hpos = DataManager.read_data(patient_hpo)
    # hpos = DataManager.fix_patient_hpos(hpos)
    # hpos = DataManager.add_custom_hpos(hpos)

    # Build patient objects.
    # !!!: This is where you can choose whether to expand HPO terms.
    print("Building patient objects...")
    patients = DataManager.make_patients(genotypes, phenotypes, geneset,
                                         hpos, ontology)

    print("Filtering patients...")
    drop_list = DataManager.read_patient_drop_list(drop_list_file)
    patients = DataManager.filter_patients(patients, drop_list, remove_ungenotyped,
                                           remove_unphenotyped, remove_unsubmitted)

    patients = PatientDatabase(patients)

    print("Running comparisons...")
    comparison = ComparisonTable(patients)

    print("Done.")
    return (
        # genotypes,
        # phenotypes,
        # patients,
        comparison,
        geneset,
        ontology,
        termset,
        # hi_genes,
        # genomedict
        )


def analyze_online(username, password, drop_list_file=None, hpo_termset_yaml=None,
                   expand_hpos=False, remove_ungenotyped=True, remove_unphenotyped=True,
                   remove_unsubmitted=True):
    print("Loading geneset...")
    geneset = make_geneset()

    print("Downloading patient data from Molgenis...")
    genotypes, phenotypes, hpos = import_chr6_data(username, password)

    print("Organizing patient HPO terms...")
    if hpo_termset_yaml is None:
        hpo_termset_yaml = c6_hpo.get_default_termset_yaml_path()
    ontology, hpos, termset = c6_hpo.make_c6_hpo_online(hpos, hpo_termset_yaml, expand_hpos)

    print("Building patient objects...")
    patients = DataManager.make_patients(genotypes, phenotypes, geneset, hpos, ontology)

    if drop_list_file is not None:
        print("Filtering patients...")
        drop_list = DataManager.read_patient_drop_list(drop_list_file)
        patients = DataManager.filter_patients(patients, drop_list, remove_ungenotyped,
                                               remove_unphenotyped, remove_unsubmitted)

    print("Building patient database...")
    patients = PatientDatabase(patients)

    print("Filtering patients with duplications...")
    patients = patients.remove_patients_by_cnv_type("Duplication x3")

    print("Comparing patients...")
    comparison = ComparisonTable(patients)

    print("Done")
    return comparison, geneset, ontology, termset
