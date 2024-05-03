"""The Chromosome 6 Project - Main Analysis Workflow.

@author: T.D. Medina
"""

from typing import List, Optional, Union

import pandas as pd

from del2phen.analysis.patient_comparison import ComparisonTable
from del2phen.analysis.data_objects import PatientDatabase, Patient, CNV
from del2phen.analysis.gene_set import GeneSet, make_geneset, read_geneset_gtf
import del2phen.analysis.hpo as c6_hpo


def make_patients(genotypes: pd.DataFrame, phenotypes: dict, geneset: GeneSet = None):
    """Construct Patient objects from data."""
    ids = set(genotypes.id) | set(phenotypes.keys())
    patients = {pid: Patient(pid) for pid in ids}
    for pid, phenotype in phenotypes.items():
        patients[pid].phenotypes = phenotype
    for genotype in genotypes.to_dict(orient="records"):
        if pd.isna(genotype["start"]) or pd.isna(genotype["stop"]):
            continue
        genotype["start"] = int(genotype["start"])
        genotype["stop"] = int(genotype["stop"])
        cnv = CNV(**genotype)
        if geneset is not None:
            cnv.set_genes(geneset)
        patients[genotype["id"]].cnvs.append(cnv)
    return patients

    # TODO: Move this to PatientDatabase and update.
    # @staticmethod
    # def make_ucsc_browser_tracks(patients, out, filter_chrom=None):
    #     """Write UCSC BED file of patient CNVs."""
    #     writer = ["browser hide all position chr6\n"]
    #     patient_list = sorted(
    #         [patient for patient in patients.values() if patient.cnvs],
    #         key=lambda x: x.cnvs[0].range.start
    #         )
    #     for i, patient in enumerate(patient_list):
    #         patient_writer = [f"track name=Patient_{i} visibility=2\n"
    #                           f"#chrom\tchromStart\tchromEnd\tname\n"]
    #         for cnv in patient.cnvs:
    #             patient_writer.append(
    #                 f"chr{cnv.chromosome}\t{cnv.range.start}\t"
    #                 f"{cnv.range.stop - 1}\t{str(cnv.copy_number).replace(' ', '')}\n")
    #         writer += patient_writer
    #     with open(out, "w") as outfile:
    #         outfile.writelines(writer)


def read_patient_drop_list(drop_list_file):
    with open(drop_list_file) as infile:
        drop_list = infile.readlines()
    drop_list = [name.strip() for name in drop_list]
    drop_list = {name for name in drop_list if not name.startswith("#")}
    return drop_list


def filter_patients(patient_dict, drop_list=None,
                    chromosomes: Optional[Union[str, List[str]]] = None,
                    copy_numbers: Optional[Union[int, List[int]]] = None,
                    keep_uncompared=False, keep_ungenotyped=False,
                    keep_unphenotyped=False):
    if drop_list is None:
        drop_list = set()
    patients_filtered = {
        name: patient for name, patient in patient_dict.items()
        if all([
            name not in drop_list,
            keep_uncompared or patient.filter_cnvs(chromosomes, copy_numbers),
            keep_ungenotyped or patient.cnvs,
            keep_unphenotyped or patient.phenotypes,
            ])
        }
    return patients_filtered


def analyze(genotypes, phenotypes, drop_list_file=None, custom_phenotype_file=None,
            chromosomes: Optional[Union[str, List[str]]] = None,
            copy_numbers: Optional[Union[int, List[int]]] = None,
            phenotype_termset_yaml=None, expand_hpo_terms=False, keep_ungenotyped=False,
            keep_unphenotyped=False, keep_uncompared=False,
            pLI_threshold=0.9, HI_threshold=10, phaplo_threshold=0.86, mode="confirm",
            **kwargs):

    # TODO: Add verbosity argument.
    # TODO: Decide on how to handle the geneset.
    # Read geneset.
    print("Loading gene set...")
    # Build GeneSet from source GTF file (gzipped) (slow, large file):
    # geneset = GeneSet("C:/Users/Ty/Documents/Chr6/hg19.ensGene.gtf.gz")

    # Build chr6 GeneSet only from source GTF file (gzipped) (slow, large file)
    # and add pLI and HI info automatically from default sources:
    geneset = make_geneset()

    # Or, load pre-made GeneSet from pickle (faster, large file):
    # with open("GeneSets/hg19.ensGene.pkl", "rb") as infile:
    #     geneset = pickle.load(infile)

    # Or, load pre-made GeneSet from bz2 pickle (less fast, small file):
    # with bz2.BZ2File("GeneSets/hg19.ensGene.pkl.bz2", "rb") as infile:
    #     geneset = cPickle.load(infile)

    print("Reading patient genotypes...")
    genotypes = pd.read_csv(genotypes, sep="\t")

    print("Reading patient phenotypes...")
    phenotypes = pd.read_csv(phenotypes, sep="\t",
                             true_values=["T", "t", 1],
                             false_values=["F", "f", 0],
                             na_values=[
            "Unsure"])

    print("Parsing phenotypes...")
    phenotypes, ontology, termset = c6_hpo.parse_phenotypes(phenotypes,
                                                            custom_phenotype_file,
                                                            phenotype_termset_yaml,
                                                            expand_hpo_terms)

    print("Building patient objects...")
    patients = make_patients(genotypes, phenotypes, geneset)

    print("Filtering patients...")
    drop_list = None if drop_list_file is None else read_patient_drop_list(drop_list_file)
    patients = filter_patients(patients, drop_list,
                               chromosomes, copy_numbers,
                               keep_uncompared, keep_ungenotyped,
                               keep_unphenotyped)
    patients = PatientDatabase(patients)

    print("Running comparisons...")
    comparison = ComparisonTable(patients,
                                 chromosomes=chromosomes,
                                 copy_numbers=copy_numbers,
                                 pLI_threshold=pLI_threshold,
                                 HI_threshold=HI_threshold,
                                 phaplo_threshold=phaplo_threshold,
                                 mode=mode)

    print("Done.")
    return (
        comparison,
        geneset,
        ontology,
        termset,
        )
