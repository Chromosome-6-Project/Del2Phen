
import importlib.resources as pkg_resources

from pronto import Ontology

from chr6_project.analysis.chr6 import ComparisonTable, DataManager
from chr6_project.analysis.data_objects import PatientDatabase
from chr6_project.analysis.gene_set import read_geneset_gtf
from chr6_project import resources


def analyze(genotypes, phenotypes, patient_hpo, geneset_gtf, drop_list_file,
            expand_hpos=False):
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
    print("Loading Human Phenotype Ontology...")
    with pkg_resources.path(resources, "hpo.obo") as ont_file:
        ontology = Ontology(ont_file)

    # Read patient genotypes.
    print("Reading patient genotype resources...")
    genotypes = DataManager.read_data(genotypes)
    genotypes = DataManager.fix_genotype_data(genotypes)
    # genotypes = trim_chromosome_names(genotypes)

    # Read patient phenotypes.
    print("Reading patient phenotype resources...")
    phenotypes = DataManager.read_data(phenotypes)

    # Read patient HPO terms.
    print("Reading patient HPO resources...")
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
        # genotypes,
        # phenotypes,
        # patients,
        comparison,
        geneset,
        ontology,
        # hi_genes,
        # genomedict
        )
