
import importlib.resources as pkg_resources

from pronto import Ontology

from chr6_project import resources


def read_hpo_ontology():
    with pkg_resources.path(resources, "hpo.obo") as ont_file:
        ontology = Ontology(ont_file)
    return ontology


def add_custom_terms(ontology):
    ontology.create_term("TUBE_FEED")
    ontology.create_term("NEURODEV_NORMAL")


def make_c6_hpo():
    ontology = read_hpo_ontology()
    add_custom_terms(ontology)
    return ontology
