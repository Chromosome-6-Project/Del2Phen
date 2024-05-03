
from collections import namedtuple
import importlib.resources as pkg_resources

import numpy as np
import pandas as pd
from pronto import Definition, Ontology, TermSet
import yaml

from del2phen import resources


def read_hpo_ontology():
    with pkg_resources.path(resources, "hpo.obo") as ont_file:
        ontology = Ontology(ont_file)
    return ontology


def get_default_termset_yaml_path():
    yaml_path = list(pkg_resources.path(resources, "default_termset.yaml").gen)[0]
    return yaml_path


def get_default_custom_phenotypes_path():
    path = list(pkg_resources.path(resources, "custom_phenotypes.tsv").gen)[0]
    return path


def read_termset_yaml(path):
    with open(path) as infile:
        termset_info = yaml.safe_load(infile)
    return termset_info


Phenoterm = namedtuple("Phenoterm", ["name", "term_id", "c6_field"])


def parse_phenotypes(phenotypes: pd.DataFrame, custom_phenotypes_file=None,
                     phenotype_termset_yaml=None, expand_hpo_terms=False):
    ontology = read_hpo_ontology()
    if custom_phenotypes_file is not None:
        custom_terms = pd.read_csv(custom_phenotypes_file, sep="\t")
        add_custom_phenotypes(ontology, custom_terms)
    records = phenotypes.set_index("id").to_dict(orient="index")
    records = {pid: {ontology[term]: response for term, response in responses.items()}
               for pid, responses in records.items()}
    if expand_hpo_terms:
        records = expand_all_patients_hpo_terms(records)
    if phenotype_termset_yaml is None:
        return records, ontology, None
    termset_info = read_termset_yaml(phenotype_termset_yaml)
    add_custom_terms(termset_info, ontology)
    termset = make_termset(termset_info, ontology)
    assign_custom_responses(termset, records, ontology)
    return records, ontology, termset


def add_custom_phenotypes(ontology, custom_phenotypes: pd.DataFrame):
    for term in custom_phenotypes.itertuples():
        ontology.create_term(term.term_id)
        ontology[term.term_id].name = term.label
        ontology[term.term_id].comment = "Custom term."


def expand_all_patients_hpo_terms(patient_hpo_data):
    expanded_hpo_data = {}
    for patient, patient_term_dict in patient_hpo_data.items():
        expanded_hpo_data[patient] = expand_patient_hpo_terms(patient_term_dict)
    return expanded_hpo_data


def expand_patient_hpo_terms(patient_term_dict):
    expanded_term_dict = patient_term_dict
    expanded = {parent for term, response in patient_term_dict.items()
                for parent in term.superclasses(with_self=False)
                if response == True}
    for hpo in expanded:
        expanded_term_dict[hpo] = True
    return expanded_term_dict


def add_custom_terms(termset_info, ontology):
    for term_id, term_values in termset_info["custom"].items():
        ontology.create_term(term_id)
        ontology[term_id].name = term_values["name"]
        # ontology[term_id].members = term_values["members"]
        ontology[term_id].comment = "Rule term."
        members = ",".join(term_values["members"])
        ontology[term_id].definition = Definition(f"{term_values['rule']}|{members}")


def make_termset(termset_info, ontology):
    termset = TermSet()
    for term_id in termset_info["single"] + list(termset_info["custom"]):
        termset.add(ontology[term_id])
    return termset


def assign_custom_responses(termset, patient_term_dict, ontology):
    terms = [term for term in termset if term.comment == "Rule term"]
    for patient, term_dict in patient_term_dict.items():
        for term in terms:
            term_dict[term] = assign_custom_response(term, term_dict, ontology)


def assign_custom_response(custom_term, term_response_dict, ontology):
    rule, members = custom_term.definition.split("|")
    members = members.split(",")
    responses = {term_response_dict[ontology[member_id]] for member_id in members}
    if rule == "any":
        if True in responses:
            return True
        if responses == {False}:
            return False
    elif rule == "none":
        if responses == {False}:
            return True
        if True in responses:
            return False
    return np.NaN
