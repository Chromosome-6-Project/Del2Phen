
import csv
import importlib.resources as pkg_resources

from pronto import Definition, Ontology, TermSet
import yaml

from chr6_project import resources


def read_hpo_ontology():
    with pkg_resources.path(resources, "hpo.obo") as ont_file:
        ontology = Ontology(ont_file)
    return ontology


def read_custom_termset_yaml(path):
    with open(path) as infile:
        termset_info = yaml.safe_load(infile)
    return termset_info


def add_custom_terms(termset_info, ontology):
    for term_id, term_values in termset_info["custom"].items():
        ontology.create_term(term_id)
        ontology[term_id].name = term_values["name"]
        # ontology[term_id].members = term_values["members"]
        ontology[term_id].comment = "Custom term"
        members = ",".join(term_values["members"])
        ontology[term_id].definition = Definition(f"{term_values['rule']}|{members}")


def make_termset(termset_info, ontology):
    termset = TermSet()
    for term_id in termset_info["single"] + list(termset_info["custom"]):
        termset.add(ontology[term_id])
    return termset


# def make_termset_and_add_to_ontology(termset_info, ontology):
#     termset = TermSet()
#     for term_id in termset_info["single"]:
#         termset.add(ontology[term_id])
#     for term_id, term_values in termset_info["custom"].items():
#         ontology.create_term(term_id)
#         ontology[term_id].name = term_values["name"]
#         ontology[term_id].members = term_values["members"]
#         ontology[term_id].comment = "Custom term"
#         ontology[term_id].description = term_values["rule"]
#         termset.add(ontology[term_id])
#     return termset


def read_patient_hpo_data(data_file, ontology):
    with open(data_file, encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter=",", quotechar='"')
        data = {line["label"]: {ontology[hp_id.replace("_", ":")]: response
                                for hp_id, response in list(line.items())[1:]}
                for line in reader}
    return data


def assign_custom_response(custom_term, term_response_dict, ontology):
    # responses = {term_response_dict[ontology[member_id]]
    #              for member_id in custom_term.members}
    # responses = {term_response_dict[ontology[member_id]]
    #              for member_id in custom_term.definition.split("|")[1].split(",")}
    rule, members = custom_term.definition.split("|")
    members = members.split(",")
    responses = {term_response_dict[ontology[member_id]] for member_id in members}
    if rule == "any":
        if "T" in responses:
            return "T"
        if responses == {"F"}:
            return "F"
    elif rule == "none":
        if responses == {"F"}:
            return "T"
        if responses == {"T"}:
            return "F"
    return "NA"


def assign_custom_responses(termset, patient_term_dict, ontology):
    terms = [term for term in termset if term.comment == "Custom term"]
    for patient, term_dict in patient_term_dict.items():
        for term in terms:
            term_dict[term] = assign_custom_response(term, term_dict, ontology)


# def add_custom_terms(ontology):
#     ontology.create_term("HP:0011470:+1")
#     ontology["HP:0011470:+1"].name = "Tube feeding in infancy"
#     ontology.create_term("NEURODEV_NORMAL")


# def read_hpo_list(path, ontology):
#     with open(path) as infile:
#         terms = infile.readlines()
#     terms = [ontology[term.strip().split("\t")[1]] for term in terms
#              if not term.startswith("#")]
#     return terms

def expand_patient_hpo_terms(patient_term_dict):
    expanded_term_dict = patient_term_dict
    expanded = {parent for term, response in patient_term_dict.items()
                for parent in term.superclasses(with_self=False)
                if response == "T"}
    for hpo in expanded:
        expanded_term_dict[hpo] = "T"
    return expanded_term_dict


def expand_all_patients_hpo_terms(patient_hpo_data):
    expanded_hpo_data = {}
    for patient, patient_term_dict in patient_hpo_data.items():
        expanded_hpo_data[patient] = expand_patient_hpo_terms(patient_term_dict)
    return expanded_hpo_data


def assign_false_subclasses(patient_term_dict):
    auto_falses = patient_term_dict
    definitives = {"T", "F"}
    non_answers = [term for term, response in patient_term_dict.items()
                   if response not in definitives]
    for term in non_answers:
        parent_responses = {patient_term_dict[parent]
                            for parent in term.superclasses(1, False)
                            if parent in patient_term_dict}
        if parent_responses == {"F"}:
            auto_falses[term] = "F"
    return auto_falses


def make_c6_hpo(patient_hpo_file, custom_termset_yaml=None, recursive_expansion=False):
    ontology = read_hpo_ontology()
    patient_hpo_data = read_patient_hpo_data(patient_hpo_file, ontology)
    if recursive_expansion:
        patient_hpo_data = expand_all_patients_hpo_terms(patient_hpo_data)
    if custom_termset_yaml is None:
        return ontology, patient_hpo_data, None
    termset_info = read_custom_termset_yaml(custom_termset_yaml)
    add_custom_terms(termset_info, ontology)
    termset = make_termset(termset_info, ontology)
    assign_custom_responses(termset, patient_hpo_data, ontology)
    return ontology, patient_hpo_data, termset


# if __name__ == "__main__":
#     import sys
#     ontology, patient_hpo_data, termset = make_c6_hpo(sys.argv[1], sys.argv[2])