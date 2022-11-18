
import molgenis.client


def import_chr6_data(username, password, remove_unsubmitted=True):
    session = molgenis.client.Session("https://chromosome6-research.gcc.rug.nl/")
    session.login(username, password)

    submit_status = session.get("c6_questionnaire", attributes="owner,status")
    submitted = {entry["owner"] for entry in submit_status if entry["status"] == "SUBMITTED"}

    genotypes = session.get("c6_array")
    genotypes = set_href_attribute_names(genotypes, "label")
    fix_contig_names(genotypes)

    phenos = session.get("c6_research_patients")
    phenos = set_href_attribute_names(phenos, "id")

    if remove_unsubmitted:
        genotypes = [entry for entry in genotypes if entry["owner"] in submitted]
        phenos = [entry for entry in phenos if entry["id"] in submitted]

    phenos, hpos = split_hpo_from_phenos(phenos)

    for entry in phenos.values():
        entry["submit_status"] = set()
    for entry in submit_status:
        if entry["owner"] in phenos:
            phenos[entry["owner"]]["submit_status"].add(entry["status"])

    return genotypes, phenos, hpos


def split_hpo_from_phenos(phenotypes):
    pheno_data = {}
    hpo_data = {}
    for entry in phenotypes:
        new_pheno_entry = dict()
        new_hpo_entry = dict()
        for k, v in entry.items():
            if k.startswith("HP_"):
                new_hpo_entry[k.replace("_", ":")] = v
            else:
                new_pheno_entry[k] = v
        pheno_data[entry["id"]] = new_pheno_entry
        hpo_data[entry["id"]] = new_hpo_entry
    return pheno_data, hpo_data


def set_href_attribute_names(data, label_name):
    new_data = []
    for entry in data:
        new_entry = {}
        for k, v in entry.items():
            if isinstance(v, dict):
                new_entry[k] = v[label_name]
            else:
                new_entry[k] = v
        new_data.append(new_entry)
    return new_data


def fix_contig_names(genotype_data):
    for entry in genotype_data:
        entry["Chromosoom"] = entry["Chromosoom"].lstrip("chrCHR")
