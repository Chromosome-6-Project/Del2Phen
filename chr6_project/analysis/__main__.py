
import sys

from chr6_project.analysis.analyze import (
    # analyze,
    analyze_online,
    )
from chr6_project.analysis.hpo import get_default_termset_yaml_path


if __name__ == "__main__":
    # c6_dir = "/home/tyler/Documents/Chr6_docs/"
    # my_comparison, my_geneset, my_ontology, termset = analyze(
    #     genotypes=f"{c6_dir}/PatientData/2022-Oct-25/c6_array_2022-10-25_18_58_17.csv",
    #     phenotypes=f"{c6_dir}/PatientData/2022-Oct-25/c6_questionnaire_2022-10-25_19_00_59.csv",
    #     patient_hpo=f"{c6_dir}/PatientData/2022-Oct-25/c6_research_patients_2022-10-25_19_49_06.csv",
    #     geneset_gtf=f"{c6_dir}/GeneSets/hg19.ensGene.chr6.gtf.gz",
    #     drop_list_file=f"{c6_dir}/PatientData/drop_list.txt",
    #     hpo_termset_yaml=f"{c6_dir}/Phenotype_Homogeneity/selected_phenotypes_hpos.yaml",
    #     expand_hpos=False
    #     )

    my_comparison, my_geneset, my_ontology, termset = analyze_online(
        username=sys.argv[1], password=sys.argv[2], drop_list_file=sys.argv[3],
        hpo_termset_yaml=get_default_termset_yaml_path()
        )
