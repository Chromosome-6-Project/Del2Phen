
from chr6_project.analysis.analyze import analyze
from chr6_project.analysis.phenotype_homogeneity import phenotype_homo_test


if __name__ == "__main__":
    my_comparison, my_geneset, my_ontology = analyze(
        genotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_array_2022-02-21_18_28_06.csv",
        phenotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_questionnaire_2022-02-21_10_18_09.csv",
        patient_hpo="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_research_patients_2022-02-21_10_20_53.csv",
        geneset_gtf="/home/tyler/Documents/Chr6_docs/GeneSets/hg19.ensGene.chr6.gtf.gz",
        drop_list_file="/home/tyler/Documents/Chr6_docs/PatientData/drop_list.txt",
        expand_hpos=False
        )
    table, selected_hpos, homos = phenotype_homo_test(my_comparison, my_ontology, .85, .2, 2)
