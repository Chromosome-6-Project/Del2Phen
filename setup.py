from setuptools import setup

setup(
    name="chr6_project",
    version="0.1.0",
    description="Chromosome 6 Project Automatic Phenotype Prediction Tool",
    url="https://github.com/Chromosome-6-Project/Automatic_Phenotype_Prediction",
    author="T.D. Medina",
    author_email="tylerdanmedina@gmail.com",
    packages=["chr6_project",
              "chr6_project.analysis",
              "chr6_project.dashboard",
              "chr6_project.molgenis_import",
              "chr6_project.resources"],
    # package_dir={"": "chr6_project"},
    package_data={"resources": ["*"], "dashboard.assets": ["*"]},
    python_requires=">=3.9",
    install_requires=[
        "dash==2.3.1",
        "dash-bootstrap-components==1.1.0",
        "molgenis-py-client==2.4.0",
        "mygene==3.2.2",
        "networkx==2.6.3",
        "numpy==1.23.2",
        "pandas==1.3.5",
        "plotly==5.7.0",
        "pronto==2.4.3",
        "PyYAML==6.0",
        "venn==0.1.3",
        ]
    )
