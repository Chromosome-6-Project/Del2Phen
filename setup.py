
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="Del2Phen",
    version="2024.7.15",
    author="T.D. Medina",
    author_email="tylerdanmedina@gmail.com",
    description="Del2Phen, a tool to predict phenotypes associated with copy-number "
                "variants based on similar CNVs in patient data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Chromosome-6-Project/Del2Phen",
    packages=["del2phen",
              "del2phen.analysis",
              "del2phen.dashboard",
              "del2phen.resources"],
    # package_dir={"": "del2phen"},
    package_data={"del2phen.resources": ["*"], "del2phen.dashboard": ["assets/*"]},
    python_requires=">=3.9, !=3.10",
    install_requires=[
        "dash",
        "dash-bootstrap-components",
        "mygene",
        "networkx",
        "numpy",
        "pandas",
        "plotly",
        "pronto",
        "PyYAML",
        "venn",
        ],
    entry_points={
        "console_scripts": [
            "del2phen = del2phen.analysis.__main__:main",
            # "del2phen-dash = del2phen.dashboard.__main__:main"
            ]
        }
    )
