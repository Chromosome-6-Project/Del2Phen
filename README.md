# Del2Phen

## Introduction

`Del2Phen` is a tool for predicting phenotypes for a given copy-number variant (CNV), based on known phenotypes from patients with similar CNVs. Patients are first grouped according to minimum similarity thresholds across multiple CNV properties, including gene and haploinsufficient genes affected by the CNV. Phenotypes are then predicted by the prevalence of phenotypes found in individuals in the group.

This project is part of The Chromosome 6 Project, a research initiative that aims to provide information on rare chromosome 6 aberrations in children to parents and healthcare professionals.


## Installation
`Del2Phen` is written in Python 3 and requires Python 3.9. We recommend using `conda` to create a Python 3.9 environment, then installing with `pip`:

```python -m pip install del2phen```


## Getting Started
At minimum, `Del2Phen` requires basic patient CNV and phenotype tabular data.

### CNV data
CNV input data must be provided as a 5-column tab-separated file. Each row should specify one CNV from one patient. Patients with multiple CNVs should be listed with the same ID across multiple lines, one CNV per line. The file must include the following headings:

- id: Patient ID of the CNV. The same ID can be used multiple times to indicate that one patient has multiple CNVs.
- chromosome: Chromosome on which the CNV is located. These chromosome names must match the names used in the geneset GTF file used.
- start: 1-indexed inclusive CNV start position.
- stop: 1-indexed inclusive CNV stop position.
- copy_number: Integer value of the CNV, e.g., a single copy deletion is 1, a single copy duplication is 3.

### Phenotype data
Patient phenotype data must be provided as a tab-separated file, where the first column contains patient IDs with the header `id`. Each additional column must have as a header either a [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/) phenotype ID (e.g., HP:0001643) or any other identifier which must then be defined as a custom phenotype term. Each entry in each column should be coercible to a Boolean or NA (t, T, true, True, f, F, false, False, 0, 1, NA) representing whether or not the patient exhibits the phenotype. Patients are not required to be present in the phenotype file to be present in the CNV file.


### Prediction
There are several ways to produce phenotype predictions with Del2Phen. The simplest way to predict phenotypes for a single query CNV is the following:

```del2phen -g cnvs.tsv -p phenotypes.tsv -cnv chr6:123456-234567:1 -op ./```

The above command will produce a table of phenotypes predicted to be found in patients with a deletion (copy number of 1) on chromosome 6 from bases 123456 to 234567. 
