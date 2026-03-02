# <sub><img width="64" height="64" alt="pure-faf logo" src="https://github.com/user-attachments/assets/4d3d9fa3-e90b-4c81-94b9-339cef5ba99d" /></sub> Snakemake workflow: `pure-faf`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/rabiafidan/pure-faf/workflows/Tests/badge.svg?branch=main)](https://github.com/rabiafidan/pure-faf/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/rabiafidan/pure-faf)

A Snakemake workflow for `VCF purification through formalin artefact filtering`

<img width="3840" height="2160" alt="pure-faf_workflow_updated" src="https://github.com/user-attachments/assets/395d9bd1-63a3-45c8-a570-49e48a104654" />


- [Snakemake workflow: `pure-faf`](#pure-faf)
  - [Documentation](#documentation)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Documentation 

See the [documentation](https://rabiafidan.github.io/pure-faf) for detailed information.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/rabiafidan/pure-faf).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

### 0) Dependencies
You need `snakemake`, `cyvcf2`, `conda` and `singularity` to run this workflow.

```bash
conda create --name pure-faf -c conda-forge -c bioconda snakemake singularity cyvcf2
conda activate pure-faf
```
If you are running this workflow on an HPC environment, make sure to also install the correct executor plugin. For example:

```
pip install snakemake-executor-plugin-slurm
```
### 1) clone the repository
```bash
git clone https://github.com/rabiafidan/pure-faf.git
```
To run the workflow from command line, change the working directory.

```bash
cd path/to/pure-faf
```

## 2) configure the pipeline
Adjust options in the default config file `config/config.yaml`.


Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```
To run the workflow with test files using **conda** and **singularity**:

```bash
snakemake --cores 2 --use-conda --use-singularity
```

## Authors

- Fatma Rabia Fidan Baş
  - The Francis Crick Institute
  - [ORCID profile](https://orcid.org/0000-0001-5877-7945)

- Brian Hanley
  - The Francis Crick Institute
  - [ORCID profile](https://orcid.org/0000-0001-5877-7945)

- Husayn Pallikonda
  - The Francis Crick Institute
  - [ORCID profile](https://orcid.org/0000-0001-5877-7945)

- Irene Lobon 
  - AstraZeneca
  - [ORCID profile](https://orcid.org/0000-0003-1170-9915)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO


