# Needler

Code and algorithm to support our paper Needler: An Algorithm to Develop a Comprehensive Targeted MS Method Capable of Monitoring the Human Proteome.

## Steps

Makefile should replicate the needed steps assuming a standard Linux environment (`make`, `wget`, etc) + [Python Poetry](https://python-poetry.org/) are available.
Alternatively, can use the provided .devcontainer which will create an environment with necessary tooling.

Makefile high-level outline:

- `make build_env` Configure the Python virtual environment with Poetry
- `make download` Download required datasets
  - Uniprot Swiss-Prot proteome fasta
  - [A deep proteome and transcriptome abundance atlas of 29 healthy human tissues](https://doi.org/10.15252/msb.20188503) proteomics reference for detectable proteins per tissue.
- `make munge` pre-process the datasets
  - extract liver proteins from Tissue Atlas study
  - extract human proteins from Uniprot
  - digest proteins into tryptic peptides
  - filter tryptic peptides
    - sized 5-30 residues
    - distinct (represented by single protein)
    - does not contain uncommon amino acid (X or U)
    - does not contain methionine (commonly oxidized and is inappropriate for MS quantification)
    - each protein represented by >= 2 peptides (common MS quantitation criteria)
  - assign predicted iRT and retention time value to all retained peptide sequences
  - produce sub-proteomes for study: liver, kinase, dub
- `make fit` run the needler algorithm to produce targeted MS methods for each proteome. Recommend invoking with a `-j <CPU_COUNT>` option to run fits in parallel. This represents a significant amount of computation and is not advised to run on a single machine as it represents potentially years of dedicated computer time.

## Data

Datasets contained within the repo itself:

- `procal_supplementaltable_S1.csv`. Supplementary data table from [PROCAL: A Set of 40 Peptide Standards for Retention Time Indexing, Column Performance Monitoring, and Collision Energy Calibration](https://doi.org/10.1002/pmic.201700263) of the 40 retention time standard peptide sequences chosen by the Kuster group.
- `prosit_predict_irt.csv` predicted indexed retention times (iRT) values for human peptide sequences using the [Prosit: proteome-wide prediction of peptide tandem mass spectra by deep learning](https://doi.org/10.1038/s41592-019-0426-7) Proteome Tools calculator. Configuring the calculator takes a bit of effort, so committing computed values here.
- `uniprot_pkinfam_202004.csv` human kinases extracted from the [Uniprot pkinfam resource](https://www.uniprot.org/docs/pkinfam)

## License

Code is available under an Apache 2 license.
