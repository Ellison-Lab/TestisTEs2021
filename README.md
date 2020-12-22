# Snakemake workflow: TestisTEs2021

This workflow generates outputs for the Ellison Lab Testes TE project.

## Get started

### prereqs

conda
bioconda in channels
cellranger >= 3.1.0
snakemake >= 5.30.1

Additionally `garnett`, which is used for cell type classification, doesn't do well with `conda` installation.
This must be installed per the recommendations at https://cole-trapnell-lab.github.io/garnett/docs_m3/#installing-garnett.

### command

```
snakemake --use-conda --profile config/profile/{your system} -kp
```

## Authors

* Matt Lawlor (@mal2017)
