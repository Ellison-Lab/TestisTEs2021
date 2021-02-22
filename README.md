# Snakemake workflow: TestisTEs2021

This workflow generates outputs for the Ellison Lab Testes TE project.

## Get started

### prereqs

conda/mamba
bioconda in channels
cellranger >= 3.1.0
snakemake >= 5.32.2
peppy == 0.30.2
pandas == 1.1.0

Additionally `garnett`, which is used for cell type classification, doesn't do well with `conda` installation.
This must be installed per the recommendations at https://cole-trapnell-lab.github.io/garnett/docs_m3/#installing-garnett.

**NOTE**
Currently, garnett installation is troublesome, so may give trouble when running on
other systems. Working on finding a good installation option.

### command

```
PYTHONHASHSEED=0 snakemake --use-conda --profile config/profile/{your system} --cores {n} -kp
```

## Authors

* Matt Lawlor (@mal2017)
