# Ellison Lab Drosophila Testes TE project

## Requirements

#### Software for data processing

For processing of raw data, the pipeline requires

* **Cellranger** v3.1.0
* **snakemake** v6.1.0
* **mamba** v0.9.1

Each of these may be installed in a 3-5 minutes with a standard internet connection.

Most other dependencies for processing of raw data are automatically handled by
**snakemake**/**mamba** during a run when called with the `--use-conda` flag (see below for running instructions).

For this reason, conda channels must be configured according to instructions available
at the [bioconda channel page](https://bioconda.github.io/user/install.html#set-up-channels).

#### Software for figures

To generate figures, R packages must be installed manually. Altogether these packages
may be installed in approximately 10 minutes via `install.packages()` for CRAN hosted packages
or `BiocManager::install()` for Bioconductor packages.

* base **R** v4.0.3
* **tidyverse** v1.3.0
* **ragg** v1.1.2
* **ggtext** v0.1.1
* **arrow** v3.0.0
* **jsonlite** v1.7.2
* **readxl** v1.3.1
* **rtracklayer** v1.50.0
* **GenomicRanges** v1.42.0
* **patchwork** v1.1.1
* **png** v0.1-7
* **grid** v4.0.3
* **magick** v2.7.1
* **extrafont** v0.17
* **ggforce** v0.3.3
* **ggrepel** v0.9.1
* **GenomicFeatures** v1.42.3
* **ggpubr** v0.4.0
* **ggpointdensity** v0.1.0
* **ComplexHeatmap** v2.6.2
* **broom** v0.7.5
* **VariantAnnotation** v1.36.0

#### System Requirements

The raw data processing pipeline has been tested
on CentOS Linux 7. Alignment and quantification steps may be resource intensive.
Grid search steps are numerous. Therefore it is recommended to run on an HPC system
to make total running time reasonable. We ran on a **slurm** cluster. Please see the **snakemake**
[documentation on system profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to
set up **snakemake** to run on other cluster submission systems.

## Installation

The primary workflow can be cloned from github via

```bash
git clone https://github.com/Ellison-Lab/TestisTEs2021
```

Some large files in the `resources/` directory can be pulled with the command below.
You may need to install **git lfs**. See https://git-lfs.github.com/ or https://anaconda.org/conda-forge/git-lfs for details.

```bash
git lfs pull
```

The primary workflow orchestrates subworkflows and figure scripts.
The project relies on multiple git submodules.
These can be pulled into the main workflow via

```bash
git submodule update --init
git submodule foreach git pull origin main
```

Finally, absolute local paths in the `config.yaml` file and `subsamples_table.csv`'s
should be updated as appropriate.

Alternatively, these submodules can be separately cloned and run separately.

* https://github.com/Ellison-Lab/gte21-custom-genome
* https://github.com/Ellison-Lab/gte21-ica-grid
* https://github.com/Ellison-Lab/gte21-scrna
* https://github.com/Ellison-Lab/gte21-total-rna-fusions
* https://github.com/Ellison-Lab/gte21-chimeric-rnaseq
* https://github.com/Ellison-Lab/gte21-tidal-xa
* https://github.com/Ellison-Lab/gte21-te-variant-expression
* https://github.com/Ellison-Lab/transposon-variants-hts
* https://github.com/Ellison-Lab/gte21-y-hetchrom

Altogether, these repositories can be downloaded in 2-3 minutes.

As before, when repos are downloaded separately absolute local paths in the `config.yaml`
 file and `subsamples_table.csv`'s should be updated as appropriate.

## Demo/Instructions

This workflow is not a single standalone tool and is therefore not suited to demo'ing with a small dataset. However the full workflow can be run within 1-3 days on a HPC cluster with modest resource requests.

A dry run can be initiated via the code below. Exclude `-n` for actual run to produce arrow format files under `results/finalized`.
These hold processed data in a space efficient file format and are easily importable to R as dataframes via the **arrow** package)

Please note that specifying the python hashseed to 0 specifically is important
for **scanpy** reproducibility, per
[this issue report](https://github.com/theislab/scanpy/issues/313).

```bash
PYTHONHASHSEED=0 snakemake --use-conda --profile {your system} --cores {n} -kpn
```
Figures can be run via the below command, where {figure name} follows the pattern
"figure1" or "supp1." This will generate PDFs and `.rds` files under `results/panels`.

```bash
PYTHONHASHSEED=0 snakemake {figure name} --use-conda --profile {your system} --cores {n} -kpn
```
