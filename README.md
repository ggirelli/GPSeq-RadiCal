# GPSeq-RadiCal

[![DOI](https://zenodo.org/badge/279701571.svg)](https://zenodo.org/badge/latestdoi/279701571)

R script for fast and streamlined **Radi**ality **Cal**culation from **GPSeq** data. GPSeq-RadiCal aims to provide the same functionalities as the [gpseqc](https://github.com/ggirelli/gpseqc) Python3 package, but with much lower memory usage and faster computation times.

## Installation

* Tested on R v3.6.3.
* Required packages: argparser, data.table, logging, outliers, pbapply, rtracklayer.

The script checks automatically for the required packages.  
When packages are missing, it provides the code to install the missing ones, one by one.

## Usage

First, prepare a tabulation-separated metadata file with four columns:

* `exid`: sequencing run ID.
* `cond`: condition description (e.g., time and/or concentration).
* `libid`: library ID.
* `fpath`: full absolute path to bed file. Bed files can be gzipped (recommended).

The bed files should be reported in order of condition strength, with the top rows being *weaker* than bottom ones. Also, the first line should contain the column headers. An example metadata file can be found [here](example_meta.tsv).

Then, to run with default parameters, execute the following command:

```
./gpseq-radical.R example_meta.tsv output_folder
```

Replacing `example_meta.tsv` with your metadata file path, and `output_folder` with the path where the script should write the output. Note that the specified output folder must not already exist.

To access the script help page, run the following:

```
./gpseq-radical.R -h
```

## Differences from `gpseqc`

* Allows for normalization in a chromosome-wise fashion too, instead of only library-wise.
* Normalization factors are calculated **after** outlier removal, not before.
* Masking is performed **before** binning, but **after** the normalization factors are calculated. In this way, reads from masked regions still count for normalization purposes, while they are masked out from chromosome-wide analysis.

## Desired features

* Re-running on the same output folder to populate with additional resolutions.
* Additional centrality estimates (currently calculates only the one selected in the GPSeq study).

## Contributing

We welcome any contributions to `GPSeq-RadiCal`. Please, refer to the [contribution guidelines](CONTRIBUTING.md) if this is your first time contributing! Also, check out our [code of conduct](CODE_OF_CONDUCT.md).

## License

```
MIT License
Copyright (c) 2020 Gabriele Girelli
```

## References

* Girelli, Gabriele, et al. "GPSeq reveals the radial organization of chromatin in the cell nucleus." Nature Biotechnology (2020): 1-10. **Genomic loci Positioning by Sequencing (GPSeq) original paper**.
