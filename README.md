# GPSeq-RadiCal
R script for fast and streamlined **Radi**ality **Cal**culation from **GPSeq** data. GPSeq-RadiCal aims to provide the same functionality as in the [gpseqc](https://github.com/ggirelli/gpseqc) Python3 package, but with much lower memory usage and faster computation times.

## Installation

* Tested on R v3.6.3.
* Required packages: argparser, data.table, logging, outliers, pbapply, rtracklayer.

The script checks automatically for the required packages.  
When packages are missing, it provides the code to install the missing ones, one by one.

## Usage

*coming soon*

## Differences from `gpseqc`

* Allows for normalization in a chromosome-wise fashion too, instead of only library-wise.
* Normalization factors are calculated **after** outlier removal, not before.

## Missing features

* Store input metadata in output folder.
* Re-running on the same output folder to populate with additional resolutions.
    - Should check for matching metadata.
* Additional centrality estimates (currently calculates only the one selected in the GPSeq study).

## Contributing

*coming soon*

## License

```
MIT License
Copyright (c) 2020 Gabriele Girelli
```
