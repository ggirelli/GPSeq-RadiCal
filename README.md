# GPSeq-RadiCal
R script for fast and streamlined **Radi**ality **Cal**culation from **GPSeq** data. GPSeq-RadiCal aims to provide the same functionalities as the [gpseqc](https://github.com/ggirelli/gpseqc) Python3 package, but with much lower memory usage and faster computation times.

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

## Desired features

* Re-running on the same output folder to populate with additional resolutions.
* Additional centrality estimates (currently calculates only the one selected in the GPSeq study).

## Contributing

We welcome any contributions to `GPSeq-RadiCal`. Please, refer to the [contribution guidelines](CONTRIBUTING.md) if this is your first time contributing! Also, check out our [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).

## License

```
MIT License
Copyright (c) 2020 Gabriele Girelli
```

## References

* Girelli, Gabriele, et al. "GPSeq reveals the radial organization of chromatin in the cell nucleus." Nature Biotechnology (2020): 1-10. **Genomic loci Positioning by Sequencing (GPSeq) original paper**.