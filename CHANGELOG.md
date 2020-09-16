# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- If duplicated regions are found in input bed file, a warning is triggered and the minimum score is used.

## [0.0.6] - 2020-07-20
### Added
- `--chrom-tag` option, to define the number of chromosomes and heterosome names.

### Changed
- Using `--chrom-tag` to select chromosomes from retrieved/read chromosome information.

### Fixed
- Now retrieving correct chromosome size information from UCSC, if not provided.
- All bins are now retained in the output, even if no sites/reads are present or the estimated centrality is NaN.

## [0.0.5] - 2020-07-16
### Added
- `--bin-bed` option to focus calculation on regions of interest (e.g., FISH probes).
- Now checking for mask bed consistency.
- `--version` option for easier version tracking.
- `--debug-info` option for easier debugging.

### Changed
- Now masking before binning (on dcasted input bed data), but after the normalization factors are calculated.

### Fixed
- Options and log filenames.

## [0.0.4] - 2020-07-15
### Added
- Additional merged rds output.

### Fixed
- Solved score masking crashes due to missing argument passing.

## [0.0.3] - 2020-07-15
### Changed
- Now using `exid` from input metadata and running one experiment at a time.
- The default `--bin-tags` parameter now include `1e6:1e5,1e5:1e4`.

## [0.0.2] - 2020-07-15
### Added
- Assert for `fpath` column in input metadata file.
- Assert for at least 2 bed files in input.
- The input metadata table is now automatically copied to the output folder.
- Help page description.
- Full detail help apge with `--more-help`.

### Changed
- The default `--bin-tags` parameter now include only `1e6:1e5,1e5`.
- It is now possible to skip input bed outlier removal by using `--bed-outlier-tag ""`.
- It is now possible to skip score rescaling by using `--score-outlier-tag ""`.

### Fixed
- Centrality estimation is now using the proper ratio of condition pairs.
- Using specified score outlier tag for rescaling purposes, previously using bed outlier tag due to a bug.

## [0.0.1] - 2020-07-15

[Unreleased]: https://github.com/ggirelli/gpseq-img-py  
[0.0.6]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.6
[0.0.5]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.5
[0.0.4]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.4
[0.0.3]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.3
[0.0.2]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.2
[0.0.1]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.1
