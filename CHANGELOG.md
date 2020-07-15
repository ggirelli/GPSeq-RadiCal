# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- The input metadata table is now automatically copied to the output folder.
- Assert for `fpath` column in input metadata file.

### Changed
- The default `--bin-tags` parameter now include only `1e6:1e5,1e5`.

### Fixed
- Centrality estimation is now using the proper ratio of condition pairs.

## [0.0.1] - 2020-07-15

[Unreleased]: https://github.com/ggirelli/gpseq-img-py  
[0.0.1]: https://github.com/ggirelli/gpseq-radical/releases/tag/v0.0.1
