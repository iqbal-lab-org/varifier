# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added new command line option `--global_align` to both command line tasks
  `vcf_eval` and `make_truth_vcf`. This is meant for small (eg covid) genomes,
  where both truth and reference are one sequence only. It is less efficient
  but more accurate.

- New tag in VCF file `VCF_QRY_VARIANT`, which is the variant with respect to
  the query sequence, as opposed to the ref sequence.


### Fixed

- various edge cases caught when there are indels and/or Ns in sequences


[Unreleased]: https://github.com/Mykrobe-tools/mykrobe/compare/v0.10.0...HEAD
https://github.com/iqbal-lab-org/varifier/compare/v0.3.1...HEAD
