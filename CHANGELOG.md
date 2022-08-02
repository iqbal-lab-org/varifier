# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add command line option `--hp_min_fix_length` to fix homopolymers of at
  least the given length, so that they match the reference sequence.

### Fixed

- Do not run `bcftools norm` when using global align, because it could
  result in inconsistency between VCF and output sequences

- Fix bug in global alignment, where alignment was forced to end at the
  end of the query sequence. Instead, allow it to stop earlier otherwise
  had some strange alignments that were incorrect (ending with a long gap and
  then a single aligned base, instead of aligned base, then a long gap)

## 0.4.0

### Added

- Added new command line option `--global_align` to both command line tasks
  `vcf_eval` and `make_truth_vcf`. This is meant for small (eg covid) genomes,
  where both truth and reference are one sequence only. It is less efficient
  but more accurate.

- Added command line options that only apply when `--global_align` is used:
  `--global_align_min_coord`, `--global_align_max_coord`,
  `--sanitise_truth_gaps`, `--use_non_acgt`,

- If using `--global_align`, writes an MSA of rev and truth sequences, and
  a FASTA of the truth sequence, but with santised gap lengths.

- New tag in VCF file `VCF_QRY_VARIANT`, which is the variant with respect to
  the query sequence, as opposed to the ref sequence.


### Fixed

- various edge cases caught when there are indels and/or Ns in sequences

- edge case caused by indels when making VCF with `--global align`.


[Unreleased]: https://github.com/iqbal-lab-org/varifier/compare/v0.4.0...HEAD
