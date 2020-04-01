# varifier

Note: full documentation is under construction

## Installation

Dependencies:

* Python 3 (tested on version 3.6.9)
* mummer installed
* `paftools.js` and `k8` in your path. See https://github.com/lh3/minimap2/tree/master/misc

Install:

```
pip3 install .
```

## Usage

To verify calls in a VCF file, you will need:

1. `test.vcf`  - the VCF file to be tested
2. `ref.fasta` - FASTA file of reference corresponding to the VCF file
3. `truth.fasta` - a truth genome FASTA file

Run:
```
varifier vcf_eval truth.fasta ref.fasta test.vcf out_dir
```

This makes a new directory called `out_dir`. The results are in the file
`summary_stats.json`.

## Tests

To run the tests, run `tox` from the root of the repository.

