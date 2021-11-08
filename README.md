# varifier

Note: full documentation is under construction

## Installation

### `conda`

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/varifier)](https://anaconda.org/bioconda/varifier)
[![bioconda version](https://anaconda.org/bioconda/varifier/badges/platforms.svg)](https://anaconda.org/bioconda/varifier)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```shell
$ conda install varifier
```

[channels]: https://bioconda.github.io/user/install.html#set-up-channels
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

### Container

Docker images are hosted at [quay.io].

#### `singularity`

Prerequisite: [`singularity`][singularity]

```shell
$ URI="docker://quay.io/iqballab/varifier"
$ singularity exec "$URI" varifier --help
```

The above will use the latest version. If you want to specify a version/commit then use a
[tag][quay.io] (or commit) like so.

```shell
$ TAG="3c8152a"
$ URI="docker://quay.io/iqballab/varifier:${TAG}"
```

#### `docker`

[![Docker Repository on Quay](https://quay.io/repository/iqballab/varifier/status "Docker Repository on Quay")](https://quay.io/repository/iqballab/varifier)

Prerequisite: [`docker`][docker]

```shhell
$ docker pull quay.io/iqballab/varifier
$ docker run quay.io/iqballab/varifier varifier --help
```

You can find all the available tags on the [quay.io repository][quay.io].

[quay.io]: https://quay.io/repository/iqballab/varifier
[singularity]: https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps
[docker]: https://docs.docker.com/v17.12/install/

### Local

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

