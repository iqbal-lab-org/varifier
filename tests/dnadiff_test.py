import filecmp
import os
import pytest

from varifier import dnadiff

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "dnadiff")


def test_make_truth_vcf():
    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    truth_fasta = os.path.join(data_dir, "run.truth.fa")
    truth_revcomp_fasta = os.path.join(data_dir, "run.truth.revcomp.fa")
    tmp_out = "tmp.dnadiff.vcf"
    if os.path.exists(tmp_out):
        os.unlink(tmp_out)
    dnadiff.make_truth_vcf(ref_fasta, truth_fasta, tmp_out)
    expect_vcf = os.path.join(data_dir, "run.expect.vcf")
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)
    dnadiff.make_truth_vcf(ref_fasta, truth_revcomp_fasta, tmp_out)
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)
