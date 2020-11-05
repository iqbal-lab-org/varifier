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

    # test dnadiff.make_truth_vcf
    dnadiff.make_truth_vcf(ref_fasta, truth_fasta, tmp_out, snps_only=False)
    expect_vcf = os.path.join(data_dir, "run.expect.vcf")
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)

    # now only with snps
    dnadiff.make_truth_vcf(ref_fasta, truth_fasta, tmp_out, snps_only=True)
    expect_vcf = os.path.join(data_dir, "run.expect.snps_only.vcf")
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)
    tmp_discarded_not_snps = "tmp.dnadiff.vcf.discarded_not_snps.vcf"
    expect_discarded_not_snps = os.path.join(data_dir, "run.expect.snps_only.discarded_not_snps.vcf")
    assert filecmp.cmp(tmp_discarded_not_snps, expect_discarded_not_snps, shallow=False)
    os.unlink(tmp_discarded_not_snps)

    # and now in rev comp
    dnadiff.make_truth_vcf(ref_fasta, truth_revcomp_fasta, tmp_out, snps_only=False)
    expect_revcomp_vcf = os.path.join(data_dir, "run.expect.revcomp.vcf")
    assert filecmp.cmp(tmp_out, expect_revcomp_vcf, shallow=False)
    os.unlink(tmp_out)

