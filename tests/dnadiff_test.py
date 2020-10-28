import filecmp
import os
import pytest
import subprocess

from varifier import dnadiff

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "dnadiff")


def test_run_dnadiff():
    ref_fasta = os.path.join(data_dir, "run_dnadiff.ref.fa")
    query_fasta = os.path.join(data_dir, "run_dnadiff.query.fa")
    expect = os.path.join(data_dir, "run_dnadiff.snps")
    tmp_out = "tmp.run_dnadiff.snps"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    dnadiff._run_dnadiff(ref_fasta, query_fasta, tmp_out, split_query=False)
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)
    dnadiff._run_dnadiff(ref_fasta, query_fasta, tmp_out, split_query=True)
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


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
