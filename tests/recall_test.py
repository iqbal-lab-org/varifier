import filecmp
import os
import pytest
import subprocess

from cluster_vcf_records import vcf_record

from varifier import recall, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "recall")


def test_vcf_file_to_dict():
    vcf_file = os.path.join(data_dir, "vcf_file_to_dict.vcf")
    expect = {
        "ref1": [vcf_record.VcfRecord("ref1\t42\t1\tT\tA\t.\tPASS\t.\tGT\t1/1")],
        "ref2": [
            vcf_record.VcfRecord("ref2\t43\t3\tT\tA,C\t.\tPASS\t.\tGT\t2/2"),
            vcf_record.VcfRecord("ref2\t44\t2\tT\tA\t.\tPASS\t.\tGT\t1/1"),
        ],
    }
    got = recall._vcf_file_to_dict(vcf_file)
    assert got == expect


def test_apply_variants_to_genome():
    ref_fasta = os.path.join(data_dir, "apply_variants_to_genome.ref.fa")
    vcf_file = os.path.join(data_dir, "apply_variants_to_genome.vcf")
    expect_fasta = os.path.join(data_dir, "apply_variants_to_genome.expect.fa")
    tmp_fasta = "tmp.apply_variants_to_genome.fa"
    subprocess.check_output(f"rm -f {tmp_fasta}", shell=True)
    recall.apply_variants_to_genome(ref_fasta, vcf_file, tmp_fasta)
    assert filecmp.cmp(tmp_fasta, expect_fasta, shallow=False)
    os.unlink(tmp_fasta)


def test_get_recall():
    ref_fasta = os.path.join(data_dir, "get_recall.ref.fa")
    truth_fasta = os.path.join(data_dir, "get_recall.truth.revcomp.fa")
    vcf_to_test = os.path.join(data_dir, "get_recall.to_test.vcf")
    tmp_out = "tmp.get_recall"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_vcf = recall.get_recall(
        ref_fasta, vcf_to_test, tmp_out, 100, debug=True, truth_fasta=truth_fasta,
    )
    expect_vcf = os.path.join(data_dir, "get_recall.expect.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    # Same again, but with a mask that removes a few variants
    mask = {"truth": set(list(range(320, 391)))}
    mask["truth"].add(180)
    got_vcf = recall.get_recall(
        ref_fasta,
        vcf_to_test,
        tmp_out,
        100,
        debug=True,
        truth_fasta=truth_fasta,
        truth_mask=mask,
    )
    expect_vcf = os.path.join(data_dir, "get_recall.expect.masked.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
