import filecmp
import os
import pytest
import subprocess
from unittest import mock

from varifier import tasks, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "tasks")


def test_make_truth_vcf():
    options = mock.Mock()
    options.truth_mask = None
    options.ref_fasta = os.path.join(data_dir, "make_truth_vcf.ref.fa")
    options.truth_fasta = os.path.join(data_dir, "make_truth_vcf.truth.fa")
    options.debug = False
    options.outdir = "tmp.tasks.make_truth_vcf"
    options.flank_length = 100
    options.max_recall_ref_len = None
    subprocess.check_output(f"rm -rf {options.outdir}", shell=True)
    tasks.make_truth_vcf.run(options)
    got_vcf = os.path.join(options.outdir, "04.truth.vcf")
    expect_vcf = os.path.join(data_dir, "make_truth_vcf.expect.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)

    subprocess.check_output(f"rm -r {options.outdir}", shell=True)
    options.truth_mask = os.path.join(data_dir, "make_truth_vcf.mask.bed")
    tasks.make_truth_vcf.run(options)
    expect_vcf = os.path.join(data_dir, "make_truth_vcf.expect.masked.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {options.outdir}", shell=True)


def test_vcf_eval():
    options = mock.Mock()
    options.vcf_in = os.path.join(data_dir, "vcf_eval.to_eval.vcf")
    options.vcf_fasta = os.path.join(data_dir, "vcf_eval.ref.fa")
    options.truth_fasta = os.path.join(data_dir, "vcf_eval.truth.fa")
    options.flank_length = 100
    options.outdir = "tmp.tasks.vcf_eval"
    options.truth_vcf = None
    options.debug = False
    options.force = False
    options.ref_mask = None
    options.truth_mask = None
    options.use_ref_calls = False
    options.max_recall_ref_len = None
    options.filter_pass = "PASS,."
    subprocess.check_output(f"rm -rf {options.outdir}", shell=True)
    tasks.vcf_eval.run(options)
    expect_json = os.path.join(data_dir, "vcf_eval.expect.summary_stats.json")
    got_json = os.path.join(options.outdir, "summary_stats.json")
    assert filecmp.cmp(expect_json, got_json, shallow=False)
    subprocess.check_output(f"rm -r {options.outdir}", shell=True)

    options.ref_mask = os.path.join(data_dir, "vcf_eval.ref_mask.bed")
    tasks.vcf_eval.run(options)
    expect_json = os.path.join(data_dir, "vcf_eval.expect.masked.summary_stats.json")
    assert filecmp.cmp(expect_json, got_json, shallow=False)
    subprocess.check_output(f"rm -r {options.outdir}", shell=True)
