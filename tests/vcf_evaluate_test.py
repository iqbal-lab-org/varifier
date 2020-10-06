import copy
import filecmp
import os
import pytest
import subprocess

from varifier import utils, vcf_evaluate

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_evaluate")


def test_add_overall_precision_and_recall_to_summary_stats():
    stats = {
        "Precision": {
            "TP": {"Count": 9, "SUM_ALLELE_MATCH_FRAC": 9.0},
            "EDIT_DIST_COUNTS": {"numerator": 18, "denominator": 20},
            "FP": {"Count": 1, "SUM_ALLELE_MATCH_FRAC": 0.99},
        },
        "Recall": {
            "TP": {"Count": 5, "SUM_ALLELE_MATCH_FRAC": 5.0},
            "EDIT_DIST_COUNTS": {"numerator": 15, "denominator": 20},
            "FN": {"Count": 1, "SUM_ALLELE_MATCH_FRAC": 0.2},
        },
    }
    expect = copy.deepcopy(stats)
    vcf_evaluate._add_overall_precision_and_recall_to_summary_stats(stats)
    assert stats != expect
    expect["Precision"]["Precision"] = 0.9
    expect["Recall"]["Recall"] = 0.83333333
    expect["Precision"]["Precision_frac"] = 0.999
    expect["Recall"]["Recall_frac"] = 0.86666667
    expect["Precision"]["Precision_edit_dist"] = 0.9
    expect["Recall"]["Recall_edit_dist"] = 0.75
    assert stats == expect


def test_filter_vcf():
    infile = os.path.join(data_dir, "filter_vcf.in.vcf")
    got_keep = "tmp.filter_vcf.keep.vcf"
    got_exclude = "tmp.filter_vcf.exclude.vcf"
    subprocess.check_output(f"rm -rf {got_keep} {got_exclude}", shell=True)
    ref_seqs = {"ref": "ATGCATGACTGCATTACTCATCATCGAATG"}
    got_counts = vcf_evaluate._filter_vcf(
        infile, got_keep, got_exclude, ref_seqs, filter_pass=None, keep_ref_calls=True
    )
    expect_counts = {
        "filter_fail": 1,
        "heterozygous": 1,
        "no_genotype": 1,
        "ref_call": 0,
        "other": 3,
    }
    assert got_counts == expect_counts
    expect_keep = os.path.join(
        data_dir, "filter_vcf.expect.no_filter_pass_keep_ref_calls.keep.vcf"
    )
    expect_exclude = os.path.join(
        data_dir, "filter_vcf.expect.no_filter_pass_exclude_ref_calls.exclude.vcf"
    )
    utils.vcf_records_are_the_same(got_keep, expect_keep)
    utils.vcf_records_are_the_same(got_exclude, expect_exclude)
    os.unlink(got_keep)
    os.unlink(got_exclude)

    got_counts = vcf_evaluate._filter_vcf(
        infile,
        got_keep,
        got_exclude,
        ref_seqs,
        filter_pass={"PASS"},
        keep_ref_calls=False,
    )
    expect_counts = {
        "filter_fail": 4,
        "heterozygous": 1,
        "no_genotype": 1,
        "ref_call": 1,
        "other": 3,
    }
    assert got_counts == expect_counts
    expect_keep = os.path.join(data_dir, "filter_vcf.expect.with_filtering.1.keep.vcf")
    expect_exclude = os.path.join(
        data_dir, "filter_vcf.expect.with_filtering.1.exclude.vcf"
    )
    utils.vcf_records_are_the_same(got_keep, expect_keep)
    utils.vcf_records_are_the_same(got_exclude, expect_exclude)
    os.unlink(got_keep)
    os.unlink(got_exclude)

    got_counts = vcf_evaluate._filter_vcf(
        infile,
        got_keep,
        got_exclude,
        ref_seqs,
        filter_pass={".", "FILTER_2", "PASS"},
        keep_ref_calls=False,
    )
    expect_counts = {
        "filter_fail": 2,
        "heterozygous": 1,
        "no_genotype": 1,
        "ref_call": 1,
        "other": 3,
    }
    assert got_counts == expect_counts
    expect_keep = os.path.join(data_dir, "filter_vcf.expect.with_filtering.2.keep.vcf")
    expect_exclude = os.path.join(
        data_dir, "filter_vcf.expect.with_filtering.2.exclude.vcf"
    )
    utils.vcf_records_are_the_same(got_keep, expect_keep)
    utils.vcf_records_are_the_same(got_exclude, expect_exclude)
    os.unlink(got_keep)
    os.unlink(got_exclude)


def test_evaluate_vcf():
    truth_fasta = os.path.join(data_dir, "evaluate_vcf.truth.fa")
    ref_fasta = os.path.join(data_dir, "evaluate_vcf.ref.fa")
    vcf_to_eval = os.path.join(data_dir, "evaluate_vcf.to_eval.vcf")
    tmp_out = "tmp.vcf_evaluate.evaluate_vcf.out"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)

    vcf_evaluate.evaluate_vcf(
        vcf_to_eval,
        ref_fasta,
        truth_fasta,
        100,
        tmp_out,
        debug=True,
        force=True,
        filter_pass={"PASS"},
    )
    summary_stats_expect_json = os.path.join(
        data_dir, "evaluate_vcf.expect.summary_stats.json"
    )
    summary_stats_got_json = os.path.join(tmp_out, "summary_stats.json")
    assert filecmp.cmp(summary_stats_got_json, summary_stats_expect_json, shallow=False)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    ref_mask_bed_file = os.path.join(data_dir, "evaluate_vcf.ref_mask.bed")
    truth_mask_bed_file = os.path.join(data_dir, "evaluate_vcf.truth_mask.bed")
    vcf_evaluate.evaluate_vcf(
        vcf_to_eval,
        ref_fasta,
        truth_fasta,
        100,
        tmp_out,
        debug=True,
        force=True,
        filter_pass={"PASS"},
        ref_mask_bed_file=ref_mask_bed_file,
        truth_mask_bed_file=truth_mask_bed_file,
    )
    summary_stats_expect_json = os.path.join(
        data_dir, "evaluate_vcf.expect.masked.summary_stats.json"
    )
    summary_stats_got_json = os.path.join(tmp_out, "summary_stats.json")
    assert filecmp.cmp(summary_stats_got_json, summary_stats_expect_json, shallow=False)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
