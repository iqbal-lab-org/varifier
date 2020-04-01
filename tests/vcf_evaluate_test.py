import copy
import filecmp
import os
import pytest
import subprocess

from varifier import vcf_evaluate

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_evaluate")


def test_add_overall_precision_and_recall_to_summary_stats():
    stats = {
        "Precision": {
            "ALL": {"TP": {"Count": 9}, "FP": {"Count": 1}},
            "FILT": {"TP": {"Count": 5}, "FP": {"Count": 0}},
        },
        "Recall": {
            "ALL": {"TP": {"Count": 4}, "FP": {"Count": 1}},
            "FILT": {"TP": {"Count": 0}, "FP": {"Count": 0}},
        },
    }
    expect = copy.deepcopy(stats)
    vcf_evaluate._add_overall_precision_and_recall_to_summary_stats(stats)
    assert stats != expect
    expect["Precision"]["ALL"]["Precision"] = 0.9
    expect["Precision"]["FILT"]["Precision"] = 1
    expect["Recall"]["ALL"]["Recall"] = 0.8
    expect["Recall"]["FILT"]["Recall"] = 0
    assert stats == expect


def test_evaluate_vcf():
    mask_bed_file = os.path.join(data_dir, "evaluate_vcf.mask.bed")
    truth_fasta = os.path.join(data_dir, "evaluate_vcf.truth.fa")
    ref_fasta = os.path.join(data_dir, "evaluate_vcf.ref.fa")
    vcf_to_eval = os.path.join(data_dir, "evaluate_vcf.to_eval.vcf")
    tmp_out = "tmp.vcf_evaluate.evaluate_vcf.out"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)

    vcf_evaluate.evaluate_vcf(
        vcf_to_eval, ref_fasta, truth_fasta, 100, tmp_out, debug=True, force=True
    )
    summary_stats_expect_json = os.path.join(
        data_dir, "evaluate_vcf.expect.summary_stats.json"
    )
    summary_stats_got_json = os.path.join(tmp_out, "summary_stats.json")
    assert filecmp.cmp(summary_stats_got_json, summary_stats_expect_json, shallow=False)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    vcf_evaluate.evaluate_vcf(
        vcf_to_eval, ref_fasta, truth_fasta, 100, tmp_out, debug=True, force=True, mask_bed_file=mask_bed_file
    )
    summary_stats_expect_json = os.path.join(
        data_dir, "evaluate_vcf.expect.masked.summary_stats.json"
    )
    summary_stats_got_json = os.path.join(tmp_out, "summary_stats.json")
    assert filecmp.cmp(summary_stats_got_json, summary_stats_expect_json, shallow=False)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
