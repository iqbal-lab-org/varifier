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
            "ALL": {
                "TP": {"Count": 9, "SUM_ALLELE_MATCH_FRAC": 9.0},
                "EDIT_DIST_COUNTS": {"numerator": 18, "denominator": 20},
                "FP": {"Count": 1, "SUM_ALLELE_MATCH_FRAC": 0.99},
            },
            "FILT": {
                "TP": {"Count": 5, "SUM_ALLELE_MATCH_FRAC": 5.0},
                "EDIT_DIST_COUNTS": {"numerator": 16, "denominator": 18},
                "FP": {"Count": 0, "SUM_ALLELE_MATCH_FRAC": 0},
            },
        },
        "Recall": {
            "ALL": {
                "TP": {"Count": 5, "SUM_ALLELE_MATCH_FRAC": 5.0},
                "EDIT_DIST_COUNTS": {"numerator": 15, "denominator": 20},
                "FN": {"Count": 1, "SUM_ALLELE_MATCH_FRAC": 0.2},
            },
            "FILT": {
                "TP": {"Count": 0, "SUM_ALLELE_MATCH_FRAC": 0.0},
                "EDIT_DIST_COUNTS": {"numerator": 14, "denominator": 15},
                "FN": {"Count": 0, "SUM_ALLELE_MATCH_FRAC": 0.0},
            },
        },
    }
    expect = copy.deepcopy(stats)
    vcf_evaluate._add_overall_precision_and_recall_to_summary_stats(stats)
    assert stats != expect
    expect["Precision"]["ALL"]["Precision"] = 0.9
    expect["Precision"]["FILT"]["Precision"] = 1
    expect["Recall"]["ALL"]["Recall"] = 0.83333333
    expect["Recall"]["FILT"]["Recall"] = 0
    expect["Precision"]["ALL"]["Precision_frac"] = 0.999
    expect["Precision"]["FILT"]["Precision_frac"] = 1.0
    expect["Recall"]["ALL"]["Recall_frac"] = 0.86666667
    expect["Recall"]["FILT"]["Recall_frac"] = 0
    expect["Precision"]["ALL"]["Precision_edit_dist"] = 0.9
    expect["Precision"]["FILT"]["Precision_edit_dist"] = 0.88888889
    expect["Recall"]["ALL"]["Recall_edit_dist"] = 0.75
    expect["Recall"]["FILT"]["Recall_edit_dist"] = 0.93333333
    assert stats == expect


def test_evaluate_vcf():
    ref_mask_bed_file = os.path.join(data_dir, "evaluate_vcf.ref_mask.bed")
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
        vcf_to_eval,
        ref_fasta,
        truth_fasta,
        100,
        tmp_out,
        debug=True,
        force=True,
        ref_mask_bed_file=ref_mask_bed_file,
    )
    summary_stats_expect_json = os.path.join(
        data_dir, "evaluate_vcf.expect.masked.summary_stats.json"
    )
    summary_stats_got_json = os.path.join(tmp_out, "summary_stats.json")
    assert filecmp.cmp(summary_stats_got_json, summary_stats_expect_json, shallow=False)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
