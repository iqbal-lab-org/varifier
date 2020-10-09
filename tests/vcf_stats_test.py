import os
import pytest

from cluster_vcf_records import vcf_record

from varifier import vcf_stats

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_stats")


def test_frs_from_vcf_record():
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tFRS\t0.42")
    assert vcf_stats._frs_from_vcf_record(record) == 0.42
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tFRS\t.")
    assert vcf_stats._frs_from_vcf_record(record) == "NA"
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/0")
    assert vcf_stats._frs_from_vcf_record(record) == "NA"
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tGT:COV\t0/0:1,3")
    assert vcf_stats._frs_from_vcf_record(record) == 0.25
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tGT:COV\t./.:1,3")
    assert vcf_stats._frs_from_vcf_record(record) == "NA"
    record = vcf_record.VcfRecord("ref\t1\t.\tA\tT\t.\tPASS\t.\tGT:COV\t0/0:0,0")
    assert vcf_stats._frs_from_vcf_record(record) == 0


def test_format_dict_to_edit_dist_scores():
    format_dict = {
        "VFR_RESULT": "CANNOT_USE_GT",
    }
    assert (None, None) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": 1,
        "VFR_ED_TA": 0,
    }
    assert (1, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 0,
    }
    assert (1, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 6,
        "VFR_ED_TA": 5,
    }
    assert (1, 6) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "FP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 0,
    }
    assert (0, 2) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "FP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 6,
        "VFR_ED_TA": 5,
    }
    assert (1, 6) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 1,
    }
    assert (0, 2) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": "NA",
        "VFR_ED_TA": 5,
        "VFR_ALLELE_MATCH_FRAC": 0.75,
    }
    assert (0.75, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)


def test_per_record_stats_from_vcf_file():
    infile = os.path.join(data_dir, "per_record_stats_from_vcf_file.vcf")
    expect = [
        {
            "CHROM": "ref1",
            "DP": 10,
            "DPF": 1.1,
            "FRS": 0.6,
            "GT_CONF": 100.0,
            "GT_CONF_PERCENTILE": 0.12,
            "POS": 1,
            "VFR_IN_MASK": 1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_RESULT": "FP",
        },
        {
            "CHROM": "ref2",
            "DP": 70,
            "DPF": 0.97,
            "FRS": 0.99,
            "GT_CONF": 200.0,
            "GT_CONF_PERCENTILE": 0.95,
            "POS": 2,
            "VFR_IN_MASK": 0,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_RESULT": "FP",
        },
    ]
    got = vcf_stats.per_record_stats_from_vcf_file(infile)
    assert got == expect


def test_summary_stats_from_per_record_stats():
    record_stats = [
        {
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.12,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.2,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_RESULT": "Partial_TP",
            "VFR_ALLELE_MATCH_FRAC": 0.25,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 2,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.3,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.4,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 3,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.5,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_IN_MASK": 1,
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
    ]
    expect = {
        "TP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.72, "SUM_EDIT_DIST": 3},
        "FP": {"Count": 4, "SUM_ALLELE_MATCH_FRAC": 1.15, "SUM_EDIT_DIST": 4},
        "EDIT_DIST_COUNTS": {"numerator": 6, "denominator": 10},
        "UNUSED": {"CONFLICT": 0, "MASKED": 1, "OTHER": 0},
    }
    got = vcf_stats.summary_stats_from_per_record_stats(record_stats)
    assert got == expect

    got = vcf_stats.summary_stats_from_per_record_stats(record_stats, for_recall=True)
    expect["FN"] = expect["FP"]
    del expect["FP"]
    assert got == expect
