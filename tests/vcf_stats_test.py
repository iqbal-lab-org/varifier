import filecmp
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
            "VFR_EDIT_DIST": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_FILTER": "PASS",
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
            "VFR_EDIT_DIST": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
        },
    ]
    got = vcf_stats.per_record_stats_from_vcf_file(infile)
    assert got == expect


def test_summary_stats_from_per_record_stats():
    record_stats = [
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.1,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.12,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.2,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "Partial_TP",
            "VFR_ALLELE_MATCH_FRAC": 0.25,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.3,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "FAIL",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.4,
            "VFR_EDIT_DIST": 1,
        },
        {
            "VFR_FILTER": "FAIL_BUT_TEST",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.5,
            "VFR_EDIT_DIST": 1,
        },
    ]
    expect = {
        "UNUSED": 1,
        "ALL": {
            "TP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.72, "SUM_EDIT_DIST": 3},
            "FP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.75, "SUM_EDIT_DIST": 3},
        },
        "FILT": {
            "TP": {"Count": 2, "SUM_ALLELE_MATCH_FRAC": 0.22, "SUM_EDIT_DIST": 2},
            "FP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.75, "SUM_EDIT_DIST": 3},
        },
    }
    got = vcf_stats.summary_stats_from_per_record_stats(record_stats)
    assert got == expect
