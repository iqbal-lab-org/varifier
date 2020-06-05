import filecmp
import os
import pytest

from cluster_vcf_records import vcf_record

from varifier import vcf_qc_annotate

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_qc_annotate")


def test_add_vfr_filter_to_record():
    record = vcf_record.VcfRecord("ref\t42\t.\tT\t.\t.\tPASS\t.\tGT\t0/0")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "NO_ALTS"

    record = vcf_record.VcfRecord(
        "ref\t42\t.\tT\tA\t.\tMISMAPPED_UNPLACEABLE\t.\tGT\t0/0"
    )
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "MISMAPPED_UNPLACEABLE"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tFOO\tBAR")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "NO_GT"
    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\t\t")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "NO_GT"

    record = vcf_record.VcfRecord("ref\t42\t.\t.\tA\t.\tPASS\t.\tGT\t1/1")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "NO_REF_SEQ"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tGT\t0/1")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "CANNOT_USE_GT"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA,*\t.\tPASS\t.\tGT\t1/1")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "PASS"
    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA,*\t.\tPASS\t.\tGT\t2/2")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "CANNOT_USE_GT"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tGT\t0/0")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "CANNOT_USE_GT"
    vcf_qc_annotate._add_vfr_filter_to_record(record, want_ref_calls=True)
    assert record.FORMAT["VFR_FILTER"] == "PASS"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tGT\t1/1")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "PASS"

    record = vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tFAIL_FILTER\t.\tGT\t1/1")
    vcf_qc_annotate._add_vfr_filter_to_record(record)
    assert record.FORMAT["VFR_FILTER"] == "FAIL_BUT_TEST"


def test_fix_cluster_filter_tag():
    cluster = [vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tPASS")]
    vcf_qc_annotate._fix_cluster_filter_tag(cluster)
    assert len(cluster) == 1
    assert cluster[0].FORMAT["VFR_FILTER"] == "PASS"

    cluster = [
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tFAIL_BUT_TEST")
    ]
    vcf_qc_annotate._fix_cluster_filter_tag(cluster)
    assert len(cluster) == 1
    assert cluster[0].FORMAT["VFR_FILTER"] == "FAIL_BUT_TEST"

    cluster = [
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tPASS"),
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tFAIL_BUT_TEST"),
    ]
    vcf_qc_annotate._fix_cluster_filter_tag(cluster)
    assert len(cluster) == 2
    assert cluster[0].FORMAT["VFR_FILTER"] == "PASS"
    assert cluster[1].FORMAT["VFR_FILTER"] == "FAIL_CONFLICT"

    cluster = [
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tFAIL_BUT_TEST"),
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tFAIL_BUT_TEST"),
    ]
    vcf_qc_annotate._fix_cluster_filter_tag(cluster)
    assert len(cluster) == 2
    assert cluster[0].FORMAT["VFR_FILTER"] == "FAIL_CONFLICT"
    assert cluster[1].FORMAT["VFR_FILTER"] == "FAIL_CONFLICT"

    cluster = [
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tPASS"),
        vcf_record.VcfRecord("ref\t42\t.\tT\tA\t.\tPASS\t.\tVFR_FILTER\tPASS"),
    ]
    vcf_qc_annotate._fix_cluster_filter_tag(cluster)
    assert len(cluster) == 2
    assert cluster[0].FORMAT["VFR_FILTER"] == "FAIL_CONFLICT"
    assert cluster[1].FORMAT["VFR_FILTER"] == "FAIL_CONFLICT"


def test_annotate_sorted_list_of_records():
    records = [
        vcf_record.VcfRecord("ref\t1\t.\tTA\tA\t.\tPASS\t.\tGT\t1/1"),
        vcf_record.VcfRecord("ref\t1\t.\tT\tA\t.\tPASS\t.\tGT\t1/1"),
        vcf_record.VcfRecord("ref\t5\t.\tG\tA\t.\tPASS\t.\tGT\t1/1"),
        vcf_record.VcfRecord("ref\t5\t.\tG\tT\t.\tPASS\t.\tFOO\tBAR"),
        vcf_record.VcfRecord("ref\t5\t.\tG\tC\t.\tPASS\t.\tGT\t0/1"),
        vcf_record.VcfRecord("ref\t10\t.\tTA\tA\t.\tFAIL\t.\tGT\t1/1"),
        vcf_record.VcfRecord("ref\t11\t.\tA\tC\t.\tPASS\t.\tGT\t1/1"),
    ]
    vcf_qc_annotate._annotate_sorted_list_of_records(records)
    expect = [
        vcf_record.VcfRecord(
            "ref\t1\t.\tTA\tA\t.\tPASS\t.\tGT:VFR_FILTER\t1/1:FAIL_CONFLICT"
        ),
        vcf_record.VcfRecord(
            "ref\t1\t.\tT\tA\t.\tPASS\t.\tGT:VFR_FILTER\t1/1:FAIL_CONFLICT"
        ),
        vcf_record.VcfRecord("ref\t5\t.\tG\tA\t.\tPASS\t.\tGT:VFR_FILTER\t1/1:PASS"),
        vcf_record.VcfRecord("ref\t5\t.\tG\tT\t.\tPASS\t.\tFOO:VFR_FILTER\tBAR:NO_GT"),
        vcf_record.VcfRecord(
            "ref\t5\t.\tG\tC\t.\tPASS\t.\tGT:VFR_FILTER\t0/1:CANNOT_USE_GT"
        ),
        vcf_record.VcfRecord(
            "ref\t10\t.\tTA\tA\t.\tFAIL\t.\tGT:VFR_FILTER\t1/1:FAIL_CONFLICT"
        ),
        vcf_record.VcfRecord("ref\t11\t.\tA\tC\t.\tPASS\t.\tGT:VFR_FILTER\t1/1:PASS"),
    ]
    assert records == expect


def test_add_qc_to_vcf():
    infile = os.path.join(data_dir, "add_qc_to_vcf.in.vcf")
    expect_want_ref = os.path.join(data_dir, "add_qc_to_vcf.expect.want_ref_calls.vcf")
    expect_not_want_ref = os.path.join(
        data_dir, "add_qc_to_vcf.expect.not_want_ref_calls.vcf"
    )
    outfile = "tmp.add_qc_to_vcf.out.vcf"
    if os.path.exists(outfile):
        os.unlink(outfile)
    vcf_qc_annotate.add_qc_to_vcf(infile, outfile, want_ref_calls=True)
    assert filecmp.cmp(outfile, expect_want_ref, shallow=False)
    os.unlink(outfile)
    vcf_qc_annotate.add_qc_to_vcf(infile, outfile, want_ref_calls=False)
    assert filecmp.cmp(outfile, expect_not_want_ref, shallow=False)
    os.unlink(outfile)
