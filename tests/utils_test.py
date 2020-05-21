import filecmp
import os
import pytest
import subprocess

from cluster_vcf_records import vcf_record

from varifier import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_read_vcf_file():
    vcf = os.path.join(data_dir, "read_vcf_file.vcf")
    vcf_reader = utils.read_vcf_file(vcf)
    header_lines = next(vcf_reader)
    got_records = list(vcf_reader)
    expect_header = [
        "#header line 1",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    ]
    expect_records = [
        (vcf_record.VcfRecord("ref\t42\t1\tT\tA\t.\tPASS\t.\tGT\t1/1"), "PASS"),
        (
            vcf_record.VcfRecord("ref\t43\t2\tT\tA\t.\tFILTER_X\t.\tGT\t1/1"),
            "FAIL_BUT_TEST",
        ),
        (vcf_record.VcfRecord("ref\t44\t3\tT\tA\t.\tPASS\t.\tFOO\tBAR"), "NO_GT"),
        (vcf_record.VcfRecord("ref\t45\t4\t.\tA\t.\tPASS\t.\tGT\t1/1"), "NO_REF_SEQ"),
        (
            vcf_record.VcfRecord(
                "ref\t46\t5\tT\tA\t.\tMISMAPPED_UNPLACEABLE\t.\tGT\t1/1"
            ),
            "MISMAPPED_UNPLACEABLE",
        ),
        (
            vcf_record.VcfRecord("ref\t47\t6\tT\tA\t.\tPASS\t.\tGT\t0/1"),
            "CANNOT_USE_GT",
        ),
        (
            vcf_record.VcfRecord("ref\t48\t7\tT\tA\t.\tPASS\t.\tGT\t0/0"),
            "CANNOT_USE_GT",
        ),
        (vcf_record.VcfRecord("ref\t49\t8\tT\tA,*\t.\tPASS\t.\tGT\t1/1"), "PASS",),
        (
            vcf_record.VcfRecord("ref\t50\t9\tT\tA,*\t.\tPASS\t.\tGT\t2/2"),
            "CANNOT_USE_GT",
        ),
    ]

    assert header_lines == expect_header
    assert got_records == expect_records


def test_mask_vcf_file():
    vcf_in = os.path.join(data_dir, "mask_vcf_file.in.vcf")
    vcf_expect = os.path.join(data_dir, "mask_vcf_file.expect.vcf")
    mask_bed_file = os.path.join(data_dir, "mask_vcf_file.in.bed")
    tmp_out = "tmp.mask_vcf_file.out.vcf"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    utils.mask_vcf_file(vcf_in, mask_bed_file, tmp_out)
    assert filecmp.cmp(tmp_out, vcf_expect, shallow=False)
    os.unlink(tmp_out)
