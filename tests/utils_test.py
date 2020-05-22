import filecmp
import os
import pytest
import subprocess

from cluster_vcf_records import vcf_record

from varifier import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_mask_vcf_file():
    vcf_in = os.path.join(data_dir, "mask_vcf_file.in.vcf")
    vcf_expect = os.path.join(data_dir, "mask_vcf_file.expect.vcf")
    mask_bed_file = os.path.join(data_dir, "mask_vcf_file.in.bed")
    tmp_out = "tmp.mask_vcf_file.out.vcf"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    utils.mask_vcf_file(vcf_in, mask_bed_file, tmp_out)
    assert filecmp.cmp(tmp_out, vcf_expect, shallow=False)
    os.unlink(tmp_out)
