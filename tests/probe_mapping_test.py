import filecmp
import os
import pytest

from varifier import vcf_qc_annotate, probe_mapping

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "probe_mapping")


def clean_files(filenames):
    for filename in filenames:
        if os.path.exists(filename):
            os.unlink(filename)


def test_annotate_vcf_with_probe_mapping():
    # This is an end-to-end test of running annotate_vcf_with_probe_mapping().
    # Input files are made by the script tests/data/probe_mapping/make_test_data.py.
    # It makes a VCF file + matching ref FASTA, and a truth reference FASTA.
    #
    # The aim was to make this as comprehensive as reasonably possible.
    # Tests calling TPs correctly, and also calling FPs correctly - particularly
    # in the positions flanking the true variants in case of off-by-one errors
    # or the minimap2/mappy mapping doing unexpected things.
    #
    # Also, test reverse-complementing the truth genome results in exactly the
    # same output. Important to test because probes then all map to the reverse
    # strand, which is a potential source of bugs.
    vcf_ref_fa = os.path.join(data_dir, "annotate_vcf_with_probe_mapping.ref.fa")
    vcf_in = os.path.join(data_dir, "annotate_vcf_with_probe_mapping.in.vcf")
    truth_ref_fa = os.path.join(data_dir, "annotate_vcf_with_probe_mapping.truth.fa")
    truth_ref_revcomp_fa = os.path.join(
        data_dir, "annotate_vcf_with_probe_mapping.truth.revcomp.fa"
    )
    tmp_vcf = "tmp.probe_mapping.annotate_vcf_with_probe_mapping.vcf"
    tmp_vcf_revcomp = f"{tmp_vcf}.revcomp"
    tmp_map = "tmp.probe_mapping.annotate_vcf_with_probe_mapping.map"
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_fa,
        100,
        tmp_vcf,
        map_outfile=tmp_map,
        use_fail_conflict=True,
    )
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_revcomp_fa,
        100,
        tmp_vcf_revcomp,
        use_fail_conflict=True,
    )
    expect_file = os.path.join(data_dir, "annotate_vcf_with_probe_mapping.expect.vcf")
    assert filecmp.cmp(tmp_vcf, expect_file, shallow=False)
    assert filecmp.cmp(tmp_vcf_revcomp, expect_file, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
