import filecmp
import os
import pytest

from varifier import probe_mapping

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
    # There is one variant that is slightly different because of how minimap
    # aligns the ref to the probe. Can't do anything about this, it's just
    # how alignemnts work.
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
    truth_mask = {"truth":{80, 81, 82}}
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_fa,
        100,
        tmp_vcf,
        map_outfile=tmp_map,
        use_fail_conflict=True,
        truth_mask=truth_mask,
    )
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_revcomp_fa,
        100,
        tmp_vcf_revcomp,
        use_fail_conflict=True,
        truth_mask=truth_mask,
    )
    expect_vcf = os.path.join(data_dir, "annotate_vcf_with_probe_mapping.expect.vcf")
    expect_rev_vcf = os.path.join(
        data_dir, "annotate_vcf_with_probe_mapping.expect.rev.vcf"
    )
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    assert filecmp.cmp(tmp_vcf_revcomp, expect_rev_vcf, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
