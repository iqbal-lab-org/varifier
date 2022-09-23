import filecmp
import os
import pytest
import random
import subprocess

from varifier import truth_variant_finding, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "truth_variant_finding")


def test_check_dependencies_in_path():
    truth_variant_finding._check_dependencies_in_path()


def test_truth_using_minimap2_paftools():
    ref_fasta = os.path.join(data_dir, "truth_using_minimap2_paftools.ref.fa")
    truth_fasta = os.path.join(data_dir, "truth_using_minimap2_paftools.truth.fa")
    truth_revcomp_fasta = os.path.join(
        data_dir, "truth_using_minimap2_paftools.truth.revcomp.fa"
    )

    tmp_vcf = "tmp.truth_using_minimap2_paftools.vcf"
    subprocess.check_output(f"rm -f {tmp_vcf}", shell=True)
    truth_variant_finding._truth_using_minimap2_paftools(
        ref_fasta, truth_fasta, tmp_vcf
    )
    expect_vcf = os.path.join(data_dir, "truth_using_minimap2_paftools.expect.vcf")
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    os.unlink(tmp_vcf)

    # The result should be the same if we reverse complement the truth reference.
    # However, these tags: "QSTART=122;QSTRAND=-" mean that the VCF file is not
    # identical (although the CHROM, REF, ALT columns are). Hence we have a
    # different expceted VCF file
    tmp_vcf_revcomp = "tmp.truth_using_minimap2_paftools.revcomp.vcf"
    truth_variant_finding._truth_using_minimap2_paftools(
        ref_fasta, truth_revcomp_fasta, tmp_vcf_revcomp
    )
    expect_vcf_revcomp = os.path.join(
        data_dir, "truth_using_minimap2_paftools.expect.revcomp.vcf"
    )
    assert filecmp.cmp(tmp_vcf_revcomp, expect_vcf_revcomp, shallow=False)
    os.unlink(tmp_vcf_revcomp)


def test_merge_vcf_files_for_probe_mapping():
    vcf_files = [
        os.path.join(data_dir, "merge_vcf_files_for_probe_mapping.in.1.vcf"),
        os.path.join(data_dir, "merge_vcf_files_for_probe_mapping.in.2.vcf"),
    ]
    tmp_vcf = "tmp.merge_vcf_files_for_probe_mapping.vcf"
    subprocess.check_output(f"rm -rf {tmp_vcf}*", shell=True)
    ref_fasta = os.path.join(data_dir, "merge_vcf_files_for_probe_mapping.ref.fa")
    truth_variant_finding._merge_vcf_files_for_probe_mapping(
        vcf_files, ref_fasta, tmp_vcf
    )
    expect_vcf = os.path.join(data_dir, "merge_vcf_files_for_probe_mapping.expect.vcf")
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    subprocess.check_output(f"rm -r {tmp_vcf}*", shell=True)


def test_filter_fps_and_long_vars_from_probe_mapped_vcf():
    vcf_in = os.path.join(
        data_dir, "filter_fps_and_long_vars_from_probe_mapped_vcf.in.vcf"
    )
    expect_vcf = os.path.join(
        data_dir, "filter_fps_and_long_vars_from_probe_mapped_vcf.expect.vcf"
    )
    tmp_vcf = "tmp.filter_fps_and_long_vars_from_probe_mapped_vcf.vcf"
    subprocess.check_output(f"rm -f {tmp_vcf}", shell=True)
    truth_variant_finding._filter_fps_and_long_vars_from_probe_mapped_vcf(
        vcf_in, tmp_vcf, None
    )
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    os.unlink(tmp_vcf)
    truth_variant_finding._filter_fps_and_long_vars_from_probe_mapped_vcf(
        vcf_in, tmp_vcf, 1
    )
    expect_vcf = os.path.join(
        data_dir, "filter_fps_and_long_vars_from_probe_mapped_vcf.expect.max_ref_1.vcf"
    )
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    os.unlink(tmp_vcf)


def test_bcftools_norm():
    ref_fasta = os.path.join(data_dir, "bcftools_norm.ref.fa")
    vcf_in = os.path.join(data_dir, "bcftools_norm.in.vcf")
    expect_vcf = os.path.join(data_dir, "bcftools_norm.expect.vcf")
    tmp_vcf = "tmp.bcftools_norm.vcf"
    subprocess.check_output(f"rm -f {tmp_vcf}", shell=True)
    truth_variant_finding._bcftools_norm(ref_fasta, vcf_in, tmp_vcf)
    # Header lines get changed, so just check the variant lines
    assert utils.vcf_records_are_the_same(tmp_vcf, expect_vcf)
    os.unlink(tmp_vcf)


def test_make_truth_vcf():
    ref_fasta = os.path.join(data_dir, "make_truth_vcf.ref.fa")
    truth_fasta = os.path.join(data_dir, "make_truth_vcf.truth.fa")
    tmp_out = "tmp.truth_variant_finding.make_truth_ref"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(ref_fasta, truth_fasta, tmp_out, 100)
    expect_vcf = os.path.join(data_dir, "make_truth_vcf.expect.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
    # Test same run again, but mask a position in the truth where there's a SNP
    truth_mask = {"truth": {59}}
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, truth_mask=truth_mask
    )
    expect_vcf = os.path.join(data_dir, "make_truth_vcf.expect.masked.vcf")
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)


# An edge case seen in Covid data. MUMmer/dnadiff etc can still call
# and indel, even though one of the genomes has an N. eg something like:
# ref:    ACGTGCAACGTA...
# truth:  ACGTGCNNNNNN...
# and it calls CAA -> C (ie delete the first two Ns!).
# Is a special case that needs fixing, because (before writing this test and
# fixing the code), the alt probe maps OK, in the sense that it maps to the
# truth genome with a hit that includes the complete allele. But the ref probe
# does not map because it has the CA in there, which can't align to the Ns
# in the truth genome. Result is alt probe looks ok, ref probe does not, and
# then is called as a TP.
#
# Update 2021-09: there is a new global alignment method, instead of using
# the default minimap2/mummer to get variants. Changed this test to run
# using the use_global_align=True option as well as the original method.
# Is easy way of checking that the new use_global_align option works.
# New method needs to have no rearrangements, and both genomes on the same
# strand, so we're not testing the rev complement tests here because they
# would fail.
def test_make_truth_vcf_handle_Ns():
    ref_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns.ref.fa")
    truth_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns.truth.fa")
    tmp_out = "tmp.truth_variant_finding.make_truth_ref_handle_Ns"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(ref_fasta, truth_fasta, tmp_out, 100)
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns.ref_v_truth_expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, use_global_align=True
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    ref_fasta_revcomp = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns.ref.revcomp.fa"
    )
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta_revcomp, truth_fasta, tmp_out, 100
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns.ref_revcomp_v_truth.expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    got_vcf = truth_variant_finding.make_truth_vcf(truth_fasta, ref_fasta, tmp_out, 100)
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns.truth_v_ref.expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(
        truth_fasta, ref_fasta, tmp_out, 100, use_global_align=True, debug=True
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    got_vcf = truth_variant_finding.make_truth_vcf(
        truth_fasta, ref_fasta_revcomp, tmp_out, 100
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns.truth_v_ref_revcomp.expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)


# Another edge case seen in covid data with Ns causing a FP indel.
# As in previous test, run make_truth_vcf both with and without the
# use_global_align option.
def test_make_truth_vcf_handle_Ns_2():
    ref_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns_2.ref.fa")
    truth_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns_2.truth.fa")
    tmp_out = "tmp.truth_variant_finding.make_truth_ref_handle_Ns_2"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(ref_fasta, truth_fasta, tmp_out, 100)
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns_2.ref_v_truth_expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    # When using global align, indels are moved to the left compared to
    # default, so different expected VCF
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, use_global_align=True
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns_2.ref_v_truth_expect.global.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)


# This is another edge case seen in covid data. Is a SNP near to a run of
# Ns. This was motivation for making the use_global_align option. Can't see
# a way to fix minimap2/nucmer method to find the SNP (they don't), but doing
# a global alignment does find it. End of covid ref genome is:
# GCTATCCCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# End of test genome is:
# GCTAACCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
#     ^
#     | SNP here, gets missed by minimap2/nucmer
def test_make_truth_vcf_hangle_Ns_3():
    ref_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns_3.ref.fa")
    truth_fasta = os.path.join(data_dir, "make_truth_vcf_handle_Ns_3.truth.fa")
    tmp_out = "tmp.truth_variant_finding.make_truth_ref_handle_Ns_3"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, use_global_align=True
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_handle_Ns_3.ref_v_truth_expect.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    # For completeness, run the same test but without using global align.
    # We'll get a VCF file with no variants, so just header lines starting with
    # a "#"
    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, use_global_align=False
    )
    with open(got_vcf) as f:
        for line in f:
            assert line.startswith("#")
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)


def test_make_truth_vcf_homopolymer_fix():
    random.seed(42)
    seq1 = "".join(random.choices(["A", "C", "G", "T"], k=100))
    seq2 = "".join(random.choices(["A", "C", "G", "T"], k=100))
    hp = "A" * 10
    ref_seq = seq1 + "C" + hp + "AT" + seq2
    truth_seq = seq1 + "C" + hp + "T" + seq2
    ref_fasta = "tmp.make_truth_vcf_homopolymer_fix.ref.fa"
    truth_fasta = "tmp.make_truth_vcf_homopolymer_fix.truth.fa"
    tmp_out = "tmp.make_truth_vcf_homopolymer_fix"
    with open(ref_fasta, "w") as f:
        print(">ref", ref_seq, sep="\n", file=f)
    with open(truth_fasta, "w") as f:
        print(">truth", truth_seq, sep="\n", file=f)
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    got_fasta = os.path.join(tmp_out, "04.qry_sanitised_gaps.fa")

    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta, truth_fasta, tmp_out, 100, use_global_align=True
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_homopolymer_fix.ref_v_truth_expect.default.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    got_seq = utils.load_one_seq_fasta_file(got_fasta)
    assert got_seq.seq == truth_seq
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta,
        truth_fasta,
        tmp_out,
        100,
        use_global_align=True,
        hp_min_fix_length=12,
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    got_seq = utils.load_one_seq_fasta_file(got_fasta)
    assert got_seq.seq == truth_seq
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    got_vcf = truth_variant_finding.make_truth_vcf(
        ref_fasta,
        truth_fasta,
        tmp_out,
        100,
        use_global_align=True,
        hp_min_fix_length=11,
    )
    expect_vcf = os.path.join(
        data_dir, "make_truth_vcf_homopolymer_fix.ref_v_truth_expect.length_11.vcf"
    )
    assert utils.vcf_records_are_the_same(got_vcf, expect_vcf)
    got_seq = utils.load_one_seq_fasta_file(got_fasta)
    assert got_seq.seq == ref_seq
    subprocess.check_output(f"rm -r {tmp_out}", shell=True)

    os.unlink(ref_fasta)
    os.unlink(truth_fasta)
