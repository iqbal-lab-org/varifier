import filecmp
import os
import pytest
import subprocess

from varifier import probe_mapping, utils

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
    truth_mask = {"truth": {80, 81, 82}}
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
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    expect_vcf_rev = os.path.join(
        data_dir, "annotate_vcf_with_probe_mapping.expect.rev.vcf"
    )
    assert filecmp.cmp(tmp_vcf_revcomp, expect_vcf_rev, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))


# Clusters of SNPs and indels are hard to evaluate when they are in separate
# records. This test is to check that it works. It was found when testing
# simulated TB data, where there was a cluster of SNPs and indels. This test
# was added when the code was changed to apply variants to the flanks of each
# probe. It is all TPs. Without adding variants to the flanks, the probe mapping
# gets messed up and about 1/3 of the calls in this test were incorrectly called
# as FP instead of TP. It's the first 2kb of the H37Rv version 3 reference
# genome with a cluster of variants at around 1000-1050
def test_annotate_with_probe_mapping_clustered_snps_and_indels():
    vcf_ref_fa = os.path.join(data_dir, "clustered_snp_indel.ref.fa")
    vcf_in = os.path.join(data_dir, "clustered_snp_indel.in.vcf")
    truth_ref_fa = os.path.join(data_dir, "clustered_snp_indel.truth.fa")
    truth_ref_revcomp_fa = os.path.join(
        data_dir, "clustered_snp_indel.truth.revcomp.fa"
    )
    tmp_vcf = "tmp.probe_mapping.clustered_snp_indel.vcf"
    tmp_vcf_revcomp = f"{tmp_vcf}.revcomp"
    tmp_map = "tmp.probe_mapping.clustered_snp_indel.map"
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_fa,
        100,
        tmp_vcf,
        map_outfile=tmp_map,
    )
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_revcomp_fa,
        100,
        tmp_vcf_revcomp,
        map_outfile=tmp_map,
    )
    expect_vcf = os.path.join(data_dir, "clustered_snp_indel.expect.vcf")
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    expect_vcf_rev = os.path.join(data_dir, "clustered_snp_indel.expect.rev.vcf")
    assert filecmp.cmp(tmp_vcf_revcomp, expect_vcf_rev, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))


# This was seen in covid data. Ref allele = ACTTTACTTG, alt=A. We need to
# check that can correctly call TP for this. But also check that
# ref allele = ACTTTACTT, alt=A gets called as FP - this was the error seen
# in covid data.
def test_deletion_wrong_lengths():
    vcf_ref_fa = os.path.join(data_dir, "deletion_wrong_lengths.ref.fa")
    vcf_in = os.path.join(data_dir, "deletion_wrong_lengths.in.vcf")
    truth_ref_fa = os.path.join(data_dir, "deletion_wrong_lengths.truth.fa")
    truth_ref_revcomp_fa = os.path.join(
        data_dir, "deletion_wrong_lengths.truth.revcomp.fa"
    )
    tmp_vcf = "tmp.probe_mapping.deletion_wrong_lengths.vcf"
    tmp_vcf_revcomp = f"{tmp_vcf}.revcomp"
    tmp_map = "tmp.probe_mapping.deletion_wrong_lengths.map"
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_fa,
        100,
        tmp_vcf,
        map_outfile=tmp_map,
    )
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_revcomp_fa,
        100,
        tmp_vcf_revcomp,
        map_outfile=tmp_map,
    )
    expect_vcf = os.path.join(data_dir, "deletion_wrong_lengths.expect.vcf")
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    expect_vcf_rev = os.path.join(data_dir, "deletion_wrong_lengths.expect.rev.vcf")
    assert filecmp.cmp(tmp_vcf_revcomp, expect_vcf_rev, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))


# Same as the previous test, but this flips the truth/ref genomes around so
# we're testing insertions instead
def test_insertion_wrong_lengths():
    vcf_ref_fa = os.path.join(data_dir, "insertion_wrong_lengths.ref.fa")
    vcf_in = os.path.join(data_dir, "insertion_wrong_lengths.in.vcf")
    truth_ref_fa = os.path.join(data_dir, "insertion_wrong_lengths.truth.fa")
    truth_ref_revcomp_fa = os.path.join(
        data_dir, "insertion_wrong_lengths.truth.revcomp.fa"
    )
    tmp_vcf = "tmp.probe_mapping.insertion_wrong_lengths.vcf"
    tmp_vcf_revcomp = f"{tmp_vcf}.revcomp"
    tmp_map = "tmp.probe_mapping.insertion_wrong_lengths.map"
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_fa,
        100,
        tmp_vcf,
        map_outfile=tmp_map,
    )
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_in,
        vcf_ref_fa,
        truth_ref_revcomp_fa,
        100,
        tmp_vcf_revcomp,
        map_outfile=tmp_map,
    )
    expect_vcf = os.path.join(data_dir, "insertion_wrong_lengths.expect.vcf")
    expect_vcf_revcomp = os.path.join(
        data_dir, "insertion_wrong_lengths.revcomp.expect.vcf"
    )
    assert filecmp.cmp(tmp_vcf, expect_vcf, shallow=False)
    assert filecmp.cmp(tmp_vcf_revcomp, expect_vcf_revcomp, shallow=False)
    clean_files((tmp_vcf, tmp_vcf_revcomp, tmp_map))


def test_qry_variant_vcf_tag():
    genome1 = os.path.join(data_dir, "qry_variant_vcf_tag.genome1.fa")
    genome2 = os.path.join(data_dir, "qry_variant_vcf_tag.genome2.fa")
    genome1_rev = os.path.join(data_dir, "qry_variant_vcf_tag.genome1.revcomp.fa")
    genome2_rev = os.path.join(data_dir, "qry_variant_vcf_tag.genome2.revcomp.fa")
    expect_1_vcf = os.path.join(data_dir, "qry_variant_vcf_tag.expect.wrt_genome1.vcf")
    expect_1_rev_vcf = os.path.join(
        data_dir, "qry_variant_vcf_tag.expect.wrt_genome1.rev.vcf"
    )
    expect_2_vcf = os.path.join(data_dir, "qry_variant_vcf_tag.expect.wrt_genome2.vcf")
    expect_2_rev_vcf = os.path.join(
        data_dir, "qry_variant_vcf_tag.expect.wrt_genome2.rev.vcf"
    )
    vcf_1_v_2 = os.path.join(data_dir, "qry_variant_vcf_tag.in.wrt_genome1.vcf")
    vcf_2_v_1 = os.path.join(data_dir, "qry_variant_vcf_tag.in.wrt_genome2.vcf")
    tmp_vcf = "tmp.probe_mapping.qry_variant_vcf_tag.wrt_genome1.vcf"
    subprocess.check_output(f"rm -f {tmp_vcf}", shell=True)
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_1_v_2,
        genome1,
        genome2,
        100,
        tmp_vcf,
    )
    assert utils.vcf_records_are_the_same(tmp_vcf, expect_1_vcf)
    os.unlink(tmp_vcf)

    tmp_vcf = "tmp.probe_mapping.qry_variant_vcf_tag.wrt_genome1.rev.vcf"
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_1_v_2,
        genome1,
        genome2_rev,
        100,
        tmp_vcf,
    )
    assert utils.vcf_records_are_the_same(tmp_vcf, expect_1_rev_vcf)
    os.unlink(tmp_vcf)

    tmp_vcf = "tmp.probe_mapping.qry_variant_vcf_tag.wrt_genome2.vcf"
    subprocess.check_output(f"rm -f {tmp_vcf}", shell=True)
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_2_v_1,
        genome2,
        genome1,
        100,
        tmp_vcf,
    )
    assert utils.vcf_records_are_the_same(tmp_vcf, expect_2_vcf)
    os.unlink(tmp_vcf)

    tmp_vcf = "tmp.probe_mapping.qry_variant_vcf_tag.wrt_genome2.rev.vcf"
    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_2_v_1,
        genome2,
        genome1_rev,
        100,
        tmp_vcf,
    )
    assert utils.vcf_records_are_the_same(tmp_vcf, expect_2_rev_vcf)
    os.unlink(tmp_vcf)
