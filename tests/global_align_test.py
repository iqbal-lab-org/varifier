import filecmp
import os
import pytest
import subprocess
from unittest import mock

from varifier import global_align, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "global_align")


def test_get_perfect_matches():
    tmp_nucmer = "tmp.get_perfect_matches"
    subprocess.check_output("rm -f {tmp_nucmer}", shell=True)
    ref_fasta = os.path.join(data_dir, "get_perfect_matches.ref.fa")
    qry_fasta = os.path.join(data_dir, "get_perfect_matches.qry.fa")
    matches = global_align.get_perfect_matches(ref_fasta, qry_fasta, tmp_nucmer)
    assert len(matches) == 2
    assert matches[0].ref_start == matches[0].qry_start == 0
    assert matches[0].ref_end == matches[0].qry_end == 59
    assert matches[1].ref_start == matches[1].qry_start == 61
    assert matches[1].ref_end == matches[1].qry_end == 99
    assert matches[0].percent_identity == matches[1].percent_identity == 100
    assert not os.path.exists(tmp_nucmer)


def test_perfect_matches_to_conservative_match_coords():
    nuc_match0 = mock.Mock()
    nuc_match0.ref_start = 30
    nuc_match0.ref_end = 40
    nuc_match0.qry_start = 30
    nuc_match0.qry_end = 40
    nuc_match0.ref_length = 100
    nuc_match0.qry_length = 100
    with pytest.raises(Exception):
        global_align.perfect_matches_to_conservative_match_coords([nuc_match0])

    nuc_match1 = mock.Mock()
    nuc_match1.ref_start = 0
    nuc_match1.ref_end = 42
    nuc_match1.qry_start = 0
    nuc_match1.qry_end = 42
    nuc_match1.ref_length = 100
    nuc_match1.qry_length = 100
    expect_match1 = {"ref_start": 0, "ref_end": 37, "qry_start": 0, "qry_end": 37}
    got = global_align.perfect_matches_to_conservative_match_coords(
        [nuc_match0, nuc_match1]
    )
    assert got == [expect_match1]

    nuc_match2 = mock.Mock()
    nuc_match2.ref_start = 80
    nuc_match2.ref_end = 99
    nuc_match2.qry_start = 80
    nuc_match2.qry_end = 99
    nuc_match2.ref_length = 100
    nuc_match2.qry_length = 100
    expect_match2 = {"ref_start": 85, "ref_end": 99, "qry_start": 85, "qry_end": 99}
    got = global_align.perfect_matches_to_conservative_match_coords(
        [nuc_match1, nuc_match2]
    )
    assert got == [expect_match1, expect_match2]

    nuc_match1.ref_start = 1
    nuc_match1.ref_end = 43
    expect_match1 = {"ref_start": 6, "ref_end": 38, "qry_start": 5, "qry_end": 37}
    got = global_align.perfect_matches_to_conservative_match_coords([nuc_match1])
    assert got == [expect_match1]

    nuc_match2.ref_start = 79
    nuc_match2.ref_end = 98
    expect_match2 = {"ref_start": 84, "ref_end": 93, "qry_start": 85, "qry_end": 94}
    got = global_align.perfect_matches_to_conservative_match_coords([nuc_match2])
    assert got == [expect_match2]


# this is based on a real example. It is possible to have an indel resulting in
# perfect matches that overlap slightly. We just remove one of the overlapping
# matches because is filled in by local alignment afterwards anyway
def test_perfect_matches_to_conservative_match_coords_handle_qry_overlap():
    # Test only keep only first match when ref coords overlap
    nuc_match0 = mock.Mock()
    nuc_match0.ref_start = 0
    nuc_match0.ref_end = 42
    nuc_match0.qry_start = 0
    nuc_match0.qry_end = 42
    nuc_match0.ref_length = 100
    nuc_match0.qry_length = 100

    nuc_match1 = mock.Mock()
    nuc_match1.ref_start = 30
    nuc_match1.ref_end = 60
    nuc_match1.qry_start = 50
    nuc_match1.qry_end = 80
    nuc_match1.ref_length = 100
    nuc_match1.qry_length = 100

    expect_match0 = {"ref_start": 0, "ref_end": 37, "qry_start": 0, "qry_end": 37}
    got = global_align.perfect_matches_to_conservative_match_coords(
        [nuc_match0, nuc_match1]
    )
    assert got == [expect_match0]

    # Test only keep only first match when qry coords overlap
    nuc_match0 = mock.Mock()
    nuc_match0.ref_start = 0
    nuc_match0.ref_end = 42
    nuc_match0.qry_start = 0
    nuc_match0.qry_end = 42
    nuc_match0.ref_length = 100
    nuc_match0.qry_length = 100

    nuc_match1 = mock.Mock()
    nuc_match1.ref_start = 60
    nuc_match1.ref_end = 90
    nuc_match1.qry_start = 33
    nuc_match1.qry_end = 72
    nuc_match1.ref_length = 100
    nuc_match1.qry_length = 100

    expect_match0 = {"ref_start": 0, "ref_end": 37, "qry_start": 0, "qry_end": 37}
    got = global_align.perfect_matches_to_conservative_match_coords(
        [nuc_match0, nuc_match1]
    )
    assert got == [expect_match0]


def test_fix_query_gaps_in_msa():
    ref = list("ACG")
    qry = list("ATG")
    got_ref, got_qry = global_align.fix_query_gaps_in_msa(ref, qry)
    assert got_ref == ref
    assert got_qry == qry

    ref = list("A-G")
    qry = list("ANC")
    got_ref, got_qry = global_align.fix_query_gaps_in_msa(ref, qry)
    assert got_ref == list("AG")
    assert got_qry == list("AC")

    ref = list("A-CG")
    qry = list("ANNC")
    got_ref, got_qry = global_align.fix_query_gaps_in_msa(ref, qry)
    assert got_ref == list("ACG")
    assert got_qry == list("ANC")

    ref = list("ACGCTACGT")
    qry = list("-NG-N-N--")
    got_ref, got_qry = global_align.fix_query_gaps_in_msa(ref, qry)
    assert got_ref == ref
    assert got_qry == list("NNGNNNNNN")

    ref = list("A--T-ACTGC-GTTTGAGTAGT")
    qry = list("ACGTNN-NGCAG-T--NN-N-T")
    got_ref, got_qry = global_align.fix_query_gaps_in_msa(ref, qry)
    assert got_ref == list("A--TACTGC-GTTTGAGTAGT")
    assert got_qry == list("ACGTNNNGCAG-TNNNNNNNT")


def test_global_align():
    ref_fasta = os.path.join(data_dir, "global_aln.ref.fa")
    qry_fasta = os.path.join(data_dir, "global_aln.qry.fa")
    tmp_nucmer = "tmp.gloabl_align.nucmer"
    subprocess.check_output(f"rm -rf {tmp_nucmer}", shell=True)
    got = global_align.global_align(ref_fasta, qry_fasta, tmp_nucmer)
    expect = (
        "AGCCCCGAGGCTATTCGT-CACCGTAGTGCGTCGACTCCCGATAGTCCT-ATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCG---",
        "AG-CCCGAGCCTATTCGNNCACCGTAGTGCGTCGACTCCCGATAGTCCTCATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCGGGT",
    )
    assert got == expect

    got = global_align.global_align(
        ref_fasta, qry_fasta, tmp_nucmer, fix_query_gap_lengths=True
    )
    expect = (
        "AGCCCCGAGGCTATTCGTCACCGTAGTGCGTCGACTCCCGATAGTCCT-ATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCG---",
        "AG-CCCGAGCCTATTCGNCACCGTAGTGCGTCGACTCCCGATAGTCCTCATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCGGGT",
    )
    assert got == expect


def test_variants_from_global_alignment():
    #          01234567--89012345
    ref_aln = "AGCTGCGC--CNATCGAT-"
    #          |  |||||||||||   |
    qry_aln = "A-TTGCGCATCTATTACTA"
    got = global_align.variants_from_global_alignment(ref_aln, qry_aln)
    expect = [
        {"ref_start": 0, "ref_allele": "AGC", "qry_allele": "AT"},
        {"ref_start": 7, "ref_allele": "C", "qry_allele": "CAT"},
        {"ref_start": 12, "ref_allele": "CGA", "qry_allele": "TAC"},
        {"ref_start": 15, "ref_allele": "T", "qry_allele": "TA"},
    ]
    assert got == expect


def test_variants_from_global_alignment_insertion_at_start():
    #             0123456789012
    ref_aln = "---ACGTGTACGTACG"
    #             |||||||||||||
    qry_aln = "GCGACGTGTACGTACG"
    got = global_align.variants_from_global_alignment(ref_aln, qry_aln)
    expect = [
        {"ref_start": 0, "ref_allele": "A", "qry_allele": "GCGA"},
    ]
    assert got == expect


def test_variants_from_global_alignment_deletion_at_start():
    #          0123456789012345
    ref_aln = "GCGACGTGTACGTACG"
    #             |||||||||||||
    qry_aln = "---ACGTGTACGTACG"
    got = global_align.variants_from_global_alignment(ref_aln, qry_aln)
    expect = [
        {"ref_start": 0, "ref_allele": "GCGA", "qry_allele": "A"},
    ]
    assert got == expect


def test_expand_combined_snps():
    var1 = {"ref_start": 0, "ref_allele": "AGC", "qry_allele": "AT"}
    var2 = {"ref_start": 2, "ref_allele": "A", "qry_allele": "T"}
    var3 = {"ref_start": 42, "ref_allele": "AG", "qry_allele": "TC"}
    var3_1 = {"ref_start": 42, "ref_allele": "A", "qry_allele": "T"}
    var3_2 = {"ref_start": 43, "ref_allele": "G", "qry_allele": "C"}
    variants_in = [var1, var2, var3]
    got = global_align.expand_combined_snps(variants_in)
    expect = [var1, var2, var3_1, var3_2]
    assert got == expect


def test_vcf_using_global_alignment():
    ref_fasta = os.path.join(data_dir, "vcf_using_global_alignment.ref.fa")
    qry_fasta = os.path.join(data_dir, "vcf_using_global_alignment.qry.fa")
    expect_vcf = os.path.join(data_dir, "vcf_using_global_alignment.vcf")
    vcf_out = "tmp.vcf_using_global_alignment.vcf"
    subprocess.check_output(f"rm -rf {vcf_out}", shell=True)
    tmp_msa = "tmp.global_aln.msa"
    subprocess.check_output(f"rm -rf {tmp_msa}", shell=True)
    tmp_qry_fasta = "tmp.global_aln.qry.fa"
    subprocess.check_output(f"rm -rf {tmp_qry_fasta}", shell=True)
    global_align.vcf_using_global_alignment(
        ref_fasta,
        qry_fasta,
        vcf_out,
        fix_query_gap_lengths=True,
        fixed_query_fasta=tmp_qry_fasta,
        msa_file=tmp_msa,
    )
    assert utils.vcf_records_are_the_same(vcf_out, expect_vcf)
    expect_qry_fa = os.path.join(
        data_dir, "vcf_using_global_alignment.expect_qry_out.fa"
    )
    expect_msa = os.path.join(data_dir, "vcf_using_global_alignment.msa")
    assert filecmp.cmp(expect_qry_fa, tmp_qry_fasta, shallow=False)
    assert filecmp.cmp(expect_msa, tmp_msa, shallow=False)
    os.unlink(tmp_qry_fasta)
    os.unlink(tmp_msa)
    os.unlink(vcf_out)

    expect_vcf = os.path.join(data_dir, "vcf_using_global_alignment.1-100.vcf")
    global_align.vcf_using_global_alignment(
        ref_fasta, qry_fasta, vcf_out, min_ref_coord=0, max_ref_coord=99
    )
    assert utils.vcf_records_are_the_same(vcf_out, expect_vcf)
    os.unlink(vcf_out)

    expect_vcf = os.path.join(data_dir, "vcf_using_global_alignment.2-99.vcf")
    global_align.vcf_using_global_alignment(
        ref_fasta, qry_fasta, vcf_out, min_ref_coord=1, max_ref_coord=98
    )
    assert utils.vcf_records_are_the_same(vcf_out, expect_vcf)
    os.unlink(vcf_out)

    expect_vcf = os.path.join(data_dir, "vcf_using_global_alignment.3-98.vcf")
    global_align.vcf_using_global_alignment(
        ref_fasta, qry_fasta, vcf_out, min_ref_coord=2, max_ref_coord=97
    )
    assert utils.vcf_records_are_the_same(vcf_out, expect_vcf)
    os.unlink(vcf_out)
