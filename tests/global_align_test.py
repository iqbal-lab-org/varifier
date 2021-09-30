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


def test_global_align():
    ref_fasta = os.path.join(data_dir, "global_aln.ref.fa")
    qry_fasta = os.path.join(data_dir, "global_aln.qry.fa")
    tmp_nucmer = "tmp.gloabl_align.nucmer"
    subprocess.check_output(f"rm -rf {tmp_nucmer}", shell=True)
    got = global_align.global_align(ref_fasta, qry_fasta, tmp_nucmer)
    expect = (
        "AGCCCCGAGGCTATTCGTCACCGTAGTGCGTCGACTCCCGATAGTCCT-ATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCG---",
        "AG-CCCGAGCCTATTCGTCACCGTAGTGCGTCGACTCCCGATAGTCCTCATGAATTAATATTTGGGCAAAAATGATTGGAGATACCGAGTTGATGGTCCCGGGT",
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
    global_align.vcf_using_global_alignment(ref_fasta, qry_fasta, vcf_out)
    assert utils.vcf_records_are_the_same(vcf_out, expect_vcf)
    os.unlink(vcf_out)
