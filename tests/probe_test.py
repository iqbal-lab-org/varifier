import collections
import pytest

from varifier import probe


def test_allele_seq():
    p = probe.Probe("aaaCGttt", 3, 4)
    assert p.allele_seq() == "CG"


def test_allele_is_in_mask():
    p = probe.Probe("ACTGACTGACTG", 3, 6)
    assert not p.allele_is_in_mask({2})
    assert p.allele_is_in_mask({2, 3})
    assert p.allele_is_in_mask({3})
    assert p.allele_is_in_mask({4})
    assert p.allele_is_in_mask({5})
    assert p.allele_is_in_mask({6})
    assert p.allele_is_in_mask({6, 7})
    assert not p.allele_is_in_mask({7})


def test_probe_allele_match_counts():
    Hit = collections.namedtuple("Hit", ["NM", "q_st", "strand", "cigar"])
    p = probe.Probe("ACGTA", 2, 2)
    cigar = [[5, 7]]  # this is 5=
    hit = Hit(0, 0, 1, cigar)
    assert p.allele_match_counts(hit) == (1, 1)

    cigar = [[2, 7], [1, 8], [2, 7]]  # this is 2=1X2=
    hit = Hit(1, 0, 1, cigar)
    assert p.allele_match_counts(hit) == (0, 1)

    cigar = [[1, 7], [1, 8], [3, 7]]  # this is 1=1X3=
    hit = Hit(1, 0, 1, cigar)
    assert p.allele_match_counts(hit) == (1, 1)

    cigar = [[3, 7], [1, 8], [1, 7]]  # this is 3=1X1=
    hit = Hit(1, 0, 1, cigar)
    assert p.allele_match_counts(hit) == (1, 1)


def test_raise_error_bad_cigar_operator():
    Hit = collections.namedtuple("Hit", ["NM", "q_st", "strand", "cigar"])
    p = probe.Probe("ACGTA", 2, 2)
    cigar = [[5, 42]]  # 42 is not a valid cigar operator
    hit = Hit(1, 0, 1, cigar)
    with pytest.raises(RuntimeError):
        foo = p.allele_match_counts(hit)


def test_padded_probe_or_ref_seq():
    Hit = collections.namedtuple("Hit", ["r_st", "q_st", "q_en", "strand", "cigar"])
    ref = "TCCGTATGG"
    ref_rev = "CCATACGGA"
    p = probe.Probe("ACGTAT", 2, 2)

    cigar = [[5, 7]]
    hit = Hit(2, 1, 6, 1, cigar)
    assert "NCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=None)
    assert "NCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=ref)
    hit = Hit(2, 1, 6, -1, cigar)
    assert "NCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=ref_rev)

    cigar = [[6, 7]]
    hit = Hit(1, 0, 6, 1, cigar)
    assert "ACGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=None)
    assert "CCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=ref)
    hit = Hit(2, 0, 5, -1, cigar)
    assert "CCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=ref_rev)

    cigar = [[2, 7], [1, 1], [3, 7]]  # 1=1I4=
    hit = Hit(1, 0, 6, 1, cigar)
    assert "ACGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=None)
    assert "CC-GTA" == p.padded_probe_or_ref_seq(hit, ref_seq=ref)
    hit = Hit(3, 0, 6, -1, cigar)
    assert "CCG-TA" == p.padded_probe_or_ref_seq(hit, ref_seq=ref_rev)

    cigar = [[1, 7], [1, 2], [5, 7]]  # 1=1D5=
    hit = Hit(1, 0, 6, 1, cigar)
    assert "A-CGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=None)
    assert "CCGTATG" == p.padded_probe_or_ref_seq(hit, ref_seq=ref)
    hit = Hit(1, 0, 6, -1, cigar)
    assert "CCGTATG" == p.padded_probe_or_ref_seq(hit, ref_seq=ref_rev)

    cigar = [[1, 7], [1, 2], [2, 8], [2, 2], [3, 7]]  # 1=1D2X2D3=
    hit = Hit(1, 0, 7, 1, cigar)
    assert "A-CG--TAT" == p.padded_probe_or_ref_seq(hit, ref_seq=None)
    assert "CCGTATGG" == p.padded_probe_or_ref_seq(hit, ref_seq=ref)
    hit = Hit(2, 0, 7, -1, cigar)
    assert "TCCGTAT" == p.padded_probe_or_ref_seq(hit, ref_seq=ref_rev)


def test_padded_seq_allele_start_end_coords():
    p = probe.Probe("ACGTA", 2, 2)
    assert (2, 2) == p.padded_seq_allele_start_end_coords("ACGTA")
    assert (2, 2) == p.padded_seq_allele_start_end_coords("ACG-TA")
    assert (3, 3) == p.padded_seq_allele_start_end_coords("-ACGTA")
    assert (3, 3) == p.padded_seq_allele_start_end_coords("A-CGTA")
    assert (3, 3) == p.padded_seq_allele_start_end_coords("AC-GTA")
    assert (4, 4) == p.padded_seq_allele_start_end_coords("--ACGTA")
    assert (4, 4) == p.padded_seq_allele_start_end_coords("-A-CGTA")

    p = probe.Probe("ACGTATC", 2, 4)
    assert (2, 4) == p.padded_seq_allele_start_end_coords("ACGTATC")
    assert (2, 5) == p.padded_seq_allele_start_end_coords("ACG-TATC")
    assert (2, 6) == p.padded_seq_allele_start_end_coords("ACG--TATC")
    assert (2, 7) == p.padded_seq_allele_start_end_coords("ACG--T-ATC")
    assert (2, 7) == p.padded_seq_allele_start_end_coords("ACG--T-A-TC")
    assert (3, 8) == p.padded_seq_allele_start_end_coords("-ACG--T-A-TC")
