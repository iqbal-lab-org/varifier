import collections
import pytest

from varifier import probe


def test_allele_seq():
    p = probe.Probe("aaaCGttt", 3, 4)
    assert p.allele_seq() == "CG"


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
