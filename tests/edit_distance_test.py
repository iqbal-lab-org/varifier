import pytest

from varifier import edit_distance


def test_needleman_wunsch():
    seq1 = "ACGTGTCACAG"
    seq2 = "AGTCTGACATG"
    aln1, aln2 = edit_distance._needleman_wunsch(seq1, seq2)
    expect1 = "ACGTGTCACA-G"
    expect2 = "A-GTCTGACATG"
    assert aln1 == expect1
    assert aln2 == expect2


def test_edit_distance_from_aln_strings():
    assert edit_distance.edit_distance_from_aln_strings("A", "A") == 0
    assert edit_distance.edit_distance_from_aln_strings("A", "C") == 1
    assert edit_distance.edit_distance_from_aln_strings("A-G", "ACG") == 1
    assert edit_distance.edit_distance_from_aln_strings("A--G", "ACTG") == 1
    assert edit_distance.edit_distance_from_aln_strings("ACG", "A-G") == 1
    assert edit_distance.edit_distance_from_aln_strings("ACTG", "A--G") == 1
    assert edit_distance.edit_distance_from_aln_strings("ACTC", "A--G") == 2
    assert edit_distance.edit_distance_from_aln_strings("A--G", "ACTC") == 2
    assert edit_distance.edit_distance_from_aln_strings("A-TG", "AC-G") == 1
    assert edit_distance.edit_distance_from_aln_strings("AC-G", "A-TG") == 1
    assert edit_distance.edit_distance_from_aln_strings("A--G", "A-CG") == 1
    assert edit_distance.edit_distance_from_aln_strings("A-CG", "A--G") == 1
    assert edit_distance.edit_distance_from_aln_strings("A--TCAC-G", "AGT--ACTG") == 2
    assert edit_distance.edit_distance_from_aln_strings("AGT--ACTG", "A--TCAC-G") == 2
    assert edit_distance.edit_distance_from_aln_strings("A--TCAC-GG", "AGT--ACTG-") == 3
    assert edit_distance.edit_distance_from_aln_strings("AGT--ACTG-", "A--TCAC-GG") == 3


def test_edit_distance_between_seqs():
    edit_distance.edit_distance_between_seqs("A", "A") == 0
    edit_distance.edit_distance_between_seqs("A", "C") == 1
    edit_distance.edit_distance_between_seqs("AG", "ACG") == 1
    edit_distance.edit_distance_between_seqs("AGTGCAT", "ACGTGCGT") == 2
    edit_distance.edit_distance_between_seqs("AGTGCAT", "AGCCTGCAT") == 1
