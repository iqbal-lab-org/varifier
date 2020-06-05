from Bio import pairwise2


def _needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap_open=-5, gap_extend=-3):
    """Returns global alignment strings from NM alignment of the
    two sequences. Dashes for gaps"""
    alignments = pairwise2.align.globalms(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )
    assert len(alignments[0][0]) == len(alignments[0][1])
    return alignments[0][0], alignments[0][1]


def edit_distance_from_aln_strings(str1, str2):
    """Input should be seqs output by _needleman_wunsch().
    Returns the edit distance between the sequences"""
    assert len(str1) == len(str2)
    edit_distance = 0
    in_gap = False

    for i, char1 in enumerate(str1):
        if char1 == "-" or str2[i] == "-":
            if not in_gap:
                in_gap = True
                edit_distance += 1
        else:
            in_gap = False

            if char1 != str2[i]:
                edit_distance += 1

    return edit_distance


def edit_distance_between_seqs(seq1, seq2):
    """Input is two strings. They are globally aligned
    and the edit distance is returned. An indel of any length
    is counted as one edit"""
    aln1, aln2 = _needleman_wunsch(seq1, seq2)
    return edit_distance_from_aln_strings(aln1, aln2)
