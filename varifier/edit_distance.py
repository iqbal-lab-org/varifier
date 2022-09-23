from Bio import pairwise2


def last_gap_end_in_string(string):
    for i, c in enumerate(reversed(string)):
        if c == "-":
            return i
    return -1


def needleman_wunsch(
    seq1,
    seq2,
    match=1,
    mismatch=-1,
    gap_open=-5,
    gap_extend=-3,
    at_genome_start=False
):
    """Returns global alignment strings from NM alignment of the
    two sequences. Dashes for gaps"""
    alignments = pairwise2.align.globalms(
        seq1,
        seq2,
        match,
        mismatch,
        gap_open,
        gap_extend,
    )
    # Alignments is a list of tuples. Each tuple has length 5. Entries:
    # 0: seq1 alignment (ie with dashes for indels)
    # 1: seq2 alignemnt
    # 2: alignment score
    # 4, 5: don't know (not using them)
    if len(alignments) == 1:
        return alignments[0][0], alignments[0][1]

    if at_genome_start:
        best_pos = last_gap_end_in_string(alignments[0][1])
    else:
        best_pos = alignments[0][1].find("-")

    best = alignments[0]

    for a in alignments[1:]:
        if at_genome_start:
            gap_pos = last_gap_end_in_string(a[1])
        else:
            gap_pos = a[1].find("-")

        if gap_pos > best_pos:
            best = a
            best_pos = gap_pos


    return best[0], best[1]


def edit_distance_from_aln_strings(str1, str2):
    """Input should be seqs output by needleman_wunsch().
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
    aln1, aln2 = needleman_wunsch(seq1, seq2)
    return edit_distance_from_aln_strings(aln1, aln2)
