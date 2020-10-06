import pyfastaq

from varifier import edit_distance


class Probe:
    def __init__(self, seq, allele_start, allele_end):
        self.seq = seq
        self.allele_start = allele_start
        self.allele_end = allele_end

    def __str__(self):
        return f"{self.seq} {self.allele_start} {self.allele_end}"

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def allele_seq(self):
        return self.seq[self.allele_start : self.allele_end + 1]

    def map_hit_includes_allele(self, map_hit):
        if map_hit.strand == -1:
            start = len(self.seq) - map_hit.q_en
            end = len(self.seq) - map_hit.q_st
        else:
            start = map_hit.q_st
            end = map_hit.q_en

        return start < self.allele_start and self.allele_end < end

    def allele_match_counts(self, map_hit):
        """Given a mappy minimap2 hit, works out how many positions in the
        alignment between the allele and the reference match.
        Returns a tuple: (matching bases, total positions).
        The minimap2 hit must have the extended cigar string, not the 'normal'
        cigar string."""
        if map_hit.NM == 0:
            l = self.allele_end - self.allele_start + 1
            return l, l

        # If there are mismatches, then we can use the extended cigar string to
        # work out how many of those mismatches are in the allele.
        # Example cigar to remind which way round I and D are:
        # read: AGT--TGATCAAGTAC
        #  ref: AGTGATGATC----AC
        # cigar: 3M2D5M4I2M
        probe_pos = map_hit.q_st
        total_positions = 0
        matches = 0

        if map_hit.strand == -1:
            map_hit.cigar.reverse()

        for length, operator in map_hit.cigar:
            if probe_pos > self.allele_end:
                break

            if operator == 7 or operator == 8:  # 7,8 are "=","X"  == match/mismatch
                for i in range(length):
                    if self.allele_start <= probe_pos <= self.allele_end:
                        if operator == 7:
                            matches += 1
                        total_positions += 1
                    probe_pos += 1
                    if probe_pos > self.allele_end:
                        break
            elif operator == 1:  # 1 = I = insertion
                if self.allele_start <= probe_pos <= self.allele_end:
                    total_positions += length
                probe_pos += length
            elif operator == 2:  # 2 = D = deletion
                if self.allele_start <= probe_pos <= self.allele_end:
                    total_positions += length
            else:
                raise RuntimeError(
                    f"Unexpected cigar operator number {operator} with length {length} from cigar"
                )

        if map_hit.strand == -1:
            map_hit.cigar.reverse()

        return matches, total_positions

    def padded_probe_or_ref_seq(self, map_hit, ref_seq=None, ref_mask=None):
        """Returns a tuple: (padded seq string, mask list of bools).
        padded seq string is the padded probe seq inferred from map_hit, or
        if ref_seq provided then the padded ref seq matching the probe.
        If ref_mask is given, should be a set of positions in the mask.
        The returned mask list of bools is same length as the returned padded
        seq string, and has True or False for whether each position is in the mask"""
        # Cigar operators:
        # 1  I  Insertion in query (pad in ref)
        # 2  D  Deletion in query (pad in query)
        # 7  =  Match
        # 8  X  Mismatch
        if map_hit.strand == -1:
            q_st = len(self.seq) - map_hit.q_en
        else:
            q_st = map_hit.q_st

        padded_seq = []
        padded_mask = []

        if ref_seq is None:
            assert ref_mask is None
            ref_seq = self.seq
            pad_operators = {2}
            non_pad_operators = {1, 7, 8}
            if map_hit.strand == -1:
                ref_seq = pyfastaq.sequences.Fasta("probe", self.seq)
                ref_seq.revcomp()
            position = q_st
        else:
            pad_operators = {1}
            non_pad_operators = {2, 7, 8}
            position = map_hit.r_st

        for operator_length, operator_type in map_hit.cigar:
            if operator_type in pad_operators:
                padded_seq.append("-" * operator_length)
                if ref_mask is not None:
                    padded_mask.extend([False] * operator_length)
            elif operator_type in non_pad_operators:
                padded_seq.append(ref_seq[position : position + operator_length])
                if ref_mask is not None:
                    for i in range(position, position + operator_length):
                        padded_mask.append(i in ref_mask)
                position += operator_length
            else:
                raise RuntimeError(
                    f"Unexpected cigar operator number {operator_type} with length {operator_length} from cigar"
                )

        if map_hit.strand == -1:
            padded_seq.extend(["N"] * map_hit.q_st)
            padded_seq = pyfastaq.sequences.Fasta("seq", "".join(padded_seq))
            padded_seq.revcomp()
            padded_seq = padded_seq.seq
            if ref_mask is not None:
                padded_mask.extend([False] * map_hit.q_st)
                padded_mask.reverse()
        else:
            padded_seq = "".join(["N"] * map_hit.q_st + padded_seq)
            if ref_mask is not None:
                padded_mask = [False] * map_hit.q_st + padded_mask

        if ref_mask is None:
            padded_mask = [False] * len(padded_seq)
        assert len(padded_seq) == len(padded_mask)
        return padded_seq, padded_mask

    def padded_seq_allele_start_end_coords(self, padded_seq):
        position = 0
        allele_start = None

        for index, base in enumerate(padded_seq):
            if base != "-":
                if position == self.allele_start:
                    allele_start = index
                if position == self.allele_end:
                    return allele_start, index

                position += 1

        return None, None

    def edit_distance_vs_ref(self, map_hit, ref_seq, ref_mask=None):
        padded_probe_seq, _ = self.padded_probe_or_ref_seq(map_hit)
        padded_ref_seq, padded_ref_mask = self.padded_probe_or_ref_seq(
            map_hit, ref_seq=ref_seq, ref_mask=ref_mask
        )
        start, end = self.padded_seq_allele_start_end_coords(padded_probe_seq)
        if start == None:
            return -1, False
        probe_allele = padded_probe_seq[start : end + 1]
        ref_allele = padded_ref_seq[start : end + 1]
        if ref_mask is None:
            in_mask = False
        else:
            in_mask = any(padded_ref_mask[start : end + 1])
        return (
            edit_distance.edit_distance_from_aln_strings(probe_allele, ref_allele),
            in_mask,
        )
