import pyfastaq


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
        return matches, total_positions
