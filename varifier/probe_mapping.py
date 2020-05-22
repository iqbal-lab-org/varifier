import collections
import os

import mappy
import pyfastaq
from cluster_vcf_records import vcf_file_read

from varifier import edit_distance, probe, vcf_qc_annotate


def _get_wanted_format(use_fail_conflict):
    wanted_format = {"PASS", "FAIL_BUT_TEST"}
    if use_fail_conflict:
        wanted_format.add("FAIL_CONFLICT")
    return wanted_format


def get_probes_and_vcf_records(
    vcf_file, ref_seqs, flank_length, use_fail_conflict=False
):
    """Input vcf_file is assumed to have been made by vcf_qc_annotate.add_qc_to_vcf(),
    so that each record has the FORMAT tag VFR_FILTER.
    For each line of the input VCF file, yields a
    tuple (vcf_record, alt probe sequence).
    vcf_file = name of VCF file.
    ref_seqs = dictionary of sequence name -> sequence.
    flank_length = number of nucleotides to add either side of variant sequence."""
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    yield header_lines
    wanted_format = _get_wanted_format(use_fail_conflict)

    for record in vcf_records:
        if record.FORMAT["VFR_FILTER"] not in wanted_format:
            yield record, None, None
            continue

        flank_start = max(0, record.POS - flank_length)
        ref_seq = ref_seqs[record.CHROM]
        if ref_seq[record.POS : record.POS + len(record.REF)] != record.REF:
            record.set_format_key_value("VFR_FILTER", "REF_STRING_MISMATCH")
            yield record, None, None
            continue

        flank_end = min(len(ref_seq) - 1, record.ref_end_pos() + flank_length)
        probe_allele_start = record.POS - flank_start

        alt_index = int(record.FORMAT["GT"].split("/")[0])
        alt_allele = record.REF if alt_index == 0 else record.ALT[alt_index - 1]
        alt_probe_allele_end = probe_allele_start + len(alt_allele) - 1
        alt_probe_seq = (
            ref_seq[flank_start : record.POS]
            + alt_allele
            + ref_seq[record.ref_end_pos() + 1 : flank_end + 1]
        )
        alt_probe = probe.Probe(alt_probe_seq, probe_allele_start, alt_probe_allele_end)
        assert alt_probe.allele_seq() == alt_allele

        ref_probe_allele_end = probe_allele_start + len(record.REF) - 1
        ref_probe_seq = (
            ref_seq[flank_start : record.POS]
            + record.REF
            + ref_seq[record.ref_end_pos() + 1 : flank_end + 1]
        )
        ref_probe = probe.Probe(ref_probe_seq, probe_allele_start, ref_probe_allele_end)
        assert ref_probe.allele_seq() == record.REF

        yield record, ref_probe, alt_probe


def probe_hits_to_best_allele_counts(probe, hits, debug_outfile=None):
    best = None, None, None
    for hit in hits:
        matches, length = probe.allele_match_counts(hit)
        if best[0] is None or matches > best[0]:
            best = matches, length, hit
    return best


def evaluate_vcf_record(
    mapper, vcf_record, ref_probe, alt_probe, map_outfile=None, use_fail_conflict=False
):
    if vcf_record.FORMAT["VFR_FILTER"] not in _get_wanted_format(use_fail_conflict):
        return

    ed = edit_distance.edit_distance_between_seqs(
        ref_probe.allele_seq(), alt_probe.allele_seq()
    )
    vcf_record.set_format_key_value("VFR_EDIT_DIST", str(ed))

    alt_hits = list(mapper.map(alt_probe.seq, MD=True))

    if map_outfile is not None:
        print("VCF", vcf_record, sep="\t", file=map_outfile)
        print("PROBE", alt_probe.seq, sep="\t", file=map_outfile)
        for hit in alt_hits:
            print(
                "PROBE_HIT",
                f"ctg={hit.ctg}",
                f"strand={hit.strand}",
                f"qstart/end={hit.q_st}/{hit.q_en}",
                f"rstart/end={hit.r_st}/{hit.r_en}",
                f"cigar={hit.cigar}",
                sep="\t",
                file=map_outfile,
            )

    alt_match, alt_allele_length, alt_best_hit = probe_hits_to_best_allele_counts(
        alt_probe, alt_hits, debug_outfile=map_outfile
    )

    if alt_match is None:
        vcf_record.set_format_key_value("VFR_RESULT", "FP_PROBE_UNMAPPED")
        return

    vcf_record.set_format_key_value("VFR_ALLELE_LEN", str(alt_allele_length))
    vcf_record.set_format_key_value("VFR_ALLELE_MATCH_COUNT", str(alt_match))
    match_frac = round(alt_match / alt_allele_length, 5) if alt_allele_length > 0 else 0
    vcf_record.set_format_key_value("VFR_ALLELE_MATCH_FRAC", str(match_frac))
    if alt_match == 0 or alt_allele_length == 0:
        result = "FP"
    elif match_frac == 1:
        ref_hits = list(mapper.map(ref_probe.seq, MD=True))
        ref_hits = [
            x
            for x in ref_hits
            if alt_best_hit.ctg == x.ctg
            and x.r_st == alt_best_hit.r_st
            and x.NM < alt_best_hit.NM
        ]
        if len(ref_hits) > 0:
            result = "FP_REF_PROBE_BETTER_MATCH"
        else:
            result = "TP"
    else:
        result = "Partial_TP"
    vcf_record.set_format_key_value("VFR_RESULT", result)


def annotate_vcf_with_probe_mapping(
    vcf_in,
    vcf_ref_fasta,
    truth_ref_fasta,
    flank_length,
    vcf_out,
    map_outfile=None,
    use_fail_conflict=False,
    use_ref_calls=False,
    debug=False,
):
    vcf_ref_seqs = {}
    pyfastaq.tasks.file_to_dict(vcf_ref_fasta, vcf_ref_seqs)
    vcf_ref_seqs = {x.split()[0]: vcf_ref_seqs[x] for x in vcf_ref_seqs}
    vcf_with_qc = vcf_out + ".debug.vcf"
    vcf_qc_annotate.add_qc_to_vcf(vcf_in, vcf_with_qc, want_ref_calls=use_ref_calls)
    probes_and_vcf_reader = get_probes_and_vcf_records(
        vcf_with_qc, vcf_ref_seqs, flank_length, use_fail_conflict=use_fail_conflict,
    )

    # Some notes on the mapper options...
    #
    # From the docs: score is the "scoring system. It is a tuple/list consisting
    # of 4, 6 or 7 positive integers. The first 4 elements specify match scoring,
    # mismatch penalty, gap open and gap extension penalty. The 5th and 6th
    # elements, if present, set long-gap open and long-gap extension penalty.
    # The 7th sets a mismatch penalty involving ambiguous bases."
    # The default mappy Python API do not work. In the tests, results in mappings
    # that make FPs turn into TPs.
    # The options 1,1,5,3 are actually the defaults from bowtie2 and seem to work.
    #
    # k=15 and w=10 are the CLI defaults. On the test data, these values result
    # in the probes near the start of the genome getting mapped, whereas those
    # probes do not get mapped using whatever the Python defaults are.
    #
    # extra_flags=0x4000000 turns on extended cigars, which we use to more easily
    # determine where the matches and mismatches are between the probe and truth
    # reference.
    mapper = mappy.Aligner(
        fn_idx_in=truth_ref_fasta,
        k=15,
        w=10,
        preset="sr",
        n_threads=1,
        extra_flags=0x4000000,
        scoring=[1, 1, 5, 3],
    )
    header_lines = next(probes_and_vcf_reader)

    if map_outfile is not None:
        f_map = open(map_outfile, "w")
    else:
        f_map = None

    new_header_lines = [
        '##FORMAT=<ID=VFR_RESULT,Number=1,Type=String,Description="FP, TP, or Partial_TP when part of the allele matches the truth reference>"',
        '##FORMAT=<ID=VFR_ALLELE_LEN,Number=1,Type=Integer,Description="Number of positions in allele that were checked if they match the truth">',
        '##FORMAT=<ID=VFR_ALLELE_MATCH_COUNT,Number=1,Type=String,Description="Number of positions in allele that match the truth">',
        '##FORMAT=<ID=VFR_ALLELE_MATCH_FRAC,Number=1,Type=String,Description="Fraction of positions in allele that match the truth">',
    ]

    with open(vcf_out, "w") as f_vcf:
        print(
            *header_lines[:-1],
            *new_header_lines,
            header_lines[-1],
            sep="\n",
            file=f_vcf,
        )

        for (vcf_record, ref_probe, alt_probe) in probes_and_vcf_reader:
            evaluate_vcf_record(
                mapper,
                vcf_record,
                ref_probe,
                alt_probe,
                map_outfile=f_map,
                use_fail_conflict=use_fail_conflict,
            )
            print(vcf_record, file=f_vcf)

    if map_outfile is not None:
        f_map.close()

    if not debug:
        os.unlink(vcf_with_qc)
