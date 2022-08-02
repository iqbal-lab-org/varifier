from operator import itemgetter, attrgetter
import os
import sys

from varifier import edit_distance, utils
import pymummer
import pyfastaq

ACGT = {"A", "C", "G", "T"}


def get_perfect_matches(ref_fasta, query_fasta, tmp_nucmer_filename, debug=False):
    # breaklen=1 stops nucmer from joining up matches to give longer
    # <100% matches. Instead, it reports individual shorter 100% hits, which
    # is what we need.
    nucmer_runner = pymummer.nucmer.Runner(
        ref_fasta, query_fasta, tmp_nucmer_filename, breaklen=1, min_id=100
    )
    nucmer_runner.run()
    file_reader = pymummer.coords_file.reader(tmp_nucmer_filename)
    matches = [x for x in file_reader if x.percent_identity == 100]
    if not debug:
        os.unlink(tmp_nucmer_filename)
    return matches


def perfect_matches_to_conservative_match_coords(matches_in, trim=5):
    matches_out = []
    for m in matches_in:
        ref_start = min(m.ref_start, m.ref_end)
        ref_end = max(m.ref_start, m.ref_end)
        qry_start = min(m.qry_start, m.qry_end)
        qry_end = max(m.qry_start, m.qry_end)
        if ref_end - ref_start < 1 + 2 * trim:
            continue

        if ref_start != 0 or qry_start != 0:
            ref_start += trim
            qry_start += trim
        if ref_end != m.ref_length - 1 or qry_end != m.qry_length - 1:
            ref_end -= trim
            qry_end -= trim

        # Can have indels that result in the perfect matches overlapping
        # a little. If this happens, ignore the current match. Is fine because
        # we fill in the gaps later with a proper alignment anyway
        if len(matches_out) > 0 and (
            ref_start <= matches_out[-1]["ref_end"] + 1
            or qry_start <= matches_out[-1]["qry_end"] + 1
        ):
            continue
        matches_out.append(
            {
                "ref_start": ref_start,
                "ref_end": ref_end,
                "qry_start": qry_start,
                "qry_end": qry_end,
            }
        )
    if len(matches_out) == 0:
        raise Exception(
            "Error making global alignment. No matches from nucmer. Cannot continue"
        )

    matches_out.sort(key=itemgetter("ref_start"))
    for i in range(len(matches_out) - 1):
        if matches_out[i]["qry_start"] >= matches_out[i + 1]["qry_start"]:
            print(
                "Cannot align two genomes with rearrangements. Found these two matches:",
                file=sys.stderr,
            )
            print(matches_out[i], file=sys.stderr)
            print(matches_out[i + 1], file=sys.stderr)
            raise NotImplementedError(
                "Error making global alignment because looks like a rearrangement"
            )
    return matches_out


def fix_query_gaps_in_msa(ref_seq, qry_seq):
    """Returns new sequences, with corrected length of gaps in the qry sequence
    alignment qry_seq, so that no gaps are also indels. eg:
    ref: ACGT-ACTGC-GTTTGAGTAGT
    qry: ACGTNNNNGCAG-T--NNNNGT
    becomes:
    ref: ACGTACTGC-GTTTGAGTAGT
    qry: ACGTNNNGCAG-TNNNNNNGT
    """
    # To fix where gap in qry is too short, we need to find dashes next to Ns,
    # and replace each dash with an N. Don't assume dashes always before Ns
    pos = 0
    while pos < len(qry_seq):
        try:
            pos = qry_seq.index("N", pos)
        except ValueError:
            break
        # If we're here, we've found the start of the next gap.
        # Replace any immediately previous dashes with an N
        i = pos - 1
        while i >= 0 and qry_seq[i] == "-":
            qry_seq[i] = "N"
            i -= 1

        # Get to the end of the gap, and replace any following dashes with Ns
        while pos < len(qry_seq) and qry_seq[pos] == "N":
            pos += 1
        while pos < len(qry_seq) and qry_seq[pos] == "-":
            qry_seq[pos] = "N"
            pos += 1

    # To fix where gap is too long, we're looking for "-" in the ref, and
    # "N" in the query. Delete all those positions from both sequences
    bad_indexes = set(
        [i for i in range(len(qry_seq)) if qry_seq[i] == "N" and ref_seq[i] == "-"]
    )
    ref_seq = [x for i, x in enumerate(ref_seq) if i not in bad_indexes]
    qry_seq = [x for i, x in enumerate(qry_seq) if i not in bad_indexes]
    return ref_seq, qry_seq


def homopolymer_char(ref_seq, qry_seq, i):
    if ref_seq[i] == qry_seq[i] and ref_seq[i] in ACGT:
        return ref_seq[i]
    elif ref_seq[i] == "-" and qry_seq[i] in ACGT:
        return qry_seq[i]
    elif qry_seq[i] == "-" and ref_seq[i] in ACGT:
        return ref_seq[i]
    else:
        return None


def is_homopolymer_block(ref_seq, qry_seq, min_length, start, end=None):
    ref_count = len([x for x in ref_seq[start:end] if x in ACGT])
    qry_count = len([x for x in qry_seq[start:end] if x in ACGT])
    return (
        ref_count > 0
        and qry_count > 0
        and (ref_count >= min_length or qry_count >= min_length)
    )


def find_homopolymer_blocks(ref_seq, qry_seq, min_poly_length):
    if min_poly_length < 3:
        raise NotImplementedError(f"min_poly_length must be at least 3")

    blocks = []
    current_char = homopolymer_char(ref_seq, qry_seq, 0)
    start = None if current_char is None else 0
    for i in range(1, len(ref_seq)):
        new_char = homopolymer_char(ref_seq, qry_seq, i)
        if new_char is None or new_char != current_char:
            if start is not None and is_homopolymer_block(
                ref_seq, qry_seq, min_poly_length, start, end=i
            ):
                blocks.append((start, i - 1, current_char))

            start = None if new_char is None else i
            current_char = new_char

    if start is not None and is_homopolymer_block(
        ref_seq, qry_seq, min_poly_length, start
    ):
        blocks.append((start, len(ref_seq) - 1, current_char))

    return blocks


def fix_homopolymer_indels_in_msa(ref_seq, qry_seq, min_poly_length):
    homopolymers = find_homopolymer_blocks(ref_seq, qry_seq, min_poly_length)
    if len(homopolymers) == 0:
        return ref_seq, qry_seq

    new_ref_seq = []
    new_qry_seq = []

    for i, (start, end, base) in enumerate(homopolymers):
        if i == 0:
            to_add_start = 0
        else:
            to_add_start = homopolymers[i - 1][1] + 1

        new_ref_seq.extend(ref_seq[to_add_start:start])
        new_qry_seq.extend(qry_seq[to_add_start:start])
        homopolymer = [x for x in ref_seq[start : end + 1] if x == base]
        new_ref_seq.extend(homopolymer)
        new_qry_seq.extend(homopolymer)

    new_ref_seq.extend(ref_seq[homopolymers[-1][1] + 1 :])
    new_qry_seq.extend(qry_seq[homopolymers[-1][1] + 1 :])
    return new_ref_seq, new_qry_seq


def global_align(
    ref_fasta,
    query_fasta,
    tmp_nucmer_filename,
    debug=False,
    fix_query_gap_lengths=False,
    hp_min_fix_length=None,
):
    try:
        ref_seq = utils.load_one_seq_fasta_file(ref_fasta)
        qry_seq = utils.load_one_seq_fasta_file(query_fasta)
    except:
        raise Exception(
            "Can only use global align on two FASTA files that each contain one sequence"
        )

    perfect_matches = get_perfect_matches(
        ref_fasta, query_fasta, tmp_nucmer_filename, debug=debug
    )
    matches = perfect_matches_to_conservative_match_coords(perfect_matches)

    if matches[0]["ref_start"] > 0:
        assert matches[0]["qry_start"] > 0
        ref_aln, qry_aln = edit_distance.needleman_wunsch(
            ref_seq[: matches[0]["ref_start"]], qry_seq[: matches[0]["qry_start"]]
        )
        aln_ref_seq = list(ref_aln)
        aln_qry_seq = list(qry_aln)
    else:
        aln_ref_seq = []
        aln_qry_seq = []

    for i, match in enumerate(matches):
        aln_ref_seq.extend(list(ref_seq[match["ref_start"] : match["ref_end"] + 1]))
        aln_qry_seq.extend(list(qry_seq[match["qry_start"] : match["qry_end"] + 1]))
        if i < len(matches) - 1:
            ref_end = matches[i + 1]["ref_start"]
            qry_end = matches[i + 1]["qry_start"]
        elif match["qry_end"] == len(qry_seq) - 1:
            assert match["ref_end"] == len(ref_seq) - 1
            break
        else:
            assert match["ref_end"] < len(ref_seq) - 1
            ref_end = len(ref_seq)
            qry_end = len(qry_seq)

        # When we're aligning the end of the genomes, we don't want to do
        # proper global align because can end up with something like this:
        # ref: GCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # qry: GCTTCTTAGGAG--------------------------------------A
        # Using penalise_seq2_end_gap=False means it does an alignment that does
        # not penalise gaps at the end of the query, and then we get this:
        # ref: GCTTCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # qry: GCTTCTTAGGAGA--------------------------------------
        penalise_end = i < len(matches) - 1
        ref_aln, qry_aln = edit_distance.needleman_wunsch(
            ref_seq[match["ref_end"] + 1 : ref_end],
            qry_seq[match["qry_end"] + 1 : qry_end],
            penalise_seq2_end_gap=penalise_end,
        )
        aln_ref_seq.extend(list(ref_aln))
        aln_qry_seq.extend(list(qry_aln))

    if fix_query_gap_lengths:
        aln_ref_seq, aln_qry_seq = fix_query_gaps_in_msa(aln_ref_seq, aln_qry_seq)

    if hp_min_fix_length is not None:
        aln_ref_seq, aln_qry_seq = fix_homopolymer_indels_in_msa(
            aln_ref_seq, aln_qry_seq, hp_min_fix_length
        )

    ref_len_check = len([x for x in aln_ref_seq if x != "-"])
    qry_len_check = len([x for x in aln_qry_seq if x != "-"])
    if (
        len(aln_ref_seq) != len(aln_qry_seq)
        or ref_len_check != len(ref_seq)
        or (
            (hp_min_fix_length is None and not fix_query_gap_lengths)
            and qry_len_check != len(qry_seq)
        )
    ):
        raise Exception(
            f"Error aligning sequences. Got unexpected length(s) after aligning sequences. Ref length={len(ref_seq)}. Query length={len(qry_seq)}. But non-gap lengths are ref={ref_len_check}, qry={qry_len_check}. Cannot continue"
        )

    return "".join(aln_ref_seq), "".join(aln_qry_seq)


def variants_from_global_alignment(ref_aln, qry_aln, ignore_non_acgt=True):
    variants = []
    assert len(ref_aln) == len(qry_aln)
    ref_coord = 0
    current_variant = None
    acgt = {"A", "C", "G", "T", "-"}

    for i, ref_char in enumerate(ref_aln):
        assert not ref_char == "-" == qry_aln[i]
        if ignore_non_acgt and (ref_char not in acgt or qry_aln[i] not in acgt):
            if current_variant is not None:
                variants.append(current_variant)
                current_variant = None
        elif ref_char != qry_aln[i]:
            if current_variant is None:
                if i > 0 and (ref_char == "-" or qry_aln[i] == "-"):
                    current_variant = {
                        "ref_start": ref_coord - 1,
                        "ref_allele": [ref_aln[i - 1]],
                        "qry_allele": [qry_aln[i - 1]],
                    }
                else:
                    current_variant = {
                        "ref_start": ref_coord,
                        "ref_allele": [],
                        "qry_allele": [],
                    }
            if ref_char != "-":
                current_variant["ref_allele"].append(ref_char)
            if qry_aln[i] != "-":
                current_variant["qry_allele"].append(qry_aln[i])
        elif current_variant is not None:
            variants.append(current_variant)
            current_variant = None

        if ref_char != "-":
            ref_coord += 1

    if current_variant is not None:
        variants.append(current_variant)
        current_variant = None

    # If we had an insertion or deletion at the start, then need to fix because
    # ref or alt allele will be empty
    if len(variants) > 0 and variants[0]["ref_start"] == 0:
        if len(variants[0]["ref_allele"]) == 0:
            new_base = qry_aln[len(variants[0]["qry_allele"])]
            variants[0]["ref_allele"] = [new_base]
            variants[0]["qry_allele"].append(new_base)
        elif len(variants[0]["qry_allele"]) == 0:
            new_base = ref_aln[len(variants[0]["ref_allele"])]
            variants[0]["ref_allele"].append(new_base)
            variants[0]["qry_allele"] = [new_base]

    for d in variants:
        d["ref_allele"] = "".join(d["ref_allele"])
        d["qry_allele"] = "".join(d["qry_allele"])

    return variants


def expand_combined_snps(variants_in):
    variants_out = []
    for variant in variants_in:
        if (
            len(variant["ref_allele"]) == len(variant["qry_allele"])
            and "N" not in variant["ref_allele"]
            and "N" not in variant["qry_allele"]
        ):
            for i in range(len(variant["ref_allele"])):
                variants_out.append(
                    {
                        "ref_start": i + variant["ref_start"],
                        "ref_allele": variant["ref_allele"][i],
                        "qry_allele": variant["qry_allele"][i],
                    }
                )
        else:
            variants_out.append(variant)
    return variants_out


def vcf_using_global_alignment(
    ref_fasta,
    query_fasta,
    vcf_out,
    debug=False,
    fix_query_gap_lengths=False,
    hp_min_fix_length=None,
    fixed_query_fasta=None,
    msa_file=None,
    min_ref_coord=0,
    max_ref_coord=float("inf"),
    ignore_non_acgt=True,
):
    ref_aln, qry_aln = global_align(
        ref_fasta,
        query_fasta,
        f"{vcf_out}.tmp.nucmer",
        debug=debug,
        fix_query_gap_lengths=fix_query_gap_lengths,
        hp_min_fix_length=hp_min_fix_length,
    )

    if msa_file is not None:
        with open(msa_file, "w") as f:
            print(ref_aln, qry_aln, sep="\n", file=f)

    if fixed_query_fasta is not None:
        no_dash_query = "".join([x for x in qry_aln if x != "-"])
        seq = pyfastaq.sequences.Fasta("gap_fixed_query", no_dash_query)
        with open(fixed_query_fasta, "w") as f:
            print(seq, file=f)

    variants = variants_from_global_alignment(
        ref_aln, qry_aln, ignore_non_acgt=ignore_non_acgt
    )
    variants = expand_combined_snps(variants)
    ref_seq = utils.load_one_seq_fasta_file(ref_fasta)
    ref_name = ref_seq.id.split()[0]
    with open(vcf_out, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        print('##FILTER=<ID=PASS,Description="All filters passed">', file=f)
        print(f"##contig=<ID={ref_name},length={len(ref_seq)}>", file=f)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
        print(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "sample",
            sep="\t",
            file=f,
        )
        for i, variant in enumerate(variants):
            end_coord = variant["ref_start"] + len(variant["ref_allele"]) - 1
            if min_ref_coord > variant["ref_start"] or end_coord > max_ref_coord:
                continue
            print(
                ref_name,
                variant["ref_start"] + 1,
                i,  # ID column
                "".join(variant["ref_allele"]),
                "".join(variant["qry_allele"]),
                ".",  # QUAL column
                "PASS",
                ".",  # FILTER column
                "GT",
                "1/1",
                sep="\t",
                file=f,
            )
