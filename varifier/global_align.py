from operator import itemgetter, attrgetter
import os
import sys

from varifier import edit_distance, utils
import pymummer


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
        if len(matches_out) > 0 and (ref_start <= matches_out[-1]["ref_end"] or  qry_start <= matches_out[-1]["qry_end"]):
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


def global_align(ref_fasta, query_fasta, tmp_nucmer_filename, debug=False):
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
        aln_ref_seq = [ref_aln]
        aln_qry_seq = [qry_aln]
    else:
        aln_ref_seq = []
        aln_qry_seq = []

    for i, match in enumerate(matches):
        aln_ref_seq.append(ref_seq[match["ref_start"] : match["ref_end"] + 1])
        aln_qry_seq.append(qry_seq[match["qry_start"] : match["qry_end"] + 1])
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

        ref_aln, qry_aln = edit_distance.needleman_wunsch(
            ref_seq[match["ref_end"] + 1 : ref_end],
            qry_seq[match["qry_end"] + 1 : qry_end],
        )
        aln_ref_seq.append(ref_aln)
        aln_qry_seq.append(qry_aln)

    aln_ref_seq = "".join(aln_ref_seq)
    aln_qry_seq = "".join(aln_qry_seq)
    ref_len_check = len([x for x in aln_ref_seq if x != "-"])
    qry_len_check = len([x for x in aln_qry_seq if x != "-"])
    if (
        len(aln_ref_seq) != len(aln_qry_seq)
        or ref_len_check != len(ref_seq)
        or qry_len_check != len(qry_seq)
    ):
        raise Exception(
            f"Error aligning sequences. Got unexpected length(s) after aligning sequences. Ref length={len(ref_seq)}. Query length={len(qry_seq)}. But non-gap lengths are ref={ref_len_check}, qry={qry_len_check}. Cannot continue"
        )

    return aln_ref_seq, aln_qry_seq


def variants_from_global_alignment(ref_aln, qry_aln):
    variants = []
    assert len(ref_aln) == len(qry_aln)
    ref_coord = 0
    current_variant = None

    for i, ref_char in enumerate(ref_aln):
        assert not ref_char == "-" == qry_aln[i]
        if "N" in [ref_char, qry_aln[i]]:
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
        if len(variant["ref_allele"]) == len(variant["qry_allele"]):
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
    min_ref_coord=0,
    max_ref_coord=float("inf"),
):
    ref_aln, qry_aln = global_align(
        ref_fasta, query_fasta, f"{vcf_out}.tmp.nucmer", debug=debug
    )
    variants = variants_from_global_alignment(ref_aln, qry_aln)
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
