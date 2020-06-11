import copy
from operator import itemgetter

from cluster_vcf_records import vcf_file_read


def _frs_from_vcf_record(record, cov_key="COV"):
    """Gets the FRS from a VCF record, if it exists. This tag is made by
    minos. If FRS tag not there, infers it from key given by cov_key, where
    the value should be alist of coverages for each allele.
    e.g. if COV=1,3 and GT=0/0, then FRS = 1/4"""

    if "FRS" in record.FORMAT:
        if record.FORMAT["FRS"] == ".":
            return "NA"
        else:
            return float(record.FORMAT["FRS"])

    if cov_key not in record.FORMAT:
        return "NA"

    genotypes = set(record.FORMAT["GT"].split("/"))
    if "." in genotypes or len(genotypes) != 1:
        return "NA"

    allele_index = int(genotypes.pop())
    coverages = [int(x) for x in record.FORMAT[cov_key].split(",")]
    total_cov = sum(coverages)
    if total_cov == 0:
        return 0
    else:
        return coverages[allele_index] / total_cov


def per_record_stats_from_vcf_file(infile):
    """Gathers stats for each record in a VCF file.
    Returns a list of dictionaries of stats. One dict per VCF line.
    List is sorted by ref seq name (CHROM), then position (POS)"""
    stats = []
    wanted_keys = [
        "DP",
        "DPF",
        "FRS",
        "GT_CONF",
        "GT_CONF_PERCENTILE",
        "VFR_ED_RA",
        "VFR_ED_TR",
        "VFR_ED_TA",
        "VFR_FILTER",
        "VFR_ALLELE_LEN",
        "VFR_ALLELE_MATCH_COUNT",
        "VFR_ALLELE_MATCH_FRAC",
        "VFR_RESULT",
    ]
    key_types = {
        "DP": int,
        "DPF": float,
        "GT_CONF": float,
        "GT_CONF_PERCENTILE": float,
        "FRS": float,
        "VFR_ED_RA": int,
        "VFR_ED_TR": int,
        "VFR_ED_TA": int,
        "VFR_ALLELE_MATCH_FRAC": float,
        "VFR_ALLELE_LEN": int,
        "VFR_ALLELE_MATCH_COUNT": int,
    }
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(infile)
    for record in vcf_records:
        record_stats = {x: record.FORMAT.get(x, "NA") for x in wanted_keys}
        record_stats["FRS"] = _frs_from_vcf_record(record)
        record_stats["CHROM"] = record.CHROM
        record_stats["POS"] = record.POS + 1
        for key, key_type in key_types.items():
            try:
                record_stats[key] = key_type(record_stats[key])
            except:
                pass

        stats.append(record_stats)

    stats.sort(key=itemgetter("CHROM", "POS"))
    return stats


def format_dict_to_edit_dist_scores(stats):
    if stats["VFR_RESULT"] == "CANNOT_USE_GT":
        return None, None

    # Can have cases where VFR_ED_TR is not present. This is the edit
    # distance between the truth and ref allele. Calculated from mapping
    # of ref probe to the truth.
    # Cases are:
    #  - FP because the alt probe does not map to the truth ref. In this
    #    case we don't map the ref probe to the truth.
    #  - FP, alt probe mapped, but the truth probe did not map
    #  - TP or Partial_TP, and we have VFR_ED_RA,VFR_ED_TA but not VFR_ED_TR.
    #    This is where alt probe mapped, but ref probe did not map, or no
    #    mapping found in same place as the alt mapping
    #    (either way, is counted as no hit).
    if stats["VFR_RESULT"].startswith("FP"):
        if stats["VFR_ED_TR"] in [0, "NA"]:
            return 0, stats["VFR_ED_RA"]
    else:
        if stats["VFR_ED_TA"] == 0:
            return stats["VFR_ED_RA"], stats["VFR_ED_RA"]
        elif stats["VFR_ED_TR"] == 0:
            return 0, stats["VFR_ED_RA"]

    # If we're here and don't know VFR_ED_TR, it's because the ref probe is very
    # different from the truth probe. It's almost certainly a big indel.
    # Only thing left we can do is return the proportion of the allele that
    # matches the truth
    if stats["VFR_ED_TR"] == "NA":
        return stats["VFR_ALLELE_MATCH_FRAC"], 1
    else:
        return stats["VFR_ED_TR"] - stats["VFR_ED_TA"], stats["VFR_ED_TR"]


def summary_stats_from_per_record_stats(per_record_stats, for_recall=False):
    """Given a list of stats made by per_record_stats_from_vcf_file(),
    returns a dictionary of summary stats. Set for_recall to True if the
    VCF was made for getting recall"""
    default_counts = {k: 0 for k in ("Count", "SUM_ALLELE_MATCH_FRAC", "SUM_EDIT_DIST")}
    stats = {"UNUSED": {"CONFLICT": 0, "OTHER": 0}}

    # By default, this is for getting the precision. Which means counting up
    # TPs and FPs. For recall, each call is an expected call from the truth.
    # This means TP is the same as for precision, but a wrong call now means
    # a FN. We expected to find the variant, but didn't.
    fp_key = "FN" if for_recall else "FP"

    for key in "ALL", "FILT":
        stats[key] = {
            "TP": copy.copy(default_counts),
            fp_key: copy.copy(default_counts),
        }
        stats[key]["EDIT_DIST_COUNTS"] = {"numerator": 0, "denominator": 0}

    for d in per_record_stats:
        if d["VFR_FILTER"] == "FAIL_CONFLICT":
            stats["UNUSED"]["CONFLICT"] += 1
        elif d["VFR_FILTER"] not in ["PASS", "FAIL_BUT_TEST"]:
            stats["UNUSED"]["OTHER"] += 1
        else:
            if d["VFR_RESULT"] == "TP":
                result = "TP"
            else:
                result = fp_key
            keys_to_update = ["ALL"]
            if d["VFR_FILTER"] == "PASS":
                keys_to_update.append("FILT")

            ed_num, ed_den = format_dict_to_edit_dist_scores(d)

            for key in keys_to_update:
                try:
                    stats[key][result]["SUM_ALLELE_MATCH_FRAC"] += d[
                        "VFR_ALLELE_MATCH_FRAC"
                    ]
                except TypeError:  # the value could be "NA"
                    pass

                stats[key][result]["SUM_EDIT_DIST"] += d["VFR_ED_RA"]
                stats[key][result]["Count"] += 1

                if ed_num is not None:
                    stats[key]["EDIT_DIST_COUNTS"]["numerator"] += ed_num
                    stats[key]["EDIT_DIST_COUNTS"]["denominator"] += ed_den

    return stats
