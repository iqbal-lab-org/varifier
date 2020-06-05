from cluster_vcf_records import vcf_file_read


def _add_vfr_filter_to_record(record, want_ref_calls=False):
    record.remove_useless_start_nucleotides()
    if len(record.ALT) == 0 or record.ALT == ["."]:
        record.set_format_key_value("VFR_FILTER", "NO_ALTS")
        return
    if "MISMAPPED_UNPLACEABLE" in record.FILTER:
        record.set_format_key_value("VFR_FILTER", "MISMAPPED_UNPLACEABLE")
        return
    if record.FORMAT is None or "GT" not in record.FORMAT:
        record.set_format_key_value("VFR_FILTER", "NO_GT")
        return
    if record.REF in [".", ""]:
        record.set_format_key_value("VFR_FILTER", "NO_REF_SEQ")
        return

    genotype = record.FORMAT["GT"]
    genotypes = genotype.split("/")
    called_alleles = set(genotypes)

    if (
        len(called_alleles) != 1
        or "." in called_alleles
        or (called_alleles == {"0"} and not want_ref_calls)
    ):
        record.set_format_key_value("VFR_FILTER", "CANNOT_USE_GT")
        return

    allele_index = int(called_alleles.pop())
    if allele_index > 0 and record.ALT[allele_index - 1] == "*":
        record.set_format_key_value("VFR_FILTER", "CANNOT_USE_GT")
        return

    if record.FILTER != {"PASS"}:
        record.set_format_key_value("VFR_FILTER", "FAIL_BUT_TEST")
        return

    record.set_format_key_value("VFR_FILTER", "PASS")


def _fix_cluster_filter_tag(records):
    if len(records) > 1:
        pass_filter = [x for x in records if x.FORMAT["VFR_FILTER"] == "PASS"]
        fail_filter = [x for x in records if x.FORMAT["VFR_FILTER"] == "FAIL_BUT_TEST"]
        for x in fail_filter:
            x.set_format_key_value("VFR_FILTER", "FAIL_CONFLICT")
        if len(pass_filter) > 1:
            for x in pass_filter:
                x.set_format_key_value("VFR_FILTER", "FAIL_CONFLICT")


def _annotate_sorted_list_of_records(records, want_ref_calls=False):
    """Annotated sorted list of VCF records. Assumes they all belong to
    the same CHROM. Changes the records in place. No copying"""
    wanted_filters = {"PASS", "FAIL_BUT_TEST"}
    cluster = []
    cluster_end = None

    for record in records:
        _add_vfr_filter_to_record(record, want_ref_calls=want_ref_calls)
        if record.FORMAT["VFR_FILTER"] not in wanted_filters:
            continue

        if len(cluster) == 0:
            cluster = [record]
            cluster_end = record.ref_end_pos()
        elif record.POS > cluster_end:
            _fix_cluster_filter_tag(cluster)
            cluster = [record]
            cluster_end = record.ref_end_pos()
        else:
            cluster.append(record)
            cluster_end = max(cluster_end, record.ref_end_pos())

    _fix_cluster_filter_tag(cluster)


def add_qc_to_vcf(infile, outfile, want_ref_calls=False):
    """Annotated VCF file with QC info needed for calculating precision and recall.
    Adds various tags to each record."""
    header_lines, vcf_records = vcf_file_read.vcf_file_to_dict(
        infile, remove_useless_start_nucleotides=True
    )
    assert header_lines[-1].startswith("#CHROM")

    with open(outfile, "w") as f:
        print(*header_lines[:-1], sep="\n", file=f)
        print(
            '##FORMAT=<ID=VFR_FILTER,Number=1,Type=String,Description="Initial filtering of VCF record. If PASS, then it is evaluated, otherwise is skipped">',
            file=f,
        )
        print(header_lines[-1], file=f)

        for chrom, records in sorted(vcf_records.items()):
            _annotate_sorted_list_of_records(records, want_ref_calls=want_ref_calls)
            print(*records, sep="\n", file=f)
