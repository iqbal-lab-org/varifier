from cluster_vcf_records import vcf_file_read


def read_vcf_file(
    vcf_file,
    homozygous_only=True,
    remove_asterisk_alts=True,
    remove_useless_start_nucleotides=True,
    discard_ref_calls=True,
):
    """Reads input VCF file, yielding vcf records, and a string containing
    either "PASS" if we can use that record, or if not then a reason
    why it can't be used."""
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    yield header_lines

    for record in vcf_records:
        record.remove_useless_start_nucleotides()
        if len(record.ALT) == 0 or record.ALT == ["."]:
            yield record, "NO_ALTS"
            continue
        if "MISMAPPED_UNPLACEABLE" in record.FILTER:
            yield record, "MISMAPPED_UNPLACEABLE"
            continue
        if record.FORMAT is None or "GT" not in record.FORMAT:
            yield record, "NO_GT"
            continue
        if record.REF in [".", ""]:
            yield record, "NO_REF_SEQ"
            continue

        genotype = record.FORMAT["GT"]
        genotypes = genotype.split("/")
        called_alleles = set(genotypes)

        if (
            len(called_alleles) != 1
            or (discard_ref_calls and called_alleles == {"0"})
            or "." in called_alleles
        ):
            yield record, "CANNOT_USE_GT"
            continue

        allele_index = int(called_alleles.pop())
        if allele_index > 0 and record.ALT[allele_index - 1] == "*":
            yield record, "CANNOT_USE_GT"
            continue

        if record.FILTER != {"PASS"}:
            yield record, "FAIL_BUT_TEST"
            continue

        yield record, "PASS"


def mask_vcf_file(vcf_in, mask_bed_file, vcf_out):
    """Removes all variants in file vcf_in where REF intersects
    an interval in mask_bed_file. Writes new vcf file vcf_out"""
    # This is a quick hacky implementation that is likely not very fast.
    # Put the coords of the mask into a set, and then for each VCF record,
    # check if the ref position(s) are in the set
    mask = {}
    with open(mask_bed_file) as f:
        for line in f:
            chrom, start, end = line.rstrip().split("\t")
            if chrom not in mask:
                mask[chrom] = set()
            for i in range(int(start), int(end)):
                mask[chrom].add(i)

    with open(vcf_in) as f_in, open(vcf_out, "w") as f_out:
        for line in f_in:
            if not line.startswith("#"):
                chrom, pos, _, ref, _ = line.split("\t", maxsplit=4)
                in_mask = False
                if chrom in mask:
                    pos = int(pos) - 1
                    for i in range(pos, pos + len(ref)):
                        if i in mask[chrom]:
                            in_mask = True
                            break

                    if in_mask:
                        continue

            print(line, end="", file=f_out)


def vcf_records_are_the_same(file1, file2):
    """Returns True if records in the two VCF files are the same.
    Ignores header lines in the files. Returns False if any lines are different"""
    _, expect_records = vcf_file_read.vcf_file_to_list(file1)
    _, got_records = vcf_file_read.vcf_file_to_list(file2)
    return got_records == expect_records
