import pyfastaq

from cluster_vcf_records import vcf_file_read


def load_mask_bed_file(mask_bed_file):
    """Loads a BED file of ref seq names, and start and end postiions.
    Returns a dictionary of ref seq name -> set of (0-based) coords in the mask."""
    mask = {}
    with pyfastaq.utils.open_file_read(mask_bed_file) as f:
        for line in f:
            chrom, start, end = line.rstrip().split("\t")
            if chrom not in mask:
                mask[chrom] = set()
            for i in range(int(start), int(end)):
                mask[chrom].add(i)
    return mask


def mask_vcf_file(vcf_in, mask_bed_file, vcf_out):
    """Removes all variants in file vcf_in where REF intersects
    an interval in mask_bed_file. Writes new vcf file vcf_out"""
    # This is a quick hacky implementation that is likely not very fast.
    # Put the coords of the mask into a set, and then for each VCF record,
    # check if the ref position(s) are in the set
    mask = load_mask_bed_file(mask_bed_file)

    with pyfastaq.utils.open_file_read(vcf_in) as f_in, open(vcf_out, "w") as f_out:
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


def file_to_dict_of_seqs(infile):
    """Given a file of sequences, returns a dictionary of
    sequence name -> pyfastaq.sequences.Fasta.
    Anything after first whitespace is removed from the names"""
    seqs = {}
    pyfastaq.tasks.file_to_dict(infile, seqs)
    seqs = {k.split()[0]: v for k, v in seqs.items()}
    for seq in seqs.values():
        seq.id = seq.id.split()[0]
    return seqs
