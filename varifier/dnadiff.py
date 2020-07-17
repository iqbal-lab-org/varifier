from operator import attrgetter
import logging
import os
import subprocess

import pyfastaq
import pymummer
from cluster_vcf_records import vcf_record
from varifier import utils

dnadiff_output_extensions = [
    "1coords",
    "1delta",
    "delta",
    "mcoords",
    "mdelta",
    "qdiff",
    "rdiff",
    "report",
    "snps",
    "unqry",
    "unref",
]


def _run_dnadiff(ref_fasta, query_fasta, outprefix):
    command = f"dnadiff -p {outprefix} {ref_fasta} {query_fasta}"
    logging.info(f"Finding variants using dnadiff with command: {command}")
    subprocess.check_output(command, shell=True)
    logging.info(f"dnadiff command finished ({command})")


def _snps_file_to_vcf(snps_file, query_fasta, outfile):
    """Loads the .snps file made by dnadiff.
    query_fasta = fasta file of query sequences.
    Writes a new VCF file unmerged records."""
    vcf_records = {}
    variants = pymummer.snp_file.get_all_variants(snps_file)
    query_seqs = utils.file_to_dict_of_seqs(query_fasta)

    for variant in variants:
        # If the variant is reversed, it means that either the ref or query had to be
        # reverse complemented when aligned by mummer. Need to do the appropriate
        # reverse (complement) fixes so the VCF has the correct REF and ALT sequences
        if variant.reverse:
            qry_seq = pyfastaq.sequences.Fasta("x", variant.qry_base)
            qry_seq.revcomp()
            variant.qry_base = "".join(reversed(qry_seq.seq))
            ref_seq = pyfastaq.sequences.Fasta("x", variant.ref_base)
            ref_seq.revcomp()
            variant.ref_base = ref_seq.seq

        if variant.var_type == pymummer.variant.SNP:
            new_record = vcf_record.VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start + 1),
                        ".",
                        variant.qry_base,
                        variant.ref_base,
                        ".",
                        ".",
                        "SVTYPE=DNADIFF_SNP",
                        "GT",
                        "1/1",
                    ]
                )
            )
        elif variant.var_type == pymummer.variant.DEL:
            # The query has sequence missing, compared to the
            # reference. We're making VCF records w.r.t. the
            # query, so this is an insertion. So need to
            # get the nucleotide before the insertion as well.
            new_record = vcf_record.VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start + 1),
                        ".",
                        query_seqs[variant.qry_name][variant.qry_start],
                        query_seqs[variant.qry_name][variant.qry_start]
                        + variant.ref_base,
                        ".",
                        ".",
                        "SVTYPE=DNADIFF_INS",
                        "GT",
                        "1/1",
                    ]
                )
            )
        elif variant.var_type == pymummer.variant.INS:
            # The ref has sequence missing, compared to the
            # query. We're making VCF records w.r.t. the
            # query, so this is a deletion. So need to
            # get the nucleotide before the deletion as well.
            new_record = vcf_record.VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start),
                        ".",
                        query_seqs[variant.qry_name][variant.qry_start - 1]
                        + variant.qry_base,
                        query_seqs[variant.qry_name][variant.qry_start - 1],
                        ".",
                        ".",
                        "SVTYPE=DNADIFF_DEL",
                        "GT",
                        "1/1",
                    ]
                )
            )
        else:
            raise Exception("Unknown variant type: " + str(variant))

        assert (
            new_record.REF
            == query_seqs[new_record.CHROM][
                new_record.POS : new_record.POS + len(new_record.REF)
            ]
        )

        if new_record.CHROM not in vcf_records:
            vcf_records[new_record.CHROM] = []

        vcf_records[new_record.CHROM].append(new_record)

    for vcf_list in vcf_records.values():
        vcf_list.sort(key=attrgetter("POS"))

    with open(outfile, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        for seq in query_seqs.values():
            print(f"##contig=<ID={seq.id},length={len(seq)}>", file=f)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", file=f)

        for key, vcf_list in sorted(vcf_records.items()):
            for record in vcf_list:
                print(record, file=f)


def make_truth_vcf(ref_fasta, truth_fasta, outfile, debug=False):
    tmp_outprefix = f"{outfile}.tmp"
    _run_dnadiff(truth_fasta, ref_fasta, tmp_outprefix)
    snps_file = f"{tmp_outprefix}.snps"
    _snps_file_to_vcf(snps_file, ref_fasta, outfile)
    if debug:
        return

    for extension in dnadiff_output_extensions:
        # not all files get written, hence try except pass
        try:
            os.unlink(tmp_outprefix + "." + extension)
        except:
            pass
