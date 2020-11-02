from operator import attrgetter
import logging
import os
import shutil
import subprocess

import pyfastaq
import pymummer
from cluster_vcf_records import vcf_record
from varifier import utils


# We only want the .snps file from the dnadiff script from MUMmer. From reading
# the docs inspecting that script, we need to run these commands:
#
# nucmer --maxmatch --delta out.delta ref.fasta query.fasta
# delta-filter -1 out.delta > out.1delta
# show-snps -rlTHC out.1delta > out.snps
#
# This is instead of just running show-snps, which runs several other commands
# in addition to making the snps file.
def _run_dnadiff_one_split(ref_fasta, query_fasta, outfile, threads=1, maxmatch=True):
    delta = f"{outfile}.tmp.delta"
    delta_1 = f"{outfile}.tmp.1delta"
    subprocess.check_output(f"rm -f {delta} {delta_1}", shell=True)
    maxmatch_opt = "--maxmatch" if maxmatch else ""
    commands = [
        f"nucmer --threads {threads} {maxmatch_opt} --delta {delta} {ref_fasta} {query_fasta}",
        f"delta-filter -1 {delta} > {delta_1}",
        f"show-snps -rlTHC {delta_1} > {outfile}",
    ]

    for command in commands:
        logging.info("Start run command: " + command)
        subprocess.check_output(command, shell=True)
        logging.info("Finish run command: " + command)

    os.unlink(delta)
    os.unlink(delta_1)


def _run_dnadiff(
    ref_fasta,
    query_fasta,
    outfile,
    split_query=False,
    debug=False,
    threads=1,
    maxmatch=True,
):
    if not split_query:
        _run_dnadiff_one_split(
            ref_fasta, query_fasta, outfile, threads=threads, maxmatch=maxmatch
        )
    else:
        tmp_snp_files = []
        seq_reader = pyfastaq.sequences.file_reader(query_fasta)
        for seq in seq_reader:
            prefix = f"{outfile}.tmp.split.{len(tmp_snp_files)}"
            tmp_fasta = f"{prefix}.fasta"
            with open(tmp_fasta, "w") as f:
                print(seq, file=f)
            snp_file = f"{prefix}.snps"
            _run_dnadiff_one_split(
                ref_fasta, tmp_fasta, snp_file, threads=threads, maxmatch=maxmatch
            )
            os.unlink(tmp_fasta)
            tmp_snp_files.append(snp_file)

        with open(outfile, "wb") as f_out:
            for snp_file in tmp_snp_files:
                with open(snp_file, "rb") as f_in:
                    shutil.copyfileobj(f_in, f_out)
                if not debug:
                    os.unlink(snp_file)


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


def make_truth_vcf(
    ref_fasta,
    truth_fasta,
    outfile,
    debug=False,
    split_ref=False,
    threads=1,
    maxmatch=True,
):
    snps_file = f"{outfile}.tmp.snps"
    _run_dnadiff(
        truth_fasta,
        ref_fasta,
        snps_file,
        split_query=split_ref,
        debug=debug,
        threads=threads,
        maxmatch=maxmatch,
    )
    _snps_file_to_vcf(snps_file, ref_fasta, outfile)
    if not debug:
        os.unlink(snps_file)
