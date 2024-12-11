import subprocess

import pyfastaq


def mafft_stdout_to_seqs(mafft_stdout, ref_name):
    # The aligned sequences are printed to stdout. We should have the
    # reference genome and the aligned genome in there (do not assume what order
    # they are output). We just want to extract the aligned seqs.
    seqs = mafft_stdout.split(">")
    # seqs looks like: ["", "ref\nACGT...", "to_align\nACGT...", ...]
    assert len(seqs) >= 3
    assert seqs[0] == ""
    ref_seq = None
    aln_seqs = {}
    seqs_to_parse = [x.split("\n", maxsplit=1) for x in seqs[1:]]
    for (name, seq_str) in seqs_to_parse:
        if name == ref_name:
            ref_seq = seq_str.replace("\n", "")
        else:
            assert name not in aln_seqs
            aln_seqs[name] = seq_str.replace("\n", "")

    return ref_seq, aln_seqs


def run_mafft(ref_seq, qry_seq):
    script = "\n".join(
        [
            f'''ref=">ref\n{ref_seq}"''',
            f'''qry=">qry\n{qry_seq}"''',
            """mafft --add <(echo "$qry") <(echo "$ref")""",
        ]
    )
    p = subprocess.run(
        ["bash"],
        input=script,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise Exception(
            f"Error running mafft. Stdout:\n{p.stdout}\n\nStderr:{p.stderr}"
        )

    ref_seq, aln_seqs = mafft_stdout_to_seqs(p.stdout, "ref")
    assert len(aln_seqs) == 1
    return list(ref_seq.upper()), list(aln_seqs["qry"].upper())
