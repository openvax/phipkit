'''
Blast phip-seq library against itself.
'''
from __future__ import (
    print_function,
    division,
    absolute_import,
)
from . common import say

import sys
import time
import tempfile
import argparse
import os
import json
import shutil
import collections
import subprocess
from tqdm import tqdm
tqdm.monitor_interval = 0  # see https://github.com/tqdm/tqdm/issues/481

import pandas

DEFAULT_BLAST_PARAMETERS = collections.OrderedDict([
    ("evalue", "1e-2"),
    ("word_size", "7"),
    ("max_target_seqs", "50"),
    ("gapopen", "11"),
    ("gapextend", "1"),
    ("threshold", "21"),
    ("comp_based_stats", "2"),
    ("num_threads", "1"),
])

DEFAULTS = {
    "seq_col": "seq_aa",
}


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    "annotations",
    nargs="+",
    help="Annotations CSV (or TSV, possibly compressed) file")
parser.add_argument(
    "--seq-col",
    default=DEFAULTS["seq_col"],
    help="Column name in annotations file for the peptide amino acid sequence. "
    "Default: %(default)s")
parser.add_argument(
    "--max-rows",
    metavar="N",
    type=int,
    help="Process only the first N entries of the annotations file. For testing.")
parser.add_argument(
    "--blast-parameter",
    action="append",
    default=[],
    nargs=2,
    metavar="X",
    help="Override default BLAST parameters with specified name and value. "
    "For example, passing '--blast-parameter evalue 1e-3' will result in the "
    "arguments '-evalue 1e-3' being passed to to blastp. Can be specified more "
    "than once. Default parameters: %s" % " ".join(
        "-%s %s" % pair for pair in DEFAULT_BLAST_PARAMETERS.items()))
parser.add_argument(
    "--reference",
    metavar="FILE.fasta",
    help="Instead of an all-against-all BLAST, align phage clone sequences "
    "to the specified reference fasta file.")
parser.add_argument(
    "--out",
    required=True,
    help="Result CSV file")


def run(argv=sys.argv[1:]):
    args = parser.parse_args(argv)

    annot_df = pandas.concat([
        pandas.read_csv(item, index_col=0, nrows=args.max_rows, sep=None)
        for item in args.annotations
    ])
    say("Read annotations file(s) of combined shape", *annot_df.shape)

    extra_blast_parameters = {}
    for (key, value) in args.blast_parameter:
        extra_blast_parameters[key] = value

    all_blast_hits_df = run_blast(
        annot_df=annot_df,
        seq_col=args.seq_col,
        extra_blast_parameters=extra_blast_parameters,
        reference_fasta_path=args.reference)

    if args.out:
        say("Writing results: %d blast alignments" % len(all_blast_hits_df))
        all_blast_hits_df.to_csv(args.out, index=False)
        say("Wrote", os.path.realpath(args.out))


def run_blast(
        annot_df,
        seq_col=DEFAULTS['seq_col'],
        extra_blast_parameters={},
        reference_fasta_path=None):
    temp_dir_handle = tempfile.TemporaryDirectory()
    temp_dir = temp_dir_handle.name
    say('Created temporary directory', temp_dir)

    say("Writing sequence fasta.")
    with open(os.path.join(temp_dir, "sequences.fna"), "w") as fd:
        for idx, row in annot_df.iterrows():
            fd.write(">%s\n" % idx)
            fd.write("%s\n" % row[seq_col])

    if reference_fasta_path:
        shutil.copy2(
            reference_fasta_path, os.path.join(temp_dir, "reference.fna"))
    else:
        shutil.copy2(
            os.path.join(temp_dir, "sequences.fna"),
            os.path.join(temp_dir, "reference.fna"))


    os.environ["BLASTDB"] = temp_dir
    makeblastdb_arguments = [
        "makeblastdb", "-in", "reference.fna", "-dbtype", "prot"
    ]
    say("Generating BLAST reference with command:")
    say(makeblastdb_arguments)
    say("(in an environment with BLASTDB set to %s)" % temp_dir)
    subprocess.check_call(makeblastdb_arguments, cwd=temp_dir)

    blast_arguments = [
        "blastp",
        "-db", "reference.fna",
        "-query", "sequences.fna",
        "-out", "results.json",
        "-outfmt", "15",
        "-use_sw_tback",
    ]
    blast_parameters = dict(DEFAULT_BLAST_PARAMETERS)
    blast_parameters.update(extra_blast_parameters)

    for (key, value) in blast_parameters.items():
        blast_arguments.append("-%s" % key)
        blast_arguments.append(value)

    say("Calling BLAST with arguments:")
    say(blast_arguments)
    start = time.time()
    subprocess.check_call(blast_arguments, cwd=temp_dir)
    say("BLAST completed in %d sec." % (time.time() - start))

    say("Reading BLAST results")
    blast_df = pandas.DataFrame.from_records(
        [x['report'] for x in json.load(
            open(os.path.join(temp_dir, "results.json"), "r"))['BlastOutput2']])
    blast_df['results'] = blast_df.results.map(lambda d: d['search'])
    blast_df['hits'] = blast_df.results.map(lambda d: d['hits'])

    blast_df.index = blast_df.results.map(lambda d: d['query_title'])
    say("Read BLAST results"), *blast_df.shape

    say("Flattening and writing BLAST hits")
    all_blast_hits_df = []

    for clone, row in tqdm(blast_df.iterrows(), total=len(blast_df)):
        for hit_dict in row.hits:
            result_dict = row.to_dict()
            result_dict['clone'] = clone
            result_dict.update(hit_dict)
            (description,) = result_dict.pop('description')
            result_dict['hit_num'] = result_dict['num']
            result_dict.update(description)
            for hsp_dict in hit_dict['hsps']:
                hsp_result_dict = dict(result_dict)
                hsp_result_dict["hsp_num"] = hsp_dict['num']
                hsp_result_dict.update(hsp_dict)
                all_blast_hits_df.append(hsp_result_dict)

    all_blast_hits_df = pandas.DataFrame(all_blast_hits_df)
    del all_blast_hits_df['hsps']
    del all_blast_hits_df['hits']
    del all_blast_hits_df['num']
    del all_blast_hits_df['results']
    del all_blast_hits_df['params']
    del all_blast_hits_df['reference']

    return all_blast_hits_df

