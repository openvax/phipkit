'''
Plot antigen calls
'''
from __future__ import (
    print_function,
    division,
    absolute_import,
)
import sys
import re
import argparse
import collections
import numpy
import statsmodels.stats.multitest
from tqdm import tqdm
tqdm.monitor_interval = 0  # see https://github.com/tqdm/tqdm/issues/481

import pandas

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

from .common import say, hits_to_dict
from .antigen_analysis import AntigenAnalysis

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "blast_results",
    metavar="BLAST_RESULTS.csv",
    help="Blast results against a proteome reference.")

parser.add_argument(
    "hits",
    metavar="HITS.csv",
    help="Hits file generated by call-hits command.")

parser.add_argument(
    "antigen_calls",
    metavar="ANTIGENS.csv",
    help="Calls from call-antigens command.")

parser.add_argument(
    "--max-antigens",
    type=int,
    metavar="N",
    help="Process only N antigens. For debugging.")

parser.add_argument(
    "--min-samples",
    type=int,
    metavar="N",
    help="Process only antigens with hits in at least N samples")

parser.add_argument(
    "--include-redundant",
    action="store_true",
    default=False,
    help="Include redundant antigens")

parser.add_argument(
    "--out",
    required=True,
    metavar="FILE.pdf",
    help="Output PDF.")


def run(argv=sys.argv[1:]):
    args = parser.parse_args(argv)

    say("Reading blast alignments")
    blast_df = pandas.read_csv(args.blast_results)
    say("Read blast alignments:\n", blast_df)

    say("Reading hits")
    hits_df = pandas.read_csv(args.hits, index_col=0)
    say("Read clone hits:\n", hits_df)

    say("Reading antigen calls")
    antigens_df = pandas.read_csv(args.antigen_calls)
    say("Read antigen calls:\n", antigens_df)

    if args.min_samples:
        say("Subselecting to antigens with >=", args.min_samples, "samples")
        samples_per_antigen = antigens_df.groupby("antigen").sample_id.nunique()
        selected_antigens = samples_per_antigen[
            samples_per_antigen  >= args.min_samples
        ].index
        say("Selected antigens", selected_antigens)
        antigens_df = antigens_df.loc[
            antigens_df.antigen.isin(selected_antigens)
        ]
        print(antigens_df)

    if args.max_antigens:
        say("Subselecting to", args.max_antigens, "antigens.")
        selected_antigens = pandas.Series(
            antigens_df.antigen.unique()).sample(n=args.max_antigens)
        say("Selected antigens", selected_antigens)
        antigens_df = antigens_df.loc[
            antigens_df.antigen.isin(selected_antigens.values)
        ]
        print(antigens_df)

    results_df = plot_antigens(
        blast_df=blast_df,
        hits_df=hits_df,
        antigens_df=antigens_df,
        out=args.out,
        include_redundant=args.include_redundant)

    say("Done. Summary")
    print(results_df)


def plot_antigens(blast_df, hits_df, antigens_df, out, include_redundant=False):
    if not include_redundant:
        antigens_df = antigens_df.loc[~antigens_df.redundant]

    analyzer = AntigenAnalysis(
        blast_df=blast_df,
        antigens_df=antigens_df,
        sample_to_hit_clones=hits_to_dict(hits_df))


    say("Generating plots")
    antigens = antigens_df.antigen.unique()
    summary_df = []
    with PdfPages(out) as pdf:
        for (i, antigen) in enumerate(tqdm(antigens)):
            fig = analyzer.plot_antigen(antigen)
            pdf.savefig(fig)
            pyplot.close()
            summary_df.append((antigen, i + 1))
    say("Wrote: ", out)

    summary_df = pandas.DataFrame(summary_df, columns=["antigen", "page"])
    return summary_df


