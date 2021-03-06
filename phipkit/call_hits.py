'''
Analyze phip-seq scores to determine hits.

The idea is that two clones with overlapping sequences that both have very
high scores likely indicate a hit (i.e. an antibody response to the
shared part of their sequences).

We first convert scores to percentile ranks. By percentile rank we mean a
number between 0 and 1 where higher scores have lower ranks. The highest score
in the scores matrix has percentile rank 0. We are looking for clones that have
low enough percentile ranks to be significant. If we draw a clone at random,
its percent rank R is distributed as Uniform(0,1).

For each sample and *pair* of phage clones that have overlapping sequences
(as determined by BLAST), we consider the percentile rank of the clones.
Call these ranks R1 and R2. Intuitively, if these are both very low then we
want to call a hit. Our null hypothesis informally is that R1 and R2 are
independent, since our key assumption is that clones that share sequence will
be correlated only when there is a real antibody response to their shared
sequence.

More formally, we define the null hypothesis to be that the joint distribution
of (R1, R2) is the 2-dimensional unit uniform distribution. We define
R = max(R1, R2) as a test statistic. Under the null hypothesis, for any constant
k, Prob[R < k] = k^2, since this is just the probability of sampling two values
both less than k from a unit uniform distribution.

We can therefore take p=R^2 to be the p-value for this pair of clones.
It gives the probability of observing two clones with percentile ranks as
low as these under the null hypothesis.

These p-values are reported as q values after Benjamini-Hochberg FDR correction.
'''
from __future__ import (
    print_function,
    division,
    absolute_import,
)
import sys
import argparse
import numpy
import statsmodels.stats.multitest
from tqdm import tqdm
tqdm.monitor_interval = 0  # see https://github.com/tqdm/tqdm/issues/481

import pandas

from .common import say

DEFAULTS = {
    'fdr': 0.01,
}

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "blast_results",
    metavar="FILE.csv",
    help="Blast results aligning clones against clones")

parser.add_argument(
    "scores",
    metavar="FILE.csv",
    help="Scores matrix. Rows are clones, columns are samples. First col "
    "gives clone ids.")

parser.add_argument(
    "--counts",
    metavar="FILE.csv",
    help="Clones x samples matrix with read counts. This is not used in "
    "the computation but will be included in the output for easier "
    "interpretation.")

parser.add_argument(
    "--fdr",
    type=float,
    metavar="X",
    default=DEFAULTS['fdr'],
    help="False discovery rate. Overlap hits with FDR < X are considered high "
    "confidence. Default: %(default)s)")

parser.add_argument(
    "--out",
    required=True,
    metavar="FILE.csv",
    help="Output file for hits")


def run(argv=sys.argv[1:]):
    args = parser.parse_args(argv)

    say("Reading blast alignments")
    blast_df = pandas.read_csv(
        args.blast_results, usecols=["clone", "title", "midline"])
    say("Read blast alignments:\n", blast_df)

    scores_df = pandas.read_csv(args.scores, index_col=0)
    say("Read scores:\n", scores_df)

    counts_df = None
    if args.counts:
        counts_df = pandas.read_csv(args.counts, index_col=0, sep=None)

    overlap_hits_df = call_hits(
        blast_df=blast_df,
        scores_df=scores_df,
        counts_df=counts_df,
        fdr=args.fdr)

    say("Writing overlap hits.")
    overlap_hits_df.to_csv(args.out, index=False)
    say("Wrote: ", args.out)


def call_hits(
        blast_df,
        scores_df,
        counts_df=None,
        fdr=DEFAULTS['fdr']):
    blast_df = blast_df.loc[
        blast_df.clone.isin(scores_df.index) &
        blast_df.title.isin(scores_df.index) &
        (blast_df.clone != blast_df.title)
    ]
    say("Subselected to non-self blast hits present in scores. New shape",
        *blast_df.shape)
    if len(blast_df) == 0:
        raise ValueError("No blast alignments")

    clone_to_others = blast_df.groupby("clone").title.unique().map(set)
    midlines = blast_df.groupby(["clone", "title"]).midline.unique()

    del blast_df
    say("Clone map:\n", clone_to_others)

    say("Calculating 2-cliques")
    cliques_set = set()
    for clone, others in tqdm(clone_to_others.items(), total=len(clone_to_others)):
        for other in others:
            if clone in clone_to_others.get(other, {}):
                tpl = tuple(sorted([clone, other]))
                cliques_set.add(tpl)

    cliques = pandas.Series(sorted(cliques_set))
    del cliques_set
    say("Done. Cliques:", len(cliques))

    say("Calculating ranks")
    ranks_df = scores_df.rank(pct=True, ascending=False, method="max")
    say("Done.")

    say("Calculating max_ranks_by_clique")
    first_clone_ranks = ranks_df.loc[cliques.str.get(0)]
    first_clone_ranks.index = cliques

    second_clone_ranks = ranks_df.loc[cliques.str.get(1)]
    second_clone_ranks.index = cliques

    max_ranks_by_clique = first_clone_ranks.where(
        first_clone_ranks > second_clone_ranks, second_clone_ranks)

    del first_clone_ranks
    del second_clone_ranks
    del cliques

    say("Done:\n", max_ranks_by_clique)
    say("Calculating overlap hits")

    p_values = max_ranks_by_clique**2

    hits_df = []

    for sample_id in tqdm(scores_df.columns):
        # Clone sets
        q_values = statsmodels.stats.multitest.multipletests(
            p_values[sample_id],
            method="fdr_bh")[1]
        take = q_values < fdr
        df = pandas.DataFrame({
            "q": q_values[take]
        }, index=p_values.index[take])

        df["p"] = p_values.loc[df.index, sample_id]
        df["max_rank"] = max_ranks_by_clique.loc[df.index, sample_id]
        df["sample_id"] = sample_id
        df["clone1"] = df.index.str.get(0)
        df["clone2"] = df.index.str.get(1)
        hits_df.append(df.reset_index(drop=True))

    hits_df = pandas.concat(hits_df, ignore_index=True)
    hits_df["midline"] = hits_df.set_index(
        ["clone1", "clone2"]).index.map(midlines).map(", ".join)
    say("Done. Overlap hits:\n", hits_df)

    # Annotate scores
    for sample_id, sub_df in hits_df.groupby("sample_id"):
        for clone in ["clone1", "clone2"]:
            hits_df.loc[
                sub_df.index,
                "%s_score" % clone
            ] = sub_df[clone].map(scores_df[sample_id])

    if counts_df is not None:
        say("Annotating counts")
        cpm_df = (
            counts_df.T * 1e6 / numpy.expand_dims(counts_df.sum(), 1)).T
        for sample_id, sub_df in hits_df.groupby("sample_id"):
            for clone in ["clone1", "clone2"]:
                hits_df.loc[
                    sub_df.index,
                    "%s_count" % clone
                ] = sub_df[clone].map(counts_df[sample_id])
                hits_df.loc[
                    sub_df.index,
                    "%s_cpm" % clone
                ] = sub_df[clone].map(cpm_df[sample_id])
                hits_df.loc[
                    sub_df.index,
                    "%s_rank" % clone
                ] = sub_df[clone].map(ranks_df[sample_id])
        say("Done annotating")

    return hits_df.sort_values("p")
