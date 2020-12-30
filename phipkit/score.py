'''
Analyze PHiP-seq read counts matrix to generate enrichment-over-beads-only scores.

Algorithm sketch:

[I] For each bead-only sample:
    [1] Bin the read counts across clones into some number of bins (default 50).
    [2] For each set of clones c associated with each bin:
            For each other sample s:
                [i] Fit a robust estimate of a Poisson distribution parameter (lambda)
                  to the read counts for clones c in sample s. We use both the median
                  and an expression involving the interquartile range, and take
                  whichever estimate for lambda is larger.
                [ii] Calculate the negative log of the survival function (1-CDF)
                  for each observed read count against the Poisson distribution
                  fit in the previous step. High values indicate outliers that
                  are higher than expected.
[II] Take the mean of the values computed above to obtain a score for each
  clone in each sample.

'''
from __future__ import (
    print_function,
    division,
    absolute_import,
)
from . common import say

import sys
import re
import argparse
import numpy
import scipy.stats
from tqdm import tqdm
tqdm.monitor_interval = 0  # see https://github.com/tqdm/tqdm/issues/481

import pandas


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    "counts",
    action="append",
    default=[],
    metavar="FILE.csv",
    help="Read counts matrix. Rows are clones, columns are samples. First col "
    "gives clone ids")

parser.add_argument(
    "--out",
    metavar="FILE.csv",
    required=True,
    help="Output file for enrichment scores")

parser.add_argument(
    "--beads-regex",
    default=".*beads.*",
    help="Samples with names matching this regex are considered beads-only "
    "(mock). Default: %(default)s.")

parser.add_argument(
    "--bins",
    type=int,
    default=50,
    help="Number of bins to use for beads only read counts."
    "Default: %(default)s.")

parser.add_argument(
    "--mask-quantile-difference-threshold",
    type=int,
    default=0.5,
    help="Advanced use. Default: %(default)s.")

parser.add_argument(
    "--clone-min-reads",
    type=int,
    metavar="N",
    default=0,
    help="Consider only clones with at least N reads in at least one sample")

parser.add_argument(
    "--sample-min-total-reads",
    type=int,
    metavar="N",
    help="Exclude samples with fewer than N total reads")

parser.add_argument(
    "--max-rows",
    metavar="N",
    type=int,
    help="Process only N clones. For testing.")


def run(argv=sys.argv[1:]):
    args = parser.parse_args(argv)
    if not args.counts:
        parser.error("A counts matrix is required")

    counts_dfs = []
    for some in args.counts:
        some_df = pandas.read_csv(
            some, sep=None, index_col=0, nrows=args.max_rows)
        print("Read counts ", some_df.shape)
        counts_dfs.append(some_df)

    counts_df = pandas.concat(counts_dfs, axis=1)

    if args.max_rows:
        counts_df = counts_df.sample(n=args.max_rows)

    print("Counts")
    print(counts_df)

    if args.clone_min_reads:
        print("Subselecting by clone-min-reads: %d" % args.clone_min_reads)
        counts_df = counts_df.loc[
            counts_df.max(1) >= args.clone_min_reads
        ].copy()
        print("Counts after filtering:")
        print(counts_df)

    if args.sample_min_total_reads:
        print("Subselecting by sample-min-reads: %d" %
              args.sample_min_total_reads)
        excluded_samples = counts_df.sum()
        excluded_samples = excluded_samples[
            excluded_samples < args.sample_min_total_reads
        ].index.values
        print("Excluding samples: ", excluded_samples)
        counts_df = counts_df[
            [c for c in counts_df.columns if c not in excluded_samples]
        ]

    scores_df = compute_scores(
        counts_df=counts_df,
        bins=args.bins,
        mask_quantile_difference_threshold=args.mask_quantile_difference_threshold,
        beads_regex=args.beads_regex)

    say("Done calculating scores. Writing results.")

    scores_df.to_csv(args.out, index=True, float_format="%0.3f")
    say("Wrote [%s]" % str(scores_df.shape), args.out)


def compute_scores(
        counts_df,
        bins=parser.get_default("bins"),
        mask_quantile_difference_threshold=parser.get_default("mask_quantile_difference_threshold"),
        beads_regex=parser.get_default("beads_regex")):
    beads_only_samples = [
        s for s in counts_df.columns if
        re.match(beads_regex, s, flags=re.IGNORECASE)
    ]
    print(
        "Beads-only regex '%s' matched %d samples: %s" % (
            beads_regex,
            len(beads_only_samples),
            " ".join(beads_only_samples)))

    if len(beads_only_samples) <= 2:
        raise ValueError("At least 2 beads only samples are needed")

    non_beads_only_samples = [
        c for c in counts_df.columns if c not in beads_only_samples
    ]

    say("Fitting empirical distributions")
    bead_sample_to_scores = {}
    for b in tqdm(beads_only_samples):
        beads_discretized = pandas.qcut(
            counts_df[b], bins, duplicates='drop')
        if beads_discretized.isnull().all():
            # Constant series: just set all values to the -1 placeholder.
            beads_discretized = beads_discretized.cat.codes
        grouped = counts_df.groupby(beads_discretized)
        survivals = pandas.DataFrame(index=counts_df.index)
        for (bin, sub_df) in grouped:
            for col in sub_df.columns:
                iqr = scipy.stats.iqr(sub_df[col].values)
                median_lambda_estimate = sub_df[col].median()

                # I determined this expression empirically by plotting the
                # interquartile range (IQR) for Poisson samples against lambda
                # across a range of lambda values.
                iqr_lambda_estimate = (iqr / 1.349)**2
                lambda_estimate = max(
                    1.0, median_lambda_estimate, iqr_lambda_estimate)
                assert not numpy.isnan(lambda_estimate)

                sub_survivals = -scipy.stats.poisson.logsf(
                    sub_df[col].values, mu=lambda_estimate)
                survivals.loc[sub_df.index, col] = numpy.nan_to_num(
                    sub_survivals, posinf=100.0)

                if bin == beads_discretized.cat.categories[-1]:
                    # The last bin is treated slightly differently: we remove
                    # clones that are ranked very highly in terms of survival
                    # but ranked low in terms of the ratio of read counts (
                    # sample to beads only). This is needed because the last
                    # bin doesn't really control for beads only read count,
                    # since the top values in the bin can be arbitrarily high.
                    sub_survival_quantiles = pandas.Series(
                        sub_survivals,
                        index=sub_df.index).rank(pct=True)
                    sub_ratio_quantiles = (sub_df[col] / sub_df[b]).rank(
                        pct=True)
                    quantile_differences = (
                        sub_survival_quantiles - sub_ratio_quantiles)

                    mask_out = quantile_differences[
                        quantile_differences > mask_quantile_difference_threshold
                    ].index.values
                    if len(mask_out) > 0:
                        print("MASKING OUT", b, col, mask_out)
                        survivals.loc[mask_out, col] = 0.0

        survivals[b] = 0.0
        bead_sample_to_scores[b] = survivals

    say("Calculating scores")
    scores_df = sum(bead_sample_to_scores.values())
    scores_df[beads_only_samples] /= len(beads_only_samples) - 1
    scores_df[non_beads_only_samples] /= len(beads_only_samples)
    return scores_df

