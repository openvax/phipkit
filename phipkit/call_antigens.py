'''
Identify epitopes in proteins using PhIP-seq hits.

The algorithm processes each protein and sample independently. For each protein
and sample, the number of hits across the protein are computed, and peaks
(local maxima) are identified iteratively. These are the putative epitopes.

To deal with the fact that often proteomes have very many similar proteins,
after identifying antigens, we annotate some antigens as "redundant". For a group
of antigens all supported by the same hit clones (i.e. the same clones align to
multiple antigens), the one with the most hits will be called non-redundant
and the rest are labelled redundant.

The result file from this analysis has the following columns:

    sample_id
        The sample name

    antigen
        Protein identifier

    start_position, end_position
        Position range in the protein that is being called an epitope

    clone_consensus
        Consensus sequence of phip-seq clones that are hits and align to the
        given position range

    antigen_sequence
        Protein sequence for the given position range

    antigen_matches_consensus
        Whether the clone consensus sequence matches the antigen sequence.
        This is experimental and probably best ignored at this point.

    priority_within_antigen
        For proteins that have multiple epitopes called, they ranked according
        to the amount of evidence for them. 1 is the highest priority, 2 is
        lower priority, etc.

    num_clones
        Number of phip-seq clones (hits) that support this epitope

    clones
        List of phip-seq clones (hits) that support this epitope

    redundant
        Whether this antigen (protein) is determined to be redundant with other
        antigens.

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

from . common import say, reconstruct_antigen_sequences

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
    "--max-samples",
    type=int,
    metavar="N",
    help="Process only N samples. For debugging.")

parser.add_argument(
    "--min-clones-per-antigen",
    default=2,
    type=int,
    metavar="N",
    help="Process only antigens with at least N hit clones. Default: %(default)s")

parser.add_argument(
    "--max-evalue",
    default=None,
    type=float,
    metavar="E",
    help="Process only alignments with evalue <= E.")

parser.add_argument(
    "--max-alignments-per-clone",
    default=None,
    type=int,
    metavar="N",
    help="Consider only at most N alignments per clone.")

parser.add_argument(
    "--out",
    required=True,
    metavar="FILE.csv",
    help="Output file.")


def run(argv=sys.argv[1:]):
    args = parser.parse_args(argv)

    say("Reading blast alignments")
    blast_df = pandas.read_csv(args.blast_results)
    say("Read blast alignments:\n", blast_df)

    if args.max_evalue:
        say("Subselecting by evalue", args.max_evalue)
        print("Old shape", blast_df.shape)
        blast_df = blast_df.loc[
            blast_df.evalue <= args.max_evalue
        ].reset_index(drop=True)
        print("New shape", blast_df.shape)

    if args.max_alignments_per_clone:
        say("Subselecting to %d alignments per clone" % (
            args.max_alignments_per_clone))
        print("Old shape", blast_df.shape)
        mask = pandas.Series(index=blast_df.index, dtype=bool)
        mask[:] = False
        for (clone, sub_blast_df) in tqdm(
                blast_df.groupby("clone"), total=blast_df.clone.nunique()):
            sub_blast_df = sub_blast_df.sort_values("evalue").head(
                args.max_alignments_per_clone)
            mask[sub_blast_df.index] = True

        blast_df = blast_df.loc[mask].reset_index(drop=True)
        print("New shape", blast_df.shape)

    hits_df = pandas.read_csv(args.hits, index_col=0)
    say("Read hits:\n", hits_df)

    if args.max_samples:
        say("Subselecting to", args.max_samples, "samples.")
        selected_samples = pandas.Series(hits_df.sample_id.unique()).sample(
            n=args.max_samples)
        say("Selected samples", selected_samples)
        hits_df = hits_df.loc[
            hits_df.sample_id.isin(selected_samples)
        ]
        print(hits_df)

    results_df = call_antigens(
        blast_df=blast_df,
        sample_to_clones=hits_to_dict(hits_df),
        min_clones_per_antigen=args.min_clones_per_antigen)

    say("Writing overlap hits.")
    results_df.to_csv(args.out, index=False)
    say("Wrote: ", args.out)


def hits_to_dict(hits_df):
    """
    Given a hits_df, return a dict of sample id -> list of hits
    """
    sample_to_clones = {}
    for sample, sub_hits_df in hits_df.groupby("sample_id"):
        sample_to_clones[sample] = sub_hits_df[
            ["clone1", "clone2"]
        ].stack().unique()
    return sample_to_clones


def find_consensus(sequences, threshold=0.7):
    """
    Given aligned sequences, return a string where each position i is the
    consensus value at position i if at least threshold fraction of strings
    have the same character that position i, and "." otherwise.
    """

    if len(sequences) == 0:
        return ""

    num_needed = len(sequences) * threshold
    max_len = max(len(s) for s in sequences)
    sequences = [s.ljust(max_len, ".") for s in sequences]
    result = []
    for i in range(max_len):
        chars = collections.Counter(s[i] for s in sequences)
        (char, count) = chars.most_common()[0]
        if char == "." and count >= num_needed: # most strings have terminated
            break
        result.append(char if count >= num_needed else ".")
    return "".join(result).strip(".")


RESULT_COLUMNS = [
    "sample_id",
    "antigen",
    "start_position",
    "end_position",
    "clone_consensus",
    "antigen_sequence",
    "antigen_matches_consensus",
    "priority_within_antigen",
    "num_clones",
    "clones"
]


def analyze_antigen_for_sample(
        sample_id, antigen, sub_blast_df, sequence):
    """
    Call epitopes for one sample and one antigen.

    Parameters
    ----------
    sample_id : string
    antigen : string
    sub_blast_df : pandas.DataFrame
        Blast hits for clones that are hits in this sample and align to the
        given antigen.
        Expected columns: program, version, search_target, clone, len, hit_num,
        id, accession, title, hsp_num, bit_score, score, evalue, identity,
        positive, query_from, query_to, hit_from, hit_to, align_len, gaps,
        qseq, hseq, midline
    sequence : string

    Returns
    -------
    pandas.DataFrame
    """

    # Each clone contributes equally overall, i.e. sum of contributions is 1
    # for all clones. Usually these are just 1, but if a clone aligns to
    # multiple places in the antigen, the contribution of each overlap is
    # less.
    clone_contributions = (1 / sub_blast_df.clone.value_counts()).to_dict()

    # Priority queue inspired implementation to find regions with maximum
    # number of clones overlapping.
    starts_and_ends = []
    for (idx, hit_from, hit_to, clone) in sub_blast_df[
            ["hit_from", "hit_to", "clone"]].itertuples():
        starts_and_ends.append(
            (idx, hit_from - 1, clone_contributions[clone]))
        starts_and_ends.append(
            (idx, hit_to, -clone_contributions[clone]))

    starts_and_ends = pandas.DataFrame(
        starts_and_ends,
        columns=["blast_idx", "position", "value"]
    ).sort_values("position").reset_index(drop=True)
    starts_and_ends["cum_value"] = starts_and_ends.value.cumsum()

    # At each iteration, we find the position ranges with the most overlapping
    # clones. There can be more than one of these since there can be ties.
    # We add an epitope for this region, and then remove the clones overlapping
    # it from further consideration. We repeat this until there are no more
    # clones.
    results = []
    priority = 1
    while len(starts_and_ends) > 0:
        sub_starts_and_ends = starts_and_ends.loc[
            starts_and_ends.cum_value == starts_and_ends.cum_value.max()
        ]
        starts = sub_starts_and_ends.position.values
        ends = starts_and_ends.loc[
            sub_starts_and_ends.index.values + 1
        ].position.values

        assert len(starts) == len(ends)

        used_by_index = (
            sub_starts_and_ends.index.to_series().map(
                lambda idx: starts_and_ends.head(idx + 1).groupby(
                    "blast_idx").value.sum() > 0)
        ).map(lambda s: s[s].index.values).values

        all_used_blast_idx = set()
        for (start, end, used_blast_idx) in zip(starts, ends, used_by_index):
            used_blast_df = sub_blast_df.loc[used_blast_idx]
            used_clones = used_blast_df.clone.unique()
            all_used_blast_idx.update(used_blast_idx)

            clone_subsequences = [
                qseq[start - (hit_from - 1):][:end - start]
                for (qseq, hit_from) in
                used_blast_df[["qseq", "hit_from"]].itertuples(index=False)
            ]
            consensus = find_consensus(clone_subsequences)
            antigen_subsequence = sequence[start : end]
            matches = bool(re.match(consensus, antigen_subsequence))

            results.append((
                sample_id,
                antigen,
                start,
                end,
                consensus,
                antigen_subsequence,
                matches,
                priority,
                len(used_clones),
                list(used_clones),
            ))

        # Drop entries for clones that were just processed.
        starts_and_ends = starts_and_ends.loc[
            ~starts_and_ends.blast_idx.isin(all_used_blast_idx)
        ].reset_index(drop=True)

        # Recompute cumsum.
        starts_and_ends["cum_value"] = starts_and_ends.value.cumsum()
        priority += 1

    result_df = pandas.DataFrame(
        results,
        columns=RESULT_COLUMNS)

    return result_df


def call_antigens(blast_df, sample_to_clones, min_clones_per_antigen=2):
    """
    blast_df : pandas.DataFrame
        Blast hits, phage clones to antigens.
        Expected columns: program, version, search_target, clone, len, hit_num,
        id, accession, title, hsp_num, bit_score, score, evalue, identity,
        positive, query_from, query_to, hit_from, hit_to, align_len, gaps,
        qseq, hseq, midline
    sample_to_clones : dict of string -> list of string
        Hits. Each sample ID mapped to a list of hit clones
    min_clones_per_antigen : int
        Number of clones required to call a hit for an antigen.
    """
    sample_to_clones = pandas.Series(sample_to_clones)
    clones = set()
    for values in sample_to_clones.values:
        clones.update(values)
    clones = sorted(clones)
    say("Calling antigens. Unique clones:", len(clones), clones)

    sub_blast_df = blast_df.loc[blast_df.clone.isin(clones)]
    clone_to_antigens = sub_blast_df.sort_values("evalue").groupby("clone").title.unique()

    say("Collecting antigen sequences")
    antigen_sequences = reconstruct_antigen_sequences(sub_blast_df)

    all_antigens = set()
    for some in clone_to_antigens.values:
        all_antigens.update(some)

    say("Total antigens", len(all_antigens))

    if len(all_antigens) == 0:
        say("No antigens. Consider increasing FDR in the call-hits command.")
        antigens_df = pandas.DataFrame(columns=RESULT_COLUMNS)
    else:
        antigens_dfs = []
        for (sample_id, sample_clones) in tqdm(sample_to_clones.items()):
            missing_clones = [
                c for c in sample_clones if c not in clone_to_antigens
            ]
            if len(missing_clones) > 0:
                say("Skipping %d unaligned hit clone for %s: %s" % (
                    len(missing_clones), sample_id, " ".join(missing_clones)))
                sample_clones = [
                    c for c in sample_clones if c not in set(missing_clones)
                ]

            sample_antigens = collections.Counter()
            for some in clone_to_antigens[sample_clones].values:
                sample_antigens.update(some)

            sample_antigens = set([
                antigen
                for (antigen, count) in sample_antigens.items()
                if count >= min_clones_per_antigen
            ])

            sample_blast_df = blast_df.loc[
                (blast_df.title.isin(sample_antigens)) &
                (blast_df.clone.isin(sample_clones))
            ]
            say("Processing %d antigens (%d clones) for sample %s" % (
                sample_blast_df.title.nunique(),
                sample_blast_df.clone.nunique(),
                sample_id))
            for (antigen, sub_sample_blast_df) in sample_blast_df.groupby("title"):
                sequence = antigen_sequences[antigen]
                results = analyze_antigen_for_sample(
                    sample_id,
                    antigen,
                    sub_sample_blast_df,
                    sequence=sequence)
                antigens_dfs.append(results)

        antigens_df = pandas.concat(antigens_dfs, ignore_index=True)

    # Annotate calls as redundant or not. First we sort the antigens by number
    # of clones supporting them. Working our way down the list, we label an
    # antigen as "redundant" if it does not add at least min_clones_per_antigen
    # additional new clones beyond those already seen.
    say("Identifying redundant antigens")
    antigen_to_clones = collections.defaultdict(set)
    for _, row in tqdm(antigens_df.iterrows(), total=len(antigens_df)):
        antigen_to_clones[row.antigen].update(row.clones)
    antigen_info = pandas.Series(antigen_to_clones).to_frame()
    antigen_info.columns = ['clones']
    antigen_info["num"] = antigen_info.clones.str.len()
    antigen_info = antigen_info.sort_values("num", ascending=False)

    redundant = []
    seen = set()
    for antigen, row in tqdm(antigen_info.iterrows(), total=len(antigen_info)):
        new = row.clones.difference(seen)
        redundant.append(len(new) < min_clones_per_antigen)
        seen.update(new)

    antigen_info["redundant"] = redundant

    say("Labeled %d antigens (%0.2f%%) as redundant" % (
        antigen_info.redundant.sum(), antigen_info.redundant.mean() * 100.0))

    antigens_df["redundant"] = antigens_df.antigen.map(antigen_info.redundant)

    return antigens_df

