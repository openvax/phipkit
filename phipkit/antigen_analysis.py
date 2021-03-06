import collections

import numpy
import pandas
import scipy

from tqdm import tqdm

from .common import say, reconstruct_antigen_sequences


def compute_coverage(antigen, sequence, blast_df):
    """
    Extract blast hits for some clones for a single antigen into a DataFrame
    indicating whether each clone aligns at each position in the antigen.

    Parameters
    ----------
    antigen : str
    sequence : str
        Antigen amino acid sequence
    blast_df : pandas.DataFrame
        Blast results for clones to consider (e.g. may be all clones or just
        those with hits).

    Returns
    -------
    pandas.DataFrame of int

    Columns are positions in the antigen (protein).
    Index are clone names

    Values indicate whether the clone aligns at that position (> 0) or not (0).
    If a clone aligns to the same antigen in multiple ways, the number of ways
    is indicated. E.g. a value of 2 at a position means the clone aligned to that
    position in two different ways.
    """
    blast_df = blast_df.loc[
        blast_df.title == antigen
    ]
    assert len(blast_df) > 0, "No hits for antigen %s" % antigen

    clones = []
    result_rows = []

    for clone, hit_rows in blast_df.groupby("clone"):
        # We take the union over all hsps for a given clone, i.e. not double
        # counting.
        result_row = numpy.zeros(len(sequence), dtype=int)
        for _, row in hit_rows.iterrows():
            pos = row.hit_from - 1
            for i in range(len(row.hseq)):
                match = row.hseq[i] not in ["-", " ", "+"] and (row.hseq[i] == row.qseq[i])
                assert match == bool(row.midline[i] not in ("+", " "))
                if match:
                    assert row.qseq[i] == sequence[pos]
                    assert row.hseq[i] == sequence[pos]
                    result_row[pos] = 1
                if row.hseq[i] != "-":
                    pos += 1
        result_rows.append(result_row)
        clones.append(clone)

    result = pandas.DataFrame(result_rows, index=clones)
    return result


class AntigenAnalysis(object):
    def __init__(self, blast_df, antigens_df, sample_to_hit_clones, sample_to_kind=None):
        self.blast_df = blast_df
        self.antigens_df = antigens_df
        self.sample_to_kind = sample_to_kind

        self.sample_to_hit_clones = dict(
            (k, list(v)) for (k, v) in sample_to_hit_clones.items())

        all_clones = set()
        for clones in self.sample_to_hit_clones.values():
            all_clones.update(clones)

        self.clone_by_sample_hits_matrix = pandas.DataFrame(
            index=sorted(all_clones),
            columns=list(self.sample_to_hit_clones),
            dtype=int)
        self.clone_by_sample_hits_matrix[:] = 0
        for (sample, clones) in self.sample_to_hit_clones.items():
            self.clone_by_sample_hits_matrix.loc[clones, sample] = 1
        say("Done.")

    def plot_antigen(
            self,
            antigen,
            heatmap=True,
            special_features={},
            max_len_to_show_sequence=100):
        import dna_features_viewer
        from matplotlib import pyplot
        import seaborn

        sequence = reconstruct_antigen_sequences(
            self.blast_df.loc[self.blast_df.title == antigen])[antigen]

        all_coverage_by_clone = compute_coverage(
            antigen,
            sequence,
            self.blast_df)

        all_coverage_total = (all_coverage_by_clone > 0).sum()

        all_samples = list(self.sample_to_hit_clones)

        # Columns are positions, rows are samples. 1 = hit at that position.
        # We build it transposed then flip it.
        hits_by_sample_and_position = pandas.DataFrame(
            index=all_coverage_by_clone.columns, dtype=int)
        for sample in all_samples:
            hits_by_sample_and_position[sample] = all_coverage_by_clone.multiply(
                self.clone_by_sample_hits_matrix.reindex(
                    all_coverage_by_clone.index).fillna(0)[sample],
                axis=0).sum()
        hits_by_sample_and_position = hits_by_sample_and_position.T

        fig, axes = pyplot.subplots(
            3 + (1 if heatmap else 0), 1,
            figsize=(10, (len(hits_by_sample_and_position) * 0.15 + 5) if heatmap else 10),
            sharex=True,
            gridspec_kw={
                "height_ratios": [0.3, 0.22, 0.22] + ([0.6] if heatmap else [])
            },
        )
        axes = list(axes)
        features = special_features.get(antigen, [])

        # Sequence
        ax = axes.pop(0)
        pyplot.sca(ax)
        record = dna_features_viewer.GraphicRecord(
            sequence=sequence, features=features)
        record.plot(ax=ax)
        if len(sequence) <= max_len_to_show_sequence:
            record.plot_sequence(ax=ax, fontdict={'size': 6})

        # Coverage
        ax = axes.pop(0)
        pyplot.sca(ax)
        pyplot.plot(
            all_coverage_total.index, all_coverage_total.values, color='black')
        pyplot.fill_between(
            all_coverage_total.index,
            all_coverage_total.values,
            color='lightblue',
            alpha=0.3)
        pyplot.title("   Coverage", loc="left", pad=0)
        pyplot.ylim(ymin=0)
        pyplot.ylabel("Clones")

        # Fraction of samples with a hit
        ax = axes.pop(0)
        pyplot.sca(ax)
        plot_series = (hits_by_sample_and_position > 0).mean(0) * 100.0
        pyplot.plot(
            plot_series.index, plot_series.values, color='black')
        pyplot.fill_between(
            plot_series.index,
            plot_series.values,
            color='red',
            alpha=0.3)
        pyplot.title("   Hits", loc="left", pad=0)
        pyplot.ylim(ymin=0)
        pyplot.ylabel("Samples (%)")

        # Heatmap
        if heatmap:
            ax = axes.pop(0)
            pyplot.sca(ax)
            plot_df = hits_by_sample_and_position.copy()
            plot_df = plot_df.astype(float)

            seaborn.heatmap(plot_df, cbar=False, ax=ax)

            # Switch to one-based numbering
            ticks = numpy.arange(19, len(sequence), 20)
            while len(ticks) > 20:
                ticks = ticks[::2]
            pyplot.xticks(ticks, ticks + 1, rotation=0)

            yticks = plot_df.index
            pyplot.yticks(numpy.arange(len(yticks)), yticks, fontsize=6)

            pyplot.xlabel("Position")
            pyplot.ylabel("Sample")

        seaborn.despine()

        title = antigen.split("OS")[0].strip()
        title += "\nHits in %d of %d clones" % (
            self.clone_by_sample_hits_matrix.reindex(
                all_coverage_by_clone.index).fillna(0)[all_samples].any(1).sum(),
            len(all_coverage_by_clone))

        pyplot.suptitle(title)

        pyplot.tight_layout()

        #import ipdb ; ipdb.set_trace()

        return fig

