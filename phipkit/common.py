import time

import numpy
from tqdm import tqdm

LAST_SAY = None


def say(*msgs):
    global LAST_SAY
    if LAST_SAY is None:
        indicator = "[%s]" % time.asctime()
    else:
        indicator = "[%0.2f sec.]" % (time.time() - LAST_SAY)
    print(indicator, *msgs)
    LAST_SAY = time.time()


def reconstruct_antigen_sequences(blast_df):
    """
    Reconstruct protein amino acid sequences from blast results.

    Parameters
    ----------
    blast_df : pandas.DataFrame

    Returns
    -------
    pandas.Series : str -> str
        Amino acid sequences for each protein

    """
    antigen_lengths = blast_df.groupby("title").hit_to.max()
    antigen_sequences = antigen_lengths.map(lambda l: bytearray(l))
    for (idx, hseq, hit_from, hit_to, title) in tqdm(
            blast_df[["hseq", "hit_from", "hit_to", "title"]].itertuples(
                index=True), total=len(blast_df)):
        hseq = hseq.replace("-", "")

        if len(hseq) != hit_to - hit_from + 1:
            raise ValueError(
                "%d != %d. %s" % (
                    len(hseq), hit_to - hit_from + 1, blast_df.loc[idx]))

        antigen_sequences[title][hit_from - 1: hit_to] = hseq.encode('ascii')
    antigen_sequences = antigen_sequences.map(lambda arr: arr.decode())
    return antigen_sequences


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