import numpy
import pandas

from phipkit.score import compute_scores


def test_basic():
    counts_df = pandas.DataFrame(index=["clone_%d" % i for i in range(1000)])

    beads_lambda = numpy.random.randint(0, 10, len(counts_df))**2
    for i in range(8):
        counts_df["beads_%d" % i] = numpy.random.poisson(
            beads_lambda, len(counts_df))

    counts_df["beads_1"] *= 1000
    counts_df.loc["clone_2", "beads_2"] = 1e6
    counts_df.loc["clone_2", "beads_3"] = 1e6

    counts_df["beads_uninformative_0"] = 0

    counts_df["sample_a"] = numpy.random.poisson(
        beads_lambda * 5, len(counts_df))

    counts_df["sample_b"] = numpy.random.poisson(
        beads_lambda * 0.5, len(counts_df))

    counts_df.loc["clone_0", "sample_a"] = 1000
    counts_df.loc["clone_1", "sample_b"] = 1000
    counts_df.loc["clone_2", "sample_b"] = 1e9
    counts_df.loc["clone_bogus"] = 1e10

    print("COUNTS")
    print(counts_df)

    results = compute_scores(counts_df=counts_df)
    print("RESULTS")
    print(results)

    for col in results:
        if col.startswith("sample_"):
            print("Sorted values for %s" % col)
            sort_series = results[col].sort_values(ascending=False)
            print(sort_series)
            print("Counts:")
            print(counts_df.loc[sort_series.head(5).index])

    max_scores = results.idxmax()
    assert max_scores["sample_a"] == "clone_0"
    assert max_scores["sample_b"] == "clone_1"