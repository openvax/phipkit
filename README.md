# phipkit
PhIP-seq data analysis toolkit

This package implements a data analysis workflow for processing the results of
phage immunoprecipitation sequencing ([PhIP-Seq](https://en.wikipedia.org/wiki/PhIP-Seq]))
experiments.

This workflow gets you from a matrix
of read counts to a list of potential epitopes for each sample. Downstream
analyses (not implemented here) may use these epitopes in case vs.
control comparisons or other analyses appropriate to the study.



## Installation

All dependencies except BLAST (required by `phipkit-blast` command) can be
installed by pip. Run this from a checkout of the repo:

```
$ pip install .
```

In order to run the `phipkit-blast` command, you will need the BLAST commands
`blastp` and `makeblastdb` in your PATH. I recommend installing BLAST using
bioconda. After [setting up bioconda](https://bioconda.github.io/user/install.html#set-up-channels),
you can just run:

```
conda install blast
```

This package has been tested using BLAST versions 2.6.0 and 2.10.1. Most likely
any recent version will work.

You can run the phipkit unit tests with:

```
$ pip install nose
$ nosetests -sv test/
```


## Assumptions and data requirements

Key assumptions for our analysis approach:
* Your PhIP-seq run included several (ideally 5 or more)
"beads-only" negative control wells, i.e. wells without any antibody added.
* Your PhIP-seq library was designed by tiling source proteins
with *overlapping* peptides, i.e. each position in the source protein is covered
by more than one phage clone  

To use this package, you will need three things:
* A library annotation file (CSV) giving the amino acid insert sequence (peptide) for each clone
* A reference proteome (FASTA) giving the amino acid sequences of the full-length
proteins that the library covers 
* A sample x clones matrix (CSV) of read counts for your experiment. This can be generated
for example using the [phip-stat](https://github.com/lasersonlab/phip-stat) package.


Note: in most cases where a CSV file is expected, you can also also pass a
tab-separated file (TSV), or a bzip2-compressed CSV or TSV. 


## Example analysis

We will walk through analyzing a toy synthetic PhIP-seq dataset. First go
into the example-data directory:

```
$ cd example-data
```

### Step 1. Blast PhIP-seq library against itself

You only need to do this once for each PhIP-seq library you run. If you mix
libraries in a single run, I recommend running this command on the combined
mixture of libraries.

To blast phage clone sequences against themselves, run:

```
$ phipkit-blast annotations.csv --out blast.all_against_all.csv
```

Example results: [blast.all_against_all.csv](example-data/blast.all_against_all.csv).

The first argument is your library annotation file. This is expected to
have phage clone names on the first column, and another column named `seq_aa`
that gives the amino acid sequence of the peptide. This column name can be
changed; see `phipkit-blast -h`.

Running this on a real library can be quite slow. I usually run it overnight.
Passing '--blast-parameter num_threads 4' can speed this up (the default is 1
thread).


### Step 2. Blast PhIP-seq library against reference proteome

There's a second blast step that you also need to do only once per library. Here
we blast the library against a set of full length protein sequences (referred
to as "antigens" later).

```
$ phipkit-blast annotations.csv --reference proteins.fasta --out blast.reference.csv
```

Example results: [blast.reference.csv](example-data/blast.reference.csv).

The proteins you use here should roughly correspond to the 
proteins used originally to generate the PhIP-seq library. They don't have to
exactly match though. For example, for an autoantigen PhIP-seq library,
it's reasonable to use the current
[human proteome](https://www.uniprot.org/proteomes/UP000005640) from Uniprot.


### Step 3. Calculate scores

This step is independent of the two BLAST steps. It calculates "scores" that
capture enrichment over beads-only negative controls. The results are conceptually
analogous to the z-scores introduced in [Yuan et al 2018](https://www.biorxiv.org/content/10.1101/285916v1),
although a different model is used based on a robust fit of a Poisson
distribution. See [score.py](phipkit/score.py) for the implementation. 

```
$ phipkit-score counts.csv --out scores.csv
```

Example results: [scores.csv](example-data/scores.csv).


### Step 4. Call hits

This step identifies phage clones (actually pairs of phage clones with
overlapping sequences) that have high enough scores to be called hits.

Our method for hit calling uses only the relative ranking of the scores, not the
values themselves. You can therefore use any alternate scoring method you like,
as long as clones with higher read counts have higher scores. The method is
implemented in [call_hits.py](phipkit/call_hits.py).

```
$ phipkit-call-hits blast.all_against_all.csv scores.csv --counts counts.csv --fdr 0.1 --out hits.csv
```

Example results: [hits.csv](example-data/hits.csv).


### Step 5. Call antigens

This step looks for regions of source proteins that have multiple overlapping
hits. It's implemented in [call_antigens.py](phipkit/call_antigens.py).

```
$ phipkit-call-antigens blast.reference.csv hits.csv --out antigens.csv
```

Example results: [antigens.csv](example-data/antigens.csv).


### Step 6 (optional). Plot antigen calls

Some rudimentary plotting support is also available.

```
$ phipkit-plot-antigens blast.reference.csv hits.csv antigens.csv --out antigens.pdf
```

Example results: [antigens.pdf](example-data/antigens.pdf).

