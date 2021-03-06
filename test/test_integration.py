import random
import tempfile
import os

import numpy
import pandas

from phipkit.blast import run_blast
from phipkit.score import compute_scores
from phipkit.call_hits import call_hits
from phipkit.call_antigens import call_antigens, hits_to_dict
from phipkit.plot_antigens import plot_antigens

random.seed(42)
numpy.random.seed(42)

TEST_PROTEINS = {
    ##### Simple cases
    "ORF6": "MFHLVDFQVTIAEILLIIMRTFKVSIWNLDYIINLIIKNLSKSLTENKYSQLDEEQPMEID",
    "ORF7b": "MIELSLIDFYLCFLAFLLFLVLIMLIIFWFSLELQDHNETCHA",
    "ORF8": "MKFLVFLGIITTVAAFHQECSLQSCTQHQPYVVDDPCPIHFYSKWYIRVGARKSAPL"
        "IELCVDEAGSKSPIQYIDIGNYTVSCLPFTINCQEPKLGSLVVRCSFYEDFLEYHDVRVVLDFI",
    "M": "MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLWLLWPV"
        "TLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMWSFNPETNILL"
        "NVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLPKEITVATSRTLSYYK"
        "LGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALLVQ",
    "S": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS"
        "NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV"
        "NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE"
        "GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT"
        "LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK"
        "CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN"
        "CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD"
        "YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC"
        "NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
        "FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP"
        "GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY"
        "ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI"
        "SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE"
        "VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC"
        "LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM"
        "QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN"
        "TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA"
        "SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA"
        "ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP"
        "LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL"
        "QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD"
        "SEPVLKGVKLHYT",

    #
    # Case 1. Clones blast to multiple antigens, but only one of the antigens is
    # actually responsible for the hits, and these can be disambiguated.
    #
    # Case 2. Cannot be disambiguated
    #
    # E with the actual hit TGTLIV
    #             CASE 1                                CASE 2
    #             *****                                 ******
    "E": "MYSFVSEETGTLIVNSVLLFLAFVVFLLVTLAILTALRLCAYCCNIVNVSLVKPSFYVYS",

    # Mutation at epitope
    #              ****!*                                ******
    "E2": "MYSFQSEQTGTLVVNSVLLFLAFKVFLLQTLAILTALRLCAYCCNQVNVKLVKPSFYVYS",

    # Insertion at epitope
    #              ****!**                                ******
    "E3": "MYSFVSEETGTLIIVNSVLLFLAFVVFLLVTLAILTALRLCAYCCNIVNVQLVKPSFYVYS",

    # Deletion at epitope
    #              ****!*
    "E4": "MYSFVSEETGTLVNSVLLFLAFVVFLLVTLAILTALRLCAYCCNIVNVNLVKPSFYVYS",

}
TEST_PROTEINS["S_redundant"] = TEST_PROTEINS["S"]


def test_integrated(save_dir=None):
    """
    Run the integration test.

    If save_dir is not None, then the input files are written to this directory.
    This can be used as a source of sample data for testing the commands.
    """
    if save_dir:
        # Write test proteins
        with open(os.path.join(save_dir, "proteins.fasta"), "w") as fd:
            for (name, sequence) in TEST_PROTEINS.items():
                fd.write(">%s\n" % name)
                fd.write(sequence)
                fd.write("\n")

    epitopes_dict = {
        "sample_a": [
            "YICGDSTE",  # S
            "LIMLIIFWF", # ORF7b
        ],
        "sample_b": [
            "TGTLIV",  # E
            "VNV.LV",  # E, E1-E4
            "SPECIAL", # not templated
        ],
    }

    overlap = 4
    length = 20
    clone_sequences = []
    for (protein, seq) in TEST_PROTEINS.items():
        start = 0
        while start < len(seq) - 1:
            take = seq[start : start + length]
            clone_sequences.append(("%s_%d" % (protein, start), take))
            start += overlap

    # Test that we can call an epitope in the middle of an alignment
    clone_sequences.extend([
        ("SPECIAL_1", "KTCPVQLWVDSTPPPGTRwqspecialeeLDGEYFTLQIRGRERFEMF".upper()),
        ("SPECIAL_2", "EEPQSDPSVEPPLSQETqwspecialeeDRRTEEENLR".upper()),
        ("SPECIAL_3", "LELKDAQqwspecialeeAHSSHLKSKK".upper()),
    ])

    annot_df = pandas.DataFrame(
        clone_sequences, columns=["name", "seq_aa"]).set_index("name")

    print("Generated sequences")
    print(annot_df)
    if save_dir:
        annot_df.to_csv(os.path.join(save_dir, "annotations.csv"))

    pairwise_blast_df = run_blast(
        annot_df,
        extra_blast_parameters={
            "word_size": "5",
            "evalue": "2",
            "max_target_seqs": "5",
        })
    print("Generated pairwise blast results:")
    print(pairwise_blast_df)

    pairwise_blast_without_self_df = pairwise_blast_df.loc[
        pairwise_blast_df.clone != pairwise_blast_df.title]
    assert (
        pairwise_blast_without_self_df.loc[
            pairwise_blast_without_self_df.clone.str.contains("SPECIAL")
        ].title.nunique()) > 2

    with tempfile.NamedTemporaryFile(mode="w") as fd:
        for (key, value) in TEST_PROTEINS.items():
            fd.file.write(">%s\n%s\n" % (key, value))
        fd.file.close()
        reference_blast_df = run_blast(
            annot_df,
            reference_fasta_path=fd.name,
            extra_blast_parameters={
                "word_size": "5",
            })
        print("Generated reference blast:")
        print(reference_blast_df)

    counts_df = pandas.DataFrame(index=annot_df.index)

    beads_lambda = numpy.random.randint(0, 10, len(counts_df))
    for i in range(8):
        counts_df["beads_%d" % i] = numpy.random.poisson(
            beads_lambda, len(counts_df)) * numpy.random.uniform(1.0, 10.0)
    annot_df["beads_lambda"] = beads_lambda

    counts_df["sample_a"] = numpy.random.poisson(
        beads_lambda * 5, len(counts_df))

    counts_df["sample_b"] = numpy.random.poisson(
        beads_lambda * 0.5, len(counts_df))

    epitopes_df = []
    for (sample, sequences) in epitopes_dict.items():
        for epitope in sequences:
            clones = annot_df.loc[annot_df.seq_aa.str.contains(epitope)].index
            print(sample, epitope, *clones)
            counts_df.loc[clones, sample] += 20
            counts_df.loc[clones, sample] *= 25
            print("Setting hit", sample, epitope, *clones)
            for clone in clones:
                epitopes_df.append((sample, epitope, clone))

    # Helps call the correct TGTLIV epitope. Would be great to remove this line.
    counts_df.loc["E_12", "sample_b"] = 0.0

    counts_df = counts_df.astype(int)

    if save_dir:
        counts_df.to_csv(os.path.join(save_dir, "counts.csv"))

    epitopes_df = pandas.DataFrame(
        epitopes_df, columns=["sample_id", "epitope", "clone"])
    sample_to_expected_hits = epitopes_df.groupby("sample_id").clone.unique()

    scores_df = compute_scores(counts_df=counts_df)
    print("Scores")
    print(scores_df)

    hits_df = call_hits(
        pairwise_blast_df, scores_df=scores_df, counts_df=counts_df, fdr=0.1)

    hits_df["expected"] = [
        row.clone1 in sample_to_expected_hits.get(row.sample_id, []) and
        row.clone2 in sample_to_expected_hits.get(row.sample_id, [])
        for _, row in hits_df.iterrows()
    ]
    print(hits_df)
    print("False discoveries:", (~hits_df.expected).sum())
    print("False discovery rate", 1 - hits_df.expected.mean())
    assert 1 - hits_df.expected.mean() < 0.5

    epitopes_df["discovered"] = [
        len(hits_df.loc[
            (hits_df.sample_id == row.sample_id) &
            (
                (hits_df.clone1 == row.clone) |
                (hits_df.clone2 == row.clone)
            )
        ]) > 0
        for _, row in epitopes_df.iterrows()
    ]
    print("Epitopes:")
    print(epitopes_df)
    print("Failed to discover", (~epitopes_df.discovered).sum())
    print("Failed to discover rate", (~epitopes_df.discovered).mean())
    assert (~epitopes_df.discovered).mean() < 0.5

    print("Calling antigens")
    antigens_df = call_antigens(reference_blast_df, hits_to_dict(hits_df))
    print(antigens_df)

    antigens_df["expected_hits"] = [
        [clone for clone in row.clones
            if clone in sample_to_expected_hits.get(row.sample_id, [])]
        for _, row in antigens_df.iterrows()
    ]
    antigens_df["unexpected_hits"] = [
        [clone for clone in row.clones
            if clone not in sample_to_expected_hits.get(row.sample_id, [])]
        for _, row in antigens_df.iterrows()
    ]

    antigens_df["expected"] = [
        " ".join([
            epitope
            for epitope in epitopes_dict.get(row.sample_id, [])
            if epitope in row.clone_consensus
        ])
        for _, row in antigens_df.iterrows()
    ]

    quality_antigens_df = antigens_df.loc[antigens_df.antigen_matches_consensus]
    fdr = (quality_antigens_df.expected == "").mean()
    print("Quality antigens:")
    print(quality_antigens_df)
    assert fdr < 0.35

    # Check if we found our minimal epitopes. Exclude the special ones, because
    # they don't have a corresponding protein.
    minimal_epitopes_df = pandas.DataFrame([
        (sample_id, clone)
        for (sample_id, clones) in epitopes_dict.items()
        for clone in clones
    ], columns=["sample_id", "epitope"])
    minimal_epitopes_df = minimal_epitopes_df.loc[
        ~minimal_epitopes_df.epitope.str.contains("SPECIAL")
    ]

    discovered_epitopes = quality_antigens_df.groupby(
        "sample_id").expected.unique()
    minimal_epitopes_df["discovered"] = [
        row.epitope in discovered_epitopes[row.sample_id]
        for _, row in minimal_epitopes_df.iterrows()
    ]
    print(minimal_epitopes_df)
    assert minimal_epitopes_df.discovered.mean() > 0.7

    # Check that S and S_redundant are called as non-redundant and redundant
    assert antigens_df.loc[antigens_df.antigen == "S"].redundant.unique() == [False]
    assert antigens_df.loc[antigens_df.antigen == "S_redundant"].redundant.unique() == [True]

    # Call using the pairwise reference and check that we recover the special
    # epitope.
    print("Calling antigens using pairwise reference")
    antigens_pairwise_df = call_antigens(
        pairwise_blast_df, hits_to_dict(hits_df))
    print(antigens_pairwise_df)
    sample_b = antigens_pairwise_df.loc[
        (antigens_pairwise_df.sample_id == "sample_b")
    ]
    assert "SPECIAL_1" in sample_b.antigen.unique()
    assert "SPECIAL_2" in sample_b.antigen.unique()
    assert "SPECIAL_3" in sample_b.antigen.unique()

    # Test plotting of antigens
    out = "/tmp/phipkit.antigen_plots.pdf"
    if os.path.exists(out):
        os.unlink(out)
    summary_df = plot_antigens(
        blast_df=reference_blast_df,
        hits_df=hits_df,
        antigens_df=antigens_df,
        out=out)
    assert os.path.exists(out)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--save-dir",
        metavar="DIR",
        help="Write input files to DIR.")

    args = parser.parse_args(sys.argv[1:])
    test_integrated(save_dir=args.save_dir)
