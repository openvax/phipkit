import time
from numpy.testing import assert_equal
from pandas.testing import assert_frame_equal
import contextlib
import numpy
import pandas
import io

from phipkit.call_antigens import find_consensus, analyze_antigen_for_sample


def test_find_consensus(repeats=1000):
    some_amino_acids = pandas.Series(["A", "S", "C", "D", "E", "F"])
    sequences_df = pandas.DataFrame(
        index=numpy.arange(1000), columns=numpy.arange(10), dtype=object)
    sequences_df[:] = "A"
    sequences_df.loc[:, 8] = "Q"
    sequences_df.loc[:, 1] = some_amino_acids.sample(
        n=len(sequences_df), replace=True).values
    print(sequences_df)
    sequences = [
        "".join(row.values)
        for _, row in sequences_df.iterrows()
    ]
    consensus = find_consensus(sequences)
    print(consensus)
    assert_equal(consensus, "A.AAAAAAQA")

    if repeats:
        start = time.time()
        for i in range(repeats):
            find_consensus(sequences)
        elapsed = time.time() - start
        calls_per_second = repeats / elapsed
        print("Total elapsed for %d calls: %0.5f sec, %0.1f calls/sec." % (
            repeats, elapsed,  calls_per_second))


EXAMPLE_BLAST_CSV_FOR_SINGLE_ANTIGEN = """
program,version,search_target,clone,len,hit_num,id,accession,title,hsp_num,bit_score,score,evalue,identity,positive,query_from,query_to,hit_from,hit_to,align_len,gaps,qseq,hseq,midline
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E_0,60,1,gnl|BL_ORD_ID|5,5,E,1,41.2022,95,1.89719e-11,20,20,1,20,1,20,20,0,MYSFVSEETGTLIVNSVLLF,MYSFVSEETGTLIVNSVLLF,MYSFVSEETGTLIVNSVLLF
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E_4,60,1,gnl|BL_ORD_ID|5,5,E,1,39.2762,90,9.78917e-11,20,20,1,20,5,24,20,0,VSEETGTLIVNSVLLFLAFV,VSEETGTLIVNSVLLFLAFV,VSEETGTLIVNSVLLFLAFV
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E_32,60,1,gnl|BL_ORD_ID|5,5,E,1,41.2022,95,1.81594e-11,20,20,1,20,33,52,20,0,ILTALRLCAYCCNIVNVSLV,ILTALRLCAYCCNIVNVSLV,ILTALRLCAYCCNIVNVSLV
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E_36,60,1,gnl|BL_ORD_ID|5,5,E,1,44.2838,103,1.13837e-12,20,20,1,20,37,56,20,0,LRLCAYCCNIVNVSLVKPSF,LRLCAYCCNIVNVSLVKPSF,LRLCAYCCNIVNVSLVKPSF
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E_40,60,1,gnl|BL_ORD_ID|5,5,E,1,43.5134,101,2.53183e-12,20,20,1,20,41,60,20,0,AYCCNIVNVSLVKPSFYVYS,AYCCNIVNVSLVKPSFYVYS,AYCCNIVNVSLVKPSFYVYS
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E2_32,60,3,gnl|BL_ORD_ID|5,5,E,1,37.7354,86,3.75643e-10,18,18,1,20,33,52,20,0,ILTALRLCAYCCNQVNVKLV,ILTALRLCAYCCNIVNVSLV,ILTALRLCAYCCN VNV LV
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E2_36,60,3,gnl|BL_ORD_ID|5,5,E,1,40.817,94,2.49408e-11,18,18,1,20,37,56,20,0,LRLCAYCCNQVNVKLVKPSF,LRLCAYCCNIVNVSLVKPSF,LRLCAYCCN VNV LVKPSF
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E2_40,60,3,gnl|BL_ORD_ID|5,5,E,1,40.0466,92,4.96867e-11,18,18,1,20,41,60,20,0,AYCCNQVNVKLVKPSFYVYS,AYCCNIVNVSLVKPSFYVYS,AYCCN VNV LVKPSFYVYS
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E3_36,60,2,gnl|BL_ORD_ID|5,5,E,1,42.743,99,5.33019e-12,19,19,1,20,36,55,20,0,ALRLCAYCCNIVNVQLVKPS,ALRLCAYCCNIVNVSLVKPS,ALRLCAYCCNIVNV LVKPS
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E3_40,60,2,gnl|BL_ORD_ID|5,5,E,1,43.1282,100,3.44009e-12,19,19,1,20,40,59,20,0,CAYCCNIVNVQLVKPSFYVY,CAYCCNIVNVSLVKPSFYVY,CAYCCNIVNV LVKPSFYVY
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E4_32,60,2,gnl|BL_ORD_ID|5,5,E,1,40.817,94,2.66331e-11,19,20,1,20,34,53,20,0,LTALRLCAYCCNIVNVNLVK,LTALRLCAYCCNIVNVSLVK,LTALRLCAYCCNIVNV+LVK
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E4_36,60,2,gnl|BL_ORD_ID|5,5,E,1,43.8986,102,1.76405e-12,19,20,1,20,38,57,20,0,RLCAYCCNIVNVNLVKPSFY,RLCAYCCNIVNVSLVKPSFY,RLCAYCCNIVNV+LVKPSFY
blastp,BLASTP 2.6.0+,{'db': 'reference.fna'},E4_40,60,2,gnl|BL_ORD_ID|5,5,E,1,40.4318,93,4.14736e-11,18,19,1,19,42,60,19,0,YCCNIVNVNLVKPSFYVYS,YCCNIVNVSLVKPSFYVYS,YCCNIVNV+LVKPSFYVYS
""".strip()


def test_analyze_antigen_for_sample(repeats=100, profile=False):
    sub_blast_df = pandas.read_csv(
        io.StringIO(EXAMPLE_BLAST_CSV_FOR_SINGLE_ANTIGEN))
    sequence = 'MYSFVSEETGTLIVNSVLLFLAFVXXXXXXXXILTALRLCAYCCNIVNVSLVKPSFYVYS'
    antigen = 'E'

    def call(sub_blast_df):
        return analyze_antigen_for_sample(
            "sample_1",
            antigen,
            sub_blast_df,
            sequence)
    results = call(sub_blast_df)
    print(results)

    expected = pandas.DataFrame([
        {
            'sample_id': 'sample_1',
            'antigen': 'E',
            'start_position': 41,
            'end_position': 52,
            'clone_consensus': 'YCCNIVNV.LV',
            'antigen_sequence': 'YCCNIVNVSLV',
            'antigen_matches_consensus': True,
            'priority_within_antigen': 1,
            'num_clones': 11,
            'clones': [
                'E_32', 'E_36', 'E_40', 'E2_32', 'E2_36', 'E2_40', 'E3_36',
                'E3_40', 'E4_32', 'E4_36', 'E4_40'],
        },
        {
            'sample_id': 'sample_1',
            'antigen': 'E',
            'start_position': 4,
            'end_position': 20,
            'clone_consensus': 'VSEETGTLIVNSVLLF',
            'antigen_sequence': 'VSEETGTLIVNSVLLF',
            'antigen_matches_consensus': True,
            'priority_within_antigen': 2,
            'num_clones': 2,
            'clones': ['E_0', 'E_4'],
        },
    ])

    assert_frame_equal(results, expected)

    if repeats:
        pprofile = None
        profiler = contextlib.nullcontext()
        if profile:
            import pprofile
            profiler = pprofile.Profile()

        # Make it bigger
        bigger_blast_df = sub_blast_df.copy()
        while len(bigger_blast_df) < 500:
            new = bigger_blast_df.copy()
            new.clone += "_B"
            bigger_blast_df = pandas.concat(
                [bigger_blast_df, new], ignore_index=True)

        start = time.time()
        with profiler:
            for i in range(repeats):
                call(bigger_blast_df)
        elapsed = time.time() - start
        if pprofile:
            profiler.print_stats()

        calls_per_second = repeats / elapsed
        print("Total elapsed for %d calls: %0.5f sec, %0.1f calls/sec." % (
            repeats, elapsed,  calls_per_second))

    #import ipdb ; ipdb.set_trace()


