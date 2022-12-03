from bioinformatics_textbook.code_challenges.ch01 import (
    construct_kmer_freq_table, count_pattern
)

def test_construct_kmer_freq_table():
    input_text = "ACGTTTCACGTTTTACGG"
    input_k = 3

    expected_freq_table = {
        "ACG": 3,
        "CGT": 2,
        "GTT": 2,
        "TTT": 3,
        "TTC": 1,
        "TCA": 1,
        "CAC": 1,
        "TTA": 1,
        "TAC": 1,
        "CGG": 1
    }

    actual_freq_table = construct_kmer_freq_table(text=input_text, k = input_k)

    assert expected_freq_table == actual_freq_table


def test_count_pattern():
    input_text = "GCGCG"
    input_pattern = "GCG"

    expected_output = 2
    actual_output = count_pattern(text=input_text, pattern=input_pattern)

    assert expected_output == actual_output
