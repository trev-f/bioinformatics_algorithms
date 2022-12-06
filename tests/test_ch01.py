from bioinformatics_textbook.code_challenges.ch01 import (
    ba1c,
    complement_dna, construct_kmer_freq_table, count_pattern, find_frequent_words,
    find_max_val_of_dict, reverse_complement_dna,
    find_starting_positions, format_starting_positions
)
import click
import pytest


@pytest.fixture
def sample_pattern_matching():
    class SamplePatternMatching:
        def __init__(self):
            self.pattern = "ATAT"
            self.genome = "GATATATGCATATACTT"
            self.positions = "1 3 9"
            self.positions_list = [1, 3, 9]
    
    yield SamplePatternMatching()


def test_find_starting_positions(sample_pattern_matching):
    pattern = sample_pattern_matching.pattern
    genome = sample_pattern_matching.genome
    expected_positions = sample_pattern_matching.positions_list

    actual_positions = find_starting_positions(pattern, genome)

    assert expected_positions == actual_positions


def test_format_starting_positions(sample_pattern_matching):
    positions_list = sample_pattern_matching.positions_list
    expected_formatted_positions = sample_pattern_matching.positions

    actual_formatted_positions = format_starting_positions(positions_list)

    assert expected_formatted_positions == actual_formatted_positions


@pytest.fixture
def sample_reverse_complement():
    class SampleReverseComplement:
        def __init__(self):
            self.dna = "AAAACCCGGT"
            self.complement = "TTTTGGGCCA"
            self.reverse_complement = "ACCGGGTTTT"

            self.dataset_path = "tests/datasets/ch01/ba1c_sample_dataset.txt"
    
    yield SampleReverseComplement()


def test_ba1c(sample_reverse_complement):
    input_file = sample_reverse_complement.dataset_path
    expected_rev_comp = sample_reverse_complement.reverse_complement

    with click.open_file(input_file, "r") as file:
        actual_rev_comp = ba1c(file)
    
    assert expected_rev_comp == actual_rev_comp


def test_reverse_complement_dna(sample_reverse_complement):
    dna = sample_reverse_complement.dna
    expected_rev_comp = sample_reverse_complement.reverse_complement

    actual_rev_comp = reverse_complement_dna(dna)

    assert expected_rev_comp == actual_rev_comp


def test_complement_dna(sample_reverse_complement):
    dna = sample_reverse_complement.dna
    expected_complement = sample_reverse_complement.complement

    actual_complement = complement_dna(dna)

    assert expected_complement == actual_complement


@pytest.fixture
def sample_text_k():
    text = "ACGTTTCACGTTTTACGG"
    k = 3

    yield text, k


def test_find_frequent_words(sample_text_k):
    input_text, input_k = sample_text_k
    expected_frequent_words = "ACG TTT"

    actual_frequent_words = find_frequent_words(text=input_text, k=input_k)

    assert expected_frequent_words == actual_frequent_words


def test_find_max_val_of_dict():
    sample_d1 = {"x": 0, "y": 24}
    expected_max_val_d1 = 24
    
    actual_max_val_d1 = find_max_val_of_dict(d=sample_d1)

    assert expected_max_val_d1 == actual_max_val_d1

    sample_d2 = {"a": -16, "y": -33}
    expected_max_val_d2 = -16

    actual_max_val_d2 = find_max_val_of_dict(d=sample_d2)

    assert expected_max_val_d2 == actual_max_val_d2


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
