from bioinformatics_textbook.code_challenges.ch01 import (
    ba1c,
    complement_dna, construct_kmer_freq_table, count_pattern, find_frequent_words,
    find_max_val_of_dict, reverse_complement_dna,
    ba1d, find_starting_positions,
    ba1e, find_clumps,
    ba1f, find_min_skew_positions, define_dna_gc_skews,
    ba1g, compute_hamming_distance, is_mismatch,
    ba1h, find_approx_occurrence_positions,
    convert_iterable_to_list_of_str, format_list_for_rosalind
)
import click
import pytest


@pytest.fixture
def sample_ba1h(fs):
    class SampleBA1H:
        def __init__(self):
            self.pattern = "AAAAA"
            self.text = "AACAAGCTGATAAACATTTAAAGAG"
            self.num_allowed_mismatches = 1
            self.approx_occurrence_positions = [0, 9, 11, 19]

            self.sample_dataset = fs.create_file(
                "sample_dataset.txt",
                contents="ATTCTGGA\nCGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC\n3\n"
            )
            self.sample_output = "6 7 26 27 78"


    yield SampleBA1H()


def test_ba1h(sample_ba1h):
    input_file = sample_ba1h.sample_dataset.path
    expected_approx_occurrence_positions = sample_ba1h.sample_output

    with click.open_file(input_file, "rb") as file:
        actual_approx_occurrence_positions = ba1h(file)
    
    assert expected_approx_occurrence_positions == actual_approx_occurrence_positions


def test_find_approx_occurrence_positions(sample_ba1h):
    pattern = sample_ba1h.pattern
    text = sample_ba1h.text
    num_allowed_mismatches = sample_ba1h.num_allowed_mismatches
    expected_approx_occurrence_positions = sample_ba1h.approx_occurrence_positions

    actual_approx_occurrence_positions = find_approx_occurrence_positions(pattern, text, num_allowed_mismatches)

    assert expected_approx_occurrence_positions == actual_approx_occurrence_positions


@pytest.fixture
def sample_ba1g(fs):
    class SampleBA1G:
        def __init__(self):
            self.dna_p = "CGAAT"
            self.dna_q = "CGGAC"
            self.hamming_distance = 2

            self.sample_dataset = fs.create_file(
                "sample_dataset.txt",
                contents="GGGCCGTTGGT\nGGACCGTTGAC\n"
            )
            self.sample_output = 3


    yield SampleBA1G()


def test_ba1g(sample_ba1g):
    input_file = sample_ba1g.sample_dataset.path
    expected_hamming_distance = sample_ba1g.sample_output

    with click.open_file(input_file, "rb") as file:
        actual_hamming_distsance = ba1g(file)
    
    assert expected_hamming_distance == actual_hamming_distsance


def test_compute_hamming_distance(sample_ba1g):
    dna_p = sample_ba1g.dna_p
    dna_q = sample_ba1g.dna_q
    expected_hamming_distance = sample_ba1g.hamming_distance

    actual_hamming_distance = compute_hamming_distance(dna_p, dna_q)

    assert expected_hamming_distance == actual_hamming_distance


def test_is_mismatch():
    # ensure matches are not called as mismatches
    assert False == is_mismatch("A", "A")
    assert False == is_mismatch("C", "C")
    assert False == is_mismatch("G", "G")
    assert False == is_mismatch("T", "T")

    # test mismatches
    assert True == is_mismatch("A", "C")
    assert True == is_mismatch("A", "G")
    assert True == is_mismatch("A", "T")
    assert True == is_mismatch("C", "A")
    assert True == is_mismatch("C", "G")
    assert True == is_mismatch("C", "T")
    assert True == is_mismatch("G", "C")
    assert True == is_mismatch("G", "A")
    assert True == is_mismatch("G", "T")
    assert True == is_mismatch("T", "C")
    assert True == is_mismatch("T", "G")
    assert True == is_mismatch("T", "A")


@pytest.fixture
def sample_ba1f(fs):
    class SampleBA1F:
        def __init__(self):
            self.genome = "CATGGGCATCGGCCATACGCC"
            self.skews = [int(skew) for skew in "0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2".split()]
            self.min_skew_positions = [21]

            self.sample_dataset = fs.create_file(
                "sample_dataset.txt",
                contents = "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG\n"
            )
            self.sample_output = "53 97"

    
    yield SampleBA1F()


def test_ba1f(sample_ba1f):
    input_file = sample_ba1f.sample_dataset.path
    expected_skew_positions = sample_ba1f.sample_output

    with click.open_file(input_file, "r") as file:
        actual_skew_positions = ba1f(file)

    assert expected_skew_positions == actual_skew_positions


def test_find_min_skew_positions(sample_ba1f):
    skews = sample_ba1f.skews
    expected_min_skew_positions = sample_ba1f.min_skew_positions

    actual_min_skew_positions = find_min_skew_positions(skews)

    assert expected_min_skew_positions == actual_min_skew_positions


def test_define_dna_gc_skews(sample_ba1f):
    genome = sample_ba1f.genome
    expected_skews = sample_ba1f.skews

    actual_skews = define_dna_gc_skews(genome)

    assert expected_skews == actual_skews


@pytest.fixture
def sample_ba1e(fs):
    class SampleBA1E:
        def __init__(self):
            self.genome = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC"
            self.k = 5
            self.L = 75
            self.t = 4

            self.fake_file = fs.create_file(
                "sample_dataset.txt",
                contents="CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC\n5 75 4\n"
            )

            self.sample_output = "CGACA GAAGA AATGT"
    
    yield SampleBA1E()


def test_ba1e(sample_ba1e):
    input_file = sample_ba1e.fake_file
    expected_clump_patterns = sample_ba1e.sample_output

    with click.open_file(input_file.path, "rb") as file:
        actual_clump_patterns = ba1e(file)

    assert expected_clump_patterns == actual_clump_patterns


def test_find_clumps(sample_ba1e):
    genome = sample_ba1e.genome
    k = sample_ba1e.k
    L = sample_ba1e.L
    t = sample_ba1e.t
    expected_clump_patterns = sample_ba1e.sample_output

    actual_clump_patterns = find_clumps(genome, k, L, t)

    assert expected_clump_patterns == actual_clump_patterns


@pytest.fixture
def sample_pattern_matching(fs):
    class SamplePatternMatching:
        def __init__(self):
            self.pattern = "ATAT"
            self.genome = "GATATATGCATATACTT"
            self.positions = "1 3 9"
            self.positions_list = [1, 3, 9]

            self.fake_file = fs.create_file("file.txt", contents="ATAT\nGATATATGCATATACTT\n")
    
    yield SamplePatternMatching()


def test_ba1d(sample_pattern_matching):
    input_file = sample_pattern_matching.fake_file
    expected_starting_positions = sample_pattern_matching.positions

    with click.open_file(input_file.path, "rb") as file:
        actual_starting_positions = ba1d(file)

    assert expected_starting_positions == actual_starting_positions


def test_find_starting_positions(sample_pattern_matching):
    pattern = sample_pattern_matching.pattern
    genome = sample_pattern_matching.genome
    expected_positions = sample_pattern_matching.positions_list

    actual_positions = find_starting_positions(pattern, genome)

    assert expected_positions == actual_positions


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


def test_format_list_for_rosalind():
    list_of_strings = ["AA", "GG"]
    expected_formatted_strings = "AA GG"

    actual_formatted_strings = format_list_for_rosalind(list_of_strings)

    assert expected_formatted_strings == actual_formatted_strings

    list_of_integers = [1, 2]
    expected_formatted_integers = "1 2"

    actual_formatted_integers = format_list_for_rosalind(list_of_integers)

    assert expected_formatted_integers == actual_formatted_integers


def test_convert_iterable_to_list_of_str():
    list_of_integers = [1, 2]
    expected_converted_integers = ["1", "2"]

    actual_converted_integers = convert_iterable_to_list_of_str(list_of_integers)

    assert expected_converted_integers == actual_converted_integers
