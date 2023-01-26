from dataclasses import dataclass

from bioinformatics_textbook.ch01.ch01 import (
    ba1d, find_starting_positions,
    ba1e, find_clumps,
    ba1f, find_min_skew_positions, define_dna_gc_skews,
    ba1g, compute_hamming_distance,
    ba1h, find_approx_occurrence_positions
)
import click
import pytest


@pytest.fixture
def sample_ba1h(fs):
    @dataclass
    class Sample:
        pattern = "AAAAA"
        text = "AACAAGCTGATAAACATTTAAAGAG"
        num_allowed_mismatches = 1
        approx_occurrence_positions = [0, 9, 11, 19]
        sample_dataset = fs.create_file(
            "sample_dataset.txt",
            contents="ATTCTGGA\nCGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC\n3\n"
        )
        sample_output = "6 7 26 27 78"

    yield Sample()


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
    @dataclass
    class Sample:
        dna_q = "CGGAC"
        hamming_distance = 2
        dna_p = "CGAAT"
        sample_dataset = fs.create_file(
            "sample_dataset.txt",
            contents="GGGCCGTTGGT\nGGACCGTTGAC\n"
        )
        sample_output = 3

    yield Sample()


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


@pytest.fixture
def sample_ba1f(fs):
    @dataclass
    class Sample:
        genome = "CATGGGCATCGGCCATACGCC"
        skews = [int(skew) for skew in "0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2".split()]
        min_skew_positions = [21]

        sample_dataset = fs.create_file(
            "sample_dataset.txt",
            contents = "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG\n"
        )
        sample_output = "53 97"
    
    yield Sample()


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
    @dataclass
    class Sample:
        genome = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC"
        k = 5
        L = 75
        t = 4
        fake_file = fs.create_file(
            "sample_dataset.txt",
            contents="CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC\n5 75 4\n"
        )
        sample_output = "CGACA GAAGA AATGT"
    
    yield Sample()


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
    @dataclass
    class Sample:
        pattern = "ATAT"
        genome = "GATATATGCATATACTT"
        positions = "1 3 9"
        positions_list = [1, 3, 9]
        fake_file = fs.create_file("file.txt", contents="ATAT\nGATATATGCATATACTT\n")
    
    yield Sample()


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
def sample_text_k():
    text = "ACGTTTCACGTTTTACGG"
    k = 3

    yield text, k
