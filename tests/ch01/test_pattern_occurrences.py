import pytest

from dataclasses import dataclass

from bioinformatics_textbook.ch01.pattern_occurrences import PatternOccurrences


@pytest.fixture
def pattern_counts():
    @dataclass
    class Sample:
        text = "GCGCG"
        pattern = "GCG"
        counts = 2

    yield Sample()


def test_count_pattern(pattern_counts):
    expected_counts = pattern_counts.counts

    actual_counts = PatternOccurrences().count_pattern(
        text=pattern_counts.text,
        pattern=pattern_counts.pattern,
    )

    assert actual_counts == expected_counts


@pytest.fixture
def sample_pattern_matching(fs):
    @dataclass
    class Sample:
        pattern = "ATAT"
        genome = "GATATATGCATATACTT"
        positions_list = [1, 3, 9]
    
    yield Sample()


def test_find_starting_positions(sample_pattern_matching):
    expected_positions = sample_pattern_matching.positions_list

    actual_positions = PatternOccurrences().find_starting_positions(
        pattern=sample_pattern_matching.pattern,
        genome=sample_pattern_matching.genome
    )

    assert actual_positions == expected_positions

