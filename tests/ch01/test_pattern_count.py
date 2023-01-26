import pytest

from dataclasses import dataclass

from bioinformatics_textbook.ch01.pattern_count import PatternCount


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

    actual_counts = PatternCount().count_pattern(
        text=pattern_counts.text,
        pattern=pattern_counts.pattern,
    )

    assert actual_counts == expected_counts
