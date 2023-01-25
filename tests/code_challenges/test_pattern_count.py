import pytest
from bioinformatics_textbook.code_challenges.pattern_count import PatternCount


@pytest.fixture
def pattern_counts():
    class CountPattern:
        def __init__(self):
            self.text = "GCGCG"
            self.pattern = "GCG"

            self.counts = 2

    
    yield CountPattern()


def test_count_pattern(pattern_counts):
    expected_counts = pattern_counts.counts

    actual_counts = PatternCount().count_pattern(
        text=pattern_counts.text,
        pattern=pattern_counts.pattern,
    )

    assert actual_counts == expected_counts
