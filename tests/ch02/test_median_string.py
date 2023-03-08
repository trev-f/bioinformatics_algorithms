from dataclasses import dataclass

import pytest

from bioinformatics_textbook.ch02.median_string import MedianString
from bioinformatics_textbook.dna import DNA


@pytest.fixture
def pattern_strings_distance():
    @dataclass
    class Sample:
        pattern = DNA('AAA')
        dnas = [
            DNA(dna) for dna in ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT',]
        ]

        distance = 5
    
    yield Sample()


def test_compute_pattern_strings_distance(pattern_strings_distance):
    expected_distance = pattern_strings_distance.distance

    actual_distance = MedianString().compute_pattern_strings_distance(
        pattern=pattern_strings_distance.pattern,
        dnas=pattern_strings_distance.dnas
    )

    assert actual_distance == expected_distance
