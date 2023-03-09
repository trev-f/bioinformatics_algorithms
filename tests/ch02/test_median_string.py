from dataclasses import dataclass

import pytest

from bioinformatics_textbook.ch02.median_string import MedianString
from bioinformatics_textbook.dna import DNA


@pytest.fixture
def median_string():
    @dataclass
    class Sample:
        k = 3
        dnas = [
            DNA(dna) for dna in ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTACGGGACAG',]
        ]

        med_string = 'GAC'

    yield Sample()


def test_find_median_string(median_string):
    expected_median_string = median_string.med_string

    actual_median_string = MedianString().find_median_strings(
        kmer_length=median_string.k,
        dnas=median_string.dnas
    )

    assert expected_median_string in actual_median_string
    

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
