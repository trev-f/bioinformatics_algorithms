from dataclasses import dataclass

import pytest

from bioinformatics_textbook.dna import DNA


@pytest.fixture
def sample_reverse_complement():
    @dataclass
    class SampleReverseComplement:
        dna = 'AAAACCCGGT'
        reverse_complement = 'ACCGGGTTTT'
    
    yield SampleReverseComplement()


def test_reverse_complement_dna(sample_reverse_complement):
    expected_rev_comp = sample_reverse_complement.reverse_complement

    actual_rev_comp = DNA(sample_reverse_complement.dna).reverse_complement()

    assert  actual_rev_comp == expected_rev_comp
