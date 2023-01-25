import pytest
from dataclasses import dataclass

from bioinformatics_textbook.ch01.reverse_complement import ReverseComplement

@pytest.fixture
def sample_reverse_complement():
    @dataclass
    class SampleReverseComplement:
        dna: str = "AAAACCCGGT"
        reverse_complement: str = "ACCGGGTTTT"

    yield SampleReverseComplement()


def test_reverse_complement_dna(sample_reverse_complement):
    expected_rev_comp = sample_reverse_complement.reverse_complement

    actual_rev_comp = ReverseComplement().reverse_complement_dna(
        dna=sample_reverse_complement.dna
    )

    assert  actual_rev_comp == expected_rev_comp