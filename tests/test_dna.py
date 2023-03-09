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


@pytest.fixture
def sample_all_possible_kmers():
    @dataclass
    class SampleAllPossibleKmers:
        k = 2

        kmers = [
            'AA', 'AC', 'AG', 'AT',
            'CA', 'CC', 'CG', 'CT',
            'GA', 'GC', 'GG', 'GT',
            'TA', 'TC', 'TG', 'TT',
        ]
    
    return SampleAllPossibleKmers


def test_generate_all_possible_kmers(sample_all_possible_kmers):
    expected_kmers = sample_all_possible_kmers.kmers

    actual_kmers = list(DNA.generate_all_possible_kmers(kmer_length=sample_all_possible_kmers.k))

    assert actual_kmers == expected_kmers
