from dataclasses import dataclass

import pytest

from bioinformatics_textbook.ch02.motif import Motif


@pytest.fixture
def k_d_motif():
    @dataclass
    class Sample:
        k = 3
        d = 1
        dnas = [
            'ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'
        ]

        k_d_motifs = set('ATA ATT GTT TTT'.split(' '))
    
    yield Sample()


def test_find_k_d_motifs(k_d_motif):
    expected_motifs = k_d_motif.k_d_motifs

    actual_motifs = Motif().find_k_d_motifs(
        kmer_length=k_d_motif.k,
        num_allowed_mismatches=k_d_motif.d,
        dnas=k_d_motif.dnas
    )

    assert actual_motifs == expected_motifs
