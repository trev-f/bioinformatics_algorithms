import pytest

from dataclasses import dataclass

from bioinformatics_textbook.ch01.frequent_words import FrequentWords


@pytest.fixture
def freq_words():
    @dataclass
    class Sample:
        text = "ACGTTTCACGTTTTACGG"
        kmer_length = 3
        most_freq_words = ["ACG", "TTT"]

    yield Sample()
       

def test_find_most_freq_words(freq_words):
    expected_most_freq_words = freq_words.most_freq_words

    actual_most_freq_words = FrequentWords().find_most_freq_words(
        text=freq_words.text,
        kmer_length=freq_words.kmer_length
    )

    assert actual_most_freq_words == expected_most_freq_words


@pytest.fixture
def freq_words_mismatches():
    @dataclass
    class Sample:
        text = "AACAAGCTGATAAACATTTAAAGAG"
        kmer_length = 5
        num_allowed_mismatches = 1
        most_freq_kmers = ["AAAAA"]

    yield Sample()


def test_find_most_freq_words_with_mismatches(freq_words_mismatches):
    expected_most_freq_kmers = freq_words_mismatches.most_freq_kmers

    actual_most_freq_kmers = FrequentWords().find_most_freq_words_with_mismatches(
        text=freq_words_mismatches.text,
        kmer_length=freq_words_mismatches.kmer_length,
        num_allowed_mismatches=freq_words_mismatches.num_allowed_mismatches
    )

    assert expected_most_freq_kmers == actual_most_freq_kmers
