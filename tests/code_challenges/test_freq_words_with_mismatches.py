import pytest
from bioinformatics_textbook.code_challenges.freq_words_with_mismatches import FrequentWords


@pytest.fixture
def freq_words_mismatches():
    class FreqWordsMismatches:
        def __init__(self):
            self.text = "AACAAGCTGATAAACATTTAAAGAG"
            self.kmer_length = 5
            self.num_allowed_mismatches = 1

            self.most_freq_kmers = ["AAAAA"]
    

    yield FreqWordsMismatches()


def test_find_most_freq_words_with_mismatches(freq_words_mismatches):
    expected_most_freq_kmers = freq_words_mismatches.most_freq_kmers

    actual_most_freq_kmers = FrequentWords().find_most_freq_words_with_mismatches(
        text=freq_words_mismatches.text,
        kmer_length=freq_words_mismatches.kmer_length,
        num_allowed_mismatches=freq_words_mismatches.num_allowed_mismatches
    )

    assert expected_most_freq_kmers == actual_most_freq_kmers
