import pytest
from bioinformatics_textbook.code_challenges.io import (
    read_text_pattern
)


@pytest.fixture
def ba1a_sample_dataset_path():
    return "tests/datasets/ch01/ba1a_sample_dataset.txt"


def test_read_text_pattern(ba1a_sample_dataset_path):
    expected_text = "GCGCG"
    expected_pattern = "GCG"

    actual_text, actual_pattern = read_text_pattern(ba1a_sample_dataset_path)

    assert expected_text == actual_text
    assert expected_pattern == actual_pattern
