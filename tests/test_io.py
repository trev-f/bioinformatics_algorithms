import click
import pytest
from bioinformatics_textbook.code_challenges.io import (
    read_last_line, read_text_pattern
)


@pytest.fixture
def ba1a_sample_dataset_path():
    return "tests/datasets/ch01/ba1a_sample_dataset.txt"


def test_read_text_pattern(ba1a_sample_dataset_path):
    expected_text = "GCGCG"
    expected_pattern = "GCG"

    with click.open_file(ba1a_sample_dataset_path) as file:
        actual_text, actual_pattern = read_text_pattern(file)

    assert expected_text == actual_text
    assert expected_pattern == actual_pattern


def test_read_last_line(fs):
    # setup fake file
    fake_file = fs.create_file("multiline_test.txt", contents="first line\nlast line")

    expected_last_line = "last line"

    with click.open_file(fake_file.path, "rb") as file:
        actual_last_line = read_last_line(file)

    assert expected_last_line == actual_last_line
