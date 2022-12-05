import click
import pytest
from bioinformatics_textbook.code_challenges.inout import (
    read_all_lines, read_not_last_line, read_last_line, read_text_pattern
)


@pytest.fixture
def ba1a_sample_dataset_path():
    return "tests/datasets/ch01/ba1a_sample_dataset.txt"


def test_read_text_pattern(ba1a_sample_dataset_path):
    expected_text = "GCGCG"
    expected_pattern = "GCG"

    with click.open_file(ba1a_sample_dataset_path, "rb") as file:
        actual_text, actual_pattern = read_text_pattern(file)

    assert expected_text == actual_text
    assert expected_pattern == actual_pattern


@pytest.fixture
def multiline_file(fs):
    fake_file = fs.create_file("multiline_test.txt", contents="first line\nlast line")

    yield fake_file


def test_read_all_lines(multiline_file):
    expected_all_lines = "first linelast line"

    with click.open_file(multiline_file.path, "r") as file:
        actual_all_lines = read_all_lines(file)
    
    assert expected_all_lines == actual_all_lines


def test_read_not_last_line(multiline_file):
    expected_not_last_line = "first line"

    with click.open_file(multiline_file.path, "rb") as file:
        actual_not_last_line = read_not_last_line(file)

    assert expected_not_last_line == actual_not_last_line


def test_read_last_line(multiline_file):
    expected_last_line = "last line"

    with click.open_file(multiline_file.path, "rb") as file:
        actual_last_line = read_last_line(file)

    assert expected_last_line == actual_last_line
