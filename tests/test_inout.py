import click
import pytest
from bioinformatics_textbook.code_challenges.inout import (
    read_all_lines, read_not_last_line, read_first_line, read_last_line, read_text_pattern,
    strip_newlines
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
def multiline_files(fs):
    class MultiLineFiles:
        def __init__(self):
            self.lf_file = fs.create_file("lf_test.txt", contents="first line\nlast line\n")
            self.cr_lf_file = fs.create_file("cr_lf_test.txt", contents="first line\r\nlast line\r\n")

            self.three_line_file = fs.create_file("three_line_file.txt", contents="first line\nsecond line\nlast line\n")


    yield MultiLineFiles()


def test_read_all_lines(multiline_files):
    lf_file = multiline_files.lf_file.path
    cr_lf_file = multiline_files.cr_lf_file.path
    expected_all_lines = "first linelast line"

    with click.open_file(lf_file, "r") as file:
        actual_all_lines_lf = read_all_lines(file)
    
    assert expected_all_lines == actual_all_lines_lf

    with click.open_file(cr_lf_file, "r") as file:
        actual_all_lines_cr_lf = read_all_lines(file)

    assert expected_all_lines == actual_all_lines_cr_lf


def test_read_not_last_line(multiline_files):
    lf_file = multiline_files.lf_file.path
    cr_lf_file = multiline_files.cr_lf_file.path
    expected_not_last_line = "first line"

    with click.open_file(lf_file, "rb") as file:
        actual_not_last_line_lf = read_not_last_line(file)

    assert expected_not_last_line == actual_not_last_line_lf

    with click.open_file(cr_lf_file, "rb") as file:
        actual_not_last_line_cr_lf = read_not_last_line(file)

    assert expected_not_last_line == actual_not_last_line_cr_lf


def test_read_first_line(multiline_files):
    input_file = multiline_files.three_line_file.path
    expected_first_line = "first line"

    with click.open_file(input_file, "rb") as file:
        actual_first_line = read_first_line(file)
    
    assert expected_first_line == actual_first_line


def test_read_last_line(multiline_files):
    lf_file = multiline_files.lf_file.path
    cr_lf_file = multiline_files.cr_lf_file.path
    expected_last_line = "last line"

    with click.open_file(lf_file, "rb") as file:
        actual_last_line_lf = read_last_line(file)

    assert expected_last_line == actual_last_line_lf

    with click.open_file(cr_lf_file, "rb") as file:
        actual_last_line_cf_lf = read_last_line(file)
    
    assert expected_last_line == actual_last_line_cf_lf


def test_strip_newlines():
    cr_text = "first\rnext\rlast\r"
    lf_text = "first\nnext\nlast\n"
    cr_lf_text = "first\r\nnext\r\nlast\r\n"
    expected_stripped_text = "firstnextlast"

    actual_stripped_text_cr = strip_newlines(cr_text)
    actual_stripped_text_lf = strip_newlines(lf_text)
    actual_stripped_text_cf_lf = strip_newlines(cr_lf_text)

    assert expected_stripped_text == actual_stripped_text_cr
    assert expected_stripped_text == actual_stripped_text_lf
    assert expected_stripped_text == actual_stripped_text_cf_lf
