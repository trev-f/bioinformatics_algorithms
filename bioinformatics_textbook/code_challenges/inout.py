import click
import re
import os


def read_text_k(input_file: click.File) -> tuple:
    # get the text from all lines except the last
    text = read_not_last_line(input_file)

    # get k from the last line
    k = int(read_last_line(input_file))

    return (text, k)


def read_text_pattern(input_file: click.File) -> tuple:
    # get the text from all lines except the last
    text = read_not_last_line(input_file)

    # get the pattern from the last line
    pattern = read_last_line(input_file)

    return (text, pattern)


def read_all_lines(input_file: click.File) -> str:
    """Read all lines of a file as a single string with no new lines

    :param input_file: Input file
    :type input_file: click.File
    :return: A string with no new lines
    :rtype: str
    """
    all_lines = input_file.read()
    all_lines_stripped = strip_newlines(all_lines)

    return all_lines_stripped


def read_not_last_line(input_file: click.File) -> str:
    """Read every line of a file except for the last line

    :param input_file: The input file. Must be opened for reading in binary mode.
    :type input_file: click.File
    :return: The lines.
    :rtype: str
    """
    input_file.seek(0, os.SEEK_END)
    position = input_file.tell() - 1

    while position > 0 and input_file.read(1) != b"\n":
        position -= 1
        input_file.seek(position, os.SEEK_SET)
    
    if position > 0:
        input_file.seek(0, os.SEEK_SET)
        not_last_lines = input_file.read(position).decode()
    
    not_last_lines_stripped = strip_newlines(not_last_lines)
    
    return not_last_lines_stripped


def read_first_line(input_file: click.File) -> str:
    """Reada the first line of a file

    :param input_file: The input file. Must be opened for reading in binary mode.
    :type input_file: click.File
    :return: The first line of the file.
    :rtype: str
    """
    first_line = input_file.readline().decode()
    first_line_stripped = strip_newlines(first_line)

    return first_line_stripped


def read_second_line(input_file: click.File) -> str:
    """Read the second line of a file

    :param input_file: The input file. Must be opened for reading in bites mode.
    :type input_file: click.File
    :return: The second line.
    :rtype: str
    """
    second_line = input_file.readlines()[1:2][0].decode()
    second_line_stripped = strip_newlines(second_line)

    return second_line_stripped


def read_last_line(input_file: click.File) -> str:
    """Read the last line of a file

    :param input_file: The input file. Must be opened for reading in binary mode.
    :type input_file: click.File
    :return: The last line.
    :rtype: str
    """
    input_file.seek(-2, os.SEEK_END)
    while input_file.read(1) != b"\n":
        input_file.seek(-2, os.SEEK_CUR)
    
    last_line = input_file.readline().decode().rstrip()
    last_line_stripped = strip_newlines(last_line)

    return last_line_stripped


def strip_newlines(text: str) -> str:
    """Strip carriage return and line feed newline characters from a text string

    :param text: A text string
    :type text: str
    :return: The text string stripped of newline characters
    :rtype: str
    """
    stripped_text = re.sub(r"\r|\n", "", text)

    return stripped_text