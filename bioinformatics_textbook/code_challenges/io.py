import click
import os


def read_text_pattern(input_file: click.File) -> tuple:
    other_lines = ""
    # get the last line
    last_line = read_last_line(input_file)

    return (other_lines, last_line)


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

    return last_line
