import click
import os


def read_text_pattern(input_file: click.File) -> tuple:
    other_lines = ""
    # get the last line
    last_line = read_last_line(input_file)

    return (other_lines, last_line)


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
    
    return not_last_lines


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
