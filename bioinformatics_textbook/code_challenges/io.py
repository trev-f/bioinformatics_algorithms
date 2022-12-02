import click
import os


def read_text_pattern(input_file: click.File) -> tuple:
    other_lines = ""
    # get the last line
    last_line = read_last_line(input_file)

    return (other_lines, last_line)


def read_last_line(input_file: click.File) -> str:
    input_file.seek(-2, os.SEEK_END)
    while input_file.read(1) != b"\n":
        input_file.seek(-2, os.SEEK_CUR)
    
    last_line = input_file.readline().decode().rstrip()

    return last_line
