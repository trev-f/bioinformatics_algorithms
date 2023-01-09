import click
import re
import os

import logging


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
    """Read the first line of a file

    :param input_file: The input file. Must be opened for reading in binary mode.
    :type input_file: click.File
    :return: The first line of the file.
    :rtype: str
    """
    # always set the file object position to the very beginning
    input_file.seek(0, os.SEEK_SET)
    
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
    # always set the file object position to the very beginning
    input_file.seek(0, os.SEEK_SET)

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


class RosalindSubmission():
    
    def __init__(self, answer: list, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.answer = answer if isinstance(answer[0], str) else self.convert_iterable_to_list_of_str(answer)
        self.logger = logger


    def format_list_for_rosalind(self) -> str:
        """Format a list as a string with elements separated by spaces as is commonly expected for solutions to problems for Rosalind.

        :param list_to_format: List to format. If elements are not strings they will be converted.
        :type list_to_format: list
        :return: String of list formatted for Rosalind.
        :rtype: str
        """
        self.logger.info("Formatting answer for submission to Rosalind")
        
        formatted_answer = " ".join(self.answer)
        return formatted_answer
    

    def convert_iterable_to_list_of_str(self, iterable) -> list:
        """Convert an iterable to a list of strings

        :param iterable: An iterable to convert
        :type iterable: _type_
        :return: A list of strings
        :rtype: list
        """
        self.logger.info("Converting iterable answer to a list of strings.")

        try:
            list_of_strings = [str(member) for member in iterable]
        except TypeError as e:
            self.logger.exception("A TypeError error occurred. Checked that the object to convert is an iterable: %s", e)
        else:
            return list_of_strings
