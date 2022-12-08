import bioinformatics_textbook.code_challenges.inout
import click
from collections import OrderedDict
import logging


def find_clumps(genome: str, pattern_length: int, window_length: int, pattern_freq_thresh: int) -> str:
    """Find k-mers that are found in clumps in the genome

    :param genome: A DNA string to search for clumps
    :type genome: str
    :param pattern_length: Length of k-mers
    :type pattern_length: int
    :param window_length: Sliding window size
    :type window_length: int
    :param pattern_freq_thresh: The minimum number of times a k-mer must appear within a window for the k-mer to form a clump
    :type pattern_freq_thresh: int
    :return: All unique k-mers that form clumps, in order of first appearance, separated by spaces
    :rtype: str
    """
    clump_patterns = []
    genome_length = len(genome)
    for i in range(genome_length - window_length + 1):
        window = genome[i: i+window_length]
        freq_table = construct_kmer_freq_table(window, pattern_length)
        for key in freq_table.keys():
            if freq_table[key] >= pattern_freq_thresh:
                clump_patterns.append(key)

    # remove duplicates from list
    # note: this is not the fastest way to get unique elements of a list,
    # but it does preserve order which helps in checking against the Rosalind sample output
    # Rosalind appears to report unique patterns in the order in which they appear first in the string
    unique_clump_patterns = list(OrderedDict.fromkeys(clump_patterns))

    # format clumps: separate each kmer by a space
    formatted_clump_patterns = format_list_for_rosalind(unique_clump_patterns)

    return formatted_clump_patterns


def ba1d(input_file: click.File) -> str:
    pattern = bioinformatics_textbook.code_challenges.inout.read_not_last_line(input_file)
    genome = bioinformatics_textbook.code_challenges.inout.read_last_line(input_file)

    starting_positions = find_starting_positions(pattern, genome)
    formatted_starting_positions = format_list_for_rosalind(starting_positions)

    return formatted_starting_positions


def find_starting_positions(pattern: str, genome: str) -> list:
    """Find all occurrences of a pattern (k-mer) in a string (genome)

    :param pattern: A k-mer sequence
    :type pattern: str
    :param genome: A DNA string (genome)
    :type genome: str
    :return: Starting positions of each occurrence of pattern in genome
    :rtype: list
    """
    k = len(pattern)
    genome_length = len(genome)

    starting_positions = []
    for i in range(genome_length - k + 1):
        window = genome[i: i + k]
        if window == pattern:
            starting_positions.append(i)

    return starting_positions


def ba1c(input_file: click.File) -> str:
    """Find the reverse complement of a string of DNA

    :param input_file: A file containing a string of DNA
    :type input_file: click.File
    :return: The reverse complement of the DNA string
    :rtype: str
    """
    dna_string = bioinformatics_textbook.code_challenges.inout.read_all_lines(input_file)

    reverse_complement = reverse_complement_dna(dna_string)

    return reverse_complement


def reverse_complement_dna(dna: str) -> str:
    """Find the reverse complement of a DNA string

    :param dna: DNA string
    :type dna: str
    :return: Reverse complement of DNA string
    :rtype: str
    """
    complement = complement_dna(dna)
    reverse_complement = complement[::-1]
    
    return reverse_complement


def complement_dna(dna: str) -> str:
    """Find the complement of a DNA string

    :param dna: DNA string
    :type dna: str
    :return: Complement of DNA string
    :rtype: str
    """
    comp_table = make_dna_complementation_table()
    complement = dna.translate(comp_table)

    return complement


def make_dna_complementation_table() -> dict:
    """Construct a translation table for complementing DNA

    :return: DNA complementation table
    :rtype: dict
    """
    comp_table = str.maketrans("ATCG", "TAGC")

    return comp_table


def ba1b(input_file: click.File) -> str:
    text, k = bioinformatics_textbook.code_challenges.inout.read_text_k(input_file)

    most_frequent_words = find_frequent_words(text, k)

    return most_frequent_words


def find_frequent_words(text: str, k: int) -> str:
    """Return the most frequent kmers from a string of text

    :param text: A string of text (typically a DNA string)
    :type text: str
    :param k: k-mer length
    :type k: int
    :return: The most frequent k-mers each separated by a space
    :rtype: str
    """
    freq_table = construct_kmer_freq_table(text, k)
    max_freq = find_max_val_of_dict(freq_table)
    
    most_freq_words = []
    for pattern in freq_table.keys():
        if freq_table[pattern] == max_freq:
            most_freq_words.append(pattern)
    
    formatted_most_freq_words = format_list_for_rosalind(most_freq_words)

    return formatted_most_freq_words


def find_max_val_of_dict(d: dict) -> float:
    """Find the max value of a dictionary

    :param d: A dictionary
    :type d: dict
    :return: The dictionary's maximum value
    :rtype: float
    """
    max_val = max(d.values())

    return max_val


def construct_kmer_freq_table(text: str, k: int) -> dict:
    """Construct a frequency table of how many times all k-mers appear in a text

    :param text: A string of text (typically a DNA string)
    :type text: str
    :param k: k-mer length
    :type k: int
    :return: Frequency table of k-mers and their counts
    :rtype: dict
    """
    freq_table = {}
    text_length = len(text)

    # slide windows of length k down the text string
    for i in range(text_length - k + 1):
        pattern = text[i: i + k]
        # if a k-mer is not present in frequency table, add it and assign a value of 1,
        # otherwise, increment the count
        freq_table[pattern] = freq_table.get(pattern, 0) + 1

    return freq_table


def ba1a(input_file: click.File) -> int:
    text, pattern = bioinformatics_textbook.code_challenges.inout.read_text_pattern(input_file)
    
    kmer_count = count_pattern(text, pattern)

    return kmer_count


def count_pattern(text: str, pattern: str) -> int:
    """Count the number of times a k-mer pattern appears as a substring of text using the sliding window method

    :param text: A string of text (typically a DNA string)
    :type text: str
    :param pattern: A k-mer of length less than or equal to that of `text`
    :type pattern: str
    :return: A count of times the k-mer pattern appears as a substring of `text`
    :rtype: int
    """
    number_kmer_appearances = 0
    
    kmer_length = len(pattern)
    number_windows = (len(text) - kmer_length) + 1
    for i in range(number_windows):
        if text[i:i + kmer_length] == pattern:
            number_kmer_appearances += 1

    return number_kmer_appearances


def format_list_for_rosalind(list_to_format: list) -> str:
    """Format a list as a string with elements separated by spaces as is commonly expected for solutions to problems for Rosalind.

    :param list_to_format: List to format. If elements are not strings they will be converted.
    :type list_to_format: list
    :return: String of list formatted for Rosalind.
    :rtype: str
    """
    if not isinstance(list_to_format[0], str):
        list_to_format = convert_iterable_to_list_of_str(list_to_format)
    formatted_list = " ".join(list_to_format)

    return formatted_list


def convert_iterable_to_list_of_str(iterable) -> list:
    """Convert an iterable to a list of strings

    :param iterable: An iterable to convert
    :type iterable: _type_
    :return: A list of strings
    :rtype: list
    """
    try:
        list_of_strings = [str(member) for member in iterable]
    except TypeError as e:
        logging.exception("A TypeError error occurred. Checked that the object to convert is an iterable: %s", e)
    else:
        return list_of_strings
