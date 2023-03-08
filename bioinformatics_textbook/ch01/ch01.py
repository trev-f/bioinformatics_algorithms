import logging
from collections import OrderedDict

import click

import bioinformatics_textbook.inout
from bioinformatics_textbook.dna import DNA


logger = logging.getLogger(__name__)

def ba1h(input_file: click.File):
    pattern = DNA(bioinformatics_textbook.inout.read_first_line(input_file))
    text = bioinformatics_textbook.inout.read_second_line(input_file)
    num_allowed_mismatches = int(bioinformatics_textbook.inout.read_last_line(input_file))

    logger.info("Pattern = %s", pattern)
    logger.info("Text = %s", text)
    logger.info("Number allowed mismatches = %s", num_allowed_mismatches)

    approx_occurrence_positions = find_approx_occurrence_positions(
        pattern, text, num_allowed_mismatches
    )

    formatted_approx_occurrence_positions = format_list_for_rosalind(approx_occurrence_positions)

    return formatted_approx_occurrence_positions


def find_approx_occurrence_positions(pattern: DNA, text: str, num_allowed_mismatches: int) -> list:
    """
    APPROACH:
        Slide k-mer along text. At each position, compute Hamming distance for k-mer and the window of text.

        Idea to improve efficiency: As soon as soon as number of allowed mismatches is reach, move to next window.

    PSEUDOCODE:
        ApproximateOccurrences <- []
        for i <- 0 to |Text| - |Pattern|
            if HammingDist(Pattern, Text[i, |Pattern|]) <= NumberAllowedMismatches
                ApproximateOccurrences.append(i)
        return ApproximateOccurrences
    """
    approx_occurrence_positions = []
    text_length = len(text)
    kmer_length = len(pattern)
    for i in range(text_length - kmer_length + 1):
        if pattern.compute_hamming_distance(text[i: i+kmer_length]) <= num_allowed_mismatches:
            approx_occurrence_positions.append(i)

    return approx_occurrence_positions


def ba1g(input_file: click.File) -> int:
    dna_p = DNA(bioinformatics_textbook.inout.read_not_last_line(input_file))
    dna_q = bioinformatics_textbook.inout.read_last_line(input_file)

    hamming_distance = dna_p.compute_hamming_distance(dna_q)

    return hamming_distance


def ba1f(input_file: click.File) -> str:
    """Find a position in a genome minimizing the GC skew

    :param input_file: A text file that defines a genome string
    :type input_file: click.File
    :return: The positions in a genome that minimize GC skew
    :rtype: str
    """
    genome = bioinformatics_textbook.inout.read_all_lines(input_file)

    gc_skews = define_dna_gc_skews(genome)
    min_skew_positions = find_min_skew_positions(gc_skews)

    return format_list_for_rosalind(min_skew_positions)


def find_min_skew_positions(skews: list) -> list:
    """Find positions of a genome where skew is minimal

    :param skews: Skew at each position of a genome (starting from 0)
    :type skews: list
    :return: Genome positions where skew is minimal
    :rtype: list
    """
    min_skew = min(skews)
    min_skew_positions = [i for i, skew in enumerate(skews) if skew == min_skew]

    return min_skew_positions


def define_dna_gc_skews(genome: str) -> list:
    """Define the GC skew of a DNA string as the as the difference between the total number of occurrences of 'G' and 'C' in a genome.

    :param genome: A DNA string
    :type genome: str
    :return: The GC skew at each position of a genome. The list starts from 0 and therefore is out of phase with the genome string.
    :rtype: list
    """
    skews = [0]
    for i, base in enumerate(genome):
        if base == "C":
            skews.append(skews[i] - 1)
        elif base == "G":
            skews.append(skews[i] + 1)
        else:
            skews.append(skews[i])

    return skews


def ba1e(input_file: click.File) -> str:
    """Find patterns forming clumps in a string

    :param input_file: An input file where the first line is a genome string and the second line integers denoting k-mer length, window size, and pattern frequency threshold.
    :type input_file: click.File
    :return: A string-separated list of k-mers that form clumps
    :rtype: str
    """
    genome = bioinformatics_textbook.inout.read_not_last_line(input_file)
    last_line = bioinformatics_textbook.inout.read_last_line(input_file)
    k, L, t = [int(element) for element in last_line.split(" ")]

    clump_patterns = find_clumps(genome, k, L, t)

    return clump_patterns


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
