import bioinformatics_textbook.code_challenges.inout
import click


def ba1d(input_file: click.File) -> str:
    pattern = bioinformatics_textbook.code_challenges.inout.read_not_last_line(input_file)
    genome = bioinformatics_textbook.code_challenges.inout.read_last_line(input_file)

    starting_positions = find_starting_positions(pattern, genome)
    formatted_starting_positions = format_starting_positions(starting_positions)

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


def format_starting_positions(starting_positions: list) -> str:
    """Format a list of starting positions to a string separated by spaces.

    :param starting_positions: Starting positions
    :type starting_positions: list
    :return: Formatted starting positions string
    :rtype: str
    """
    starting_positions = [str(position) for position in starting_positions]
    formatted_starting_positions = " ".join(starting_positions)

    return formatted_starting_positions


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

    return " ".join(most_freq_words)


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
