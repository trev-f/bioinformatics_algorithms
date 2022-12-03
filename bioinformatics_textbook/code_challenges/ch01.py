import bioinformatics_textbook.code_challenges.inout
import click


def ba1b(input_file: click.File) -> str:
    text, k = bioinformatics_textbook.code_challenges.inout.read_text_k(input_file)

    most_frequent_words = find_frequent_words(text, k)

    return most_frequent_words


def find_frequent_words(text: str, k: int) -> str:

    """
    BetterFrequentWords(Text, k)
    FrequentPatterns ← an array of strings of length 0
    freqMap ← FrequencyTable(Text, k)
    max ← MaxMap(freqMap)
    for all strings Pattern in freqMap
        if freqMap[pattern] = max
            append Pattern to frequentPatterns
    return frequentPatterns
    """
    return ""


def construct_kmer_freq_table(text: str, k: int) -> dict:

    """
    FrequencyTable(Text, k)
    freqMap ← empty map
    n ← |Text|
    for i ← 0 to n − k
        Pattern ← Text(i, k)
        if freqMap[Pattern] doesn't exist
            freqMap[Pattern]← 1
        else
           freqMap[pattern] ←freqMap[pattern]+1 
    return freqMap
    """
    freq_table = {}
    text_length = len(text)

    for i in range(text_length - k + 1):
        pattern = text[i: i + k]
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
