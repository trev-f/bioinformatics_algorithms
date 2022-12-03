import bioinformatics_textbook.code_challenges.inout
import click


def ba1b(input_file: click.File) -> str:
    text, k = bioinformatics_textbook.code_challenges.inout.read_text_k(input_file)

    most_frequent_words = find_frequent_words(text, k)

    return most_frequent_words


def find_frequent_words(text: str, k: int) -> str:

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
