import logging

import click

from bioinformatics_textbook.inout import RosalindDataset
from bioinformatics_textbook.dna import DNA


class FrequentWords:

    def __init__(self, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.logger = logger


    def find_most_freq_words(self, text: str, kmer_length: int) -> list:
        """Find the most frequent k-mers in a string of text

        :param text: A string of text (typically a DNA string)
        :type text: str
        :param kmer_length: k-mer length
        :type kmer_length: int
        :return: The most frequent k-mers in the text
        :rtype: list
        """
        freq_table = self._construct_kmer_freq_table(text=text, kmer_length=kmer_length)
        max_freq = self._find_max_val_of_dict(d=freq_table)

        most_freq_words = []
        for pattern in freq_table.keys():
            if freq_table[pattern] == max_freq:
                most_freq_words.append(pattern)
        
        return most_freq_words

    
    def find_most_freq_words_with_mismatches(self, text: str, kmer_length: int, num_allowed_mismatches: int) -> list:
        """Find the most frequent k-mers with up to a number of allowed mismatches in a string of text

        :param text: A string of text (typically a DNA string)
        :type text: str
        :param kmer_length: k-mer length
        :type kmer_length: int
        :param num_allowed_mismatches: The maximum Hamming distance (the number of allowed mismatches)
        :type num_allowed_mismatches: int
        :return: The most frequent k-mers in the text with at most the allowed number of mismatches
        :rtype: list
        """
        self.logger.info("Find most frequent words with mismatches.")

        # construct frequency table of k-mers with mismatches
        freq_table = {}
        for i in range(self._compute_number_sliding_windows(text=text, kmer_length=kmer_length)):
            kmer = DNA(text[i: i+kmer_length])
            neighborhood = kmer.generate_d_neighborhood(num_allowed_mismatches=num_allowed_mismatches)
            for neighbor in neighborhood:
                # if a k-mer is not present in frequency table, add it and assign a value of 1,
                # otherwise, increment the count
                freq_table[neighbor] = freq_table.get(neighbor, 0) + 1
        
        # select most frequent words
        most_freq_words = []
        max_freq = self._find_max_val_of_dict(d=freq_table)
        for kmer in freq_table.keys():
            if freq_table[kmer] == max_freq:
                most_freq_words.append(kmer)

        return most_freq_words
    

    def find_most_freq_words_with_mismatches_and_rc(self, text: str, kmer_length: int, num_allowed_mismatches: int) -> list:
        self.logger.info("Find most frequent words with mismatches and reverse complements.")

        # construct frequency table of k-mers with mismatches and reverse complements
        freq_table = {}
        for i in range(self._compute_number_sliding_windows(text=text, kmer_length=kmer_length)):
            kmer = DNA(text[i: i+kmer_length])
            neighborhood = kmer.generate_d_neighborhood(num_allowed_mismatches=num_allowed_mismatches)
            for neighbor in neighborhood:
                freq_table[neighbor] = freq_table.get(neighbor, 0) + 1
                rc = DNA(neighbor).reverse_complement()
                freq_table[rc] = freq_table.get(rc, 0) + 1

        # compute most frequent k-mers with reverse compliments
        most_freq_words = []
        max_freq = self._find_max_val_of_dict(d=freq_table)
        for kmer in freq_table.keys():
            if freq_table[kmer] == max_freq:
                most_freq_words.append(kmer)

        return most_freq_words

    def _construct_kmer_freq_table(self, text: str, kmer_length: int) -> dict:
        """Construct a frequency table of how many times all k-mers appear in a text

        :param text: A string of text (typically a DNA string)
        :type text: str
        :param k: k-mer length
        :type k: int
        :return: Frequency table of k-mers and their counts
        :rtype: dict
        """
        freq_table = {}
        # slide windows of length k down the text string
        for i in range(self._compute_number_sliding_windows(text, kmer_length)):
            pattern = text[i: i + kmer_length]
            # if a k-mer is not present in frequency table, add it and assign a value of 1,
            # otherwise, increment the count
            freq_table[pattern] = freq_table.get(pattern, 0) + 1

        return freq_table


    def _slice_suffix(self, pattern: str) -> str:
        """Get the suffix of a string

        :param pattern: The string to slice
        :type pattern: str
        :return: The suffix
        :rtype: str
        """
        return pattern[1:]
    

    def _slice_first_nucleotide(self, pattern: str) -> str:
        """Get the first nucleotide of a DNA string

        :param pattern: The DNA string to slice
        :type pattern: str
        :return: The first nucleotide base
        :rtype: str
        """
        return pattern[:1]


    def _find_max_val_of_dict(self, d: dict) -> float:
        """Find the max value of a dictionary

        :param d: A dictionary
        :type d: dict
        :return: The dictionary's maximum value
        :rtype: float
        """
        max_val = max(d.values())

        return max_val
    

    def _compute_number_sliding_windows(self, text: str, kmer_length: int) -> int:
        """Convenience method to compute the number of sliding windows down a string of text. The return value will most commonly be the input of a range constructor.

        :param text: A string of text (typically a DNA string).
        :type text: str
        :param kmer_length: k-mer length
        :type kmer_length: int
        :return: The number of sliding windows to use along text.
        :rtype: int
        """
        return len(text) - kmer_length + 1


class PatternHammingDist(RosalindDataset):

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with pattern and max allowed Hamming distance")

        self.pattern = self._read_first_line()
        self.hamming_dist = int(self._read_last_line())

        self._log_init()


    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("Pattern: %s", self.pattern)
        self.logger.info("Max allowed Hamming distance: %s", self.hamming_dist)


class TextKmerLengthHammingDist(RosalindDataset):

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with a DNA string, k-mer length, and max allowed Hamming distance.")

        self.text = self._read_first_line()
        k_hamming = self._read_last_line()
        self.kmer_length = int(k_hamming.split(" ")[0])
        self.hamming_dist = int(k_hamming.split(" ")[1])

        self._log_init()


    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("Text: %s+...", self.text[:10])
        self.logger.info("k-mer length: %s", self.kmer_length)
        self.logger.info("Max allowed Hamming distance: %s", self.hamming_dist)


class TextKmerLength(RosalindDataset):
    
    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)):
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with a DNA string and k-mer length.")

        self.text = self._read_not_last_line()
        self.kmer_length = int(self._read_last_line())

        self._log_init()
    

    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("Text: %s+...", self.text[:10])
        self.logger.info("k-mer length: %s", self.kmer_length)
