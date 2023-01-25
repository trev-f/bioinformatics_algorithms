import click

import logging

from bioinformatics_textbook.ch01.inout import RosalindDataset
from bioinformatics_textbook.ch01.ch01 import (
    compute_hamming_distance, reverse_complement_dna
)


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
            kmer = text[i: i+kmer_length]
            neighborhood = self.find_neighbors(kmer, num_allowed_mismatches)
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
            kmer = text[i: i+kmer_length]
            neighborhood = self.find_neighbors(kmer, num_allowed_mismatches)
            for neighbor in neighborhood:
                freq_table[neighbor] = freq_table.get(neighbor, 0) + 1
                rc = reverse_complement_dna(dna=neighbor)
                freq_table[rc] = freq_table.get(rc, 0) + 1

        # compute most frequent k-mers with reverse compliments
        most_freq_words = []
        max_freq = self._find_max_val_of_dict(d=freq_table)
        for kmer in freq_table.keys():
            if freq_table[kmer] == max_freq:
                most_freq_words.append(kmer)

        return most_freq_words


    def find_neighbors(self, pattern: str, num_allowed_mismatches: int) -> list:
        """Find the set of all k-mers whose Hamming distance from a pattern does not exceed a set maximum.

        :param pattern: A DNA string k-mer
        :type pattern: str
        :param num_allowed_mismatches: The maximum allowed Hamming distance.
        :type num_allowed_mismatches: int
        :return: k-mers that are in the specified neighborhood of pattern.
        :rtype: list
        """
        nucleotides = {"A", "C", "G", "T"}
        
        # if there are no allowed mismatches the pattern's neighborhood is itself
        if num_allowed_mismatches == 0:
            return [pattern]

        # if the pattern is only a single nucleotide, its neighborhood is the set of all nucleotides
        # this  serves as the base condition for the recursion
        if len(pattern) == 1:
            return list(nucleotides)
        
        # making a neighborhood a set instead of a list means duplicates do not have to be explicitly removed
        neighborhood = set()
        suffix_neighbors = self.find_neighbors(self._slice_suffix(pattern), num_allowed_mismatches)
        for suffix_neighbor in suffix_neighbors:
            # the suffix neighbors must have a Hamming distance to pattern that is greater than or equal to the max allowed Hamming distance
            if compute_hamming_distance(self._slice_suffix(pattern), suffix_neighbor) < num_allowed_mismatches:
                for nucleotide in nucleotides:
                    neighborhood.add(nucleotide + suffix_neighbor)
            else:
                neighborhood.add(self._slice_first_nucleotide(pattern) + suffix_neighbor)
        
        return list(neighborhood)
    

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
    

    def _construct_kmer_rc_neighborhood(self, kmer_neighborhood: list) -> list:
        """Construct a k-mer mismatch neighborhood that also includes reverse complements

        :param kmer_neighborhood: The d-neighborhood of a k-mer
        :type kmer_neighborhood: list
        :return: The d-neighborhood of a k-mer with its reverse complements. Duplicates removed.
        :rtype: list
        """
        rc_neighborhood = [reverse_complement_dna(kmer) for kmer in kmer_neighborhood]

        # coerce to set to remove duplicates then back to list
        return list(set(kmer_neighborhood + rc_neighborhood))



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
