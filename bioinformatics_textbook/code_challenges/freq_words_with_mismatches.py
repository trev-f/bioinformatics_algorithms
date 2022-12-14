import click

import logging

from bioinformatics_textbook.code_challenges.inout import RosalindDataset
from bioinformatics_textbook.code_challenges.ch01 import compute_hamming_distance


class FrequentWords:

    def __init__(self, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.logger = logger

    
    def find_most_freq_words_with_mismatches(self, text: str, kmer_length: int, num_allowed_mismatches: int) -> list:
        self.logger.info("Find most frequent words.")
        # initialize empty collections and 
        patterns = []
        freq_table = {}
        n = len(text)
        for i in range(n - kmer_length + 1):
            pattern = text[i: i+kmer_length]
            neighborhood = self.find_neighbors(pattern, num_allowed_mismatches)
            for neighbor in neighborhood:
                # if a k-mer is not present in frequency table, add it and assign a value of 1,
                # otherwise, increment the count
                freq_table[neighbor] = freq_table.get(neighbor, 0) + 1
        
        max_freq = max(freq_table.values())
        for kmer in freq_table.keys():
            if freq_table[kmer] == max_freq:
                patterns.append(kmer)

        return patterns


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
        suffix_neighbors = self.find_neighbors(self.slice_suffix(pattern), num_allowed_mismatches)
        for suffix_neighbor in suffix_neighbors:
            # the suffix neighbors must have a Hamming distance to pattern that is greater than or equal to the max allowed Hamming distance
            if compute_hamming_distance(self.slice_suffix(pattern), suffix_neighbor) < num_allowed_mismatches:
                for nucleotide in nucleotides:
                    neighborhood.add(nucleotide + suffix_neighbor)
            else:
                neighborhood.add(self.slice_first_nucleotide(pattern) + suffix_neighbor)
        
        return list(neighborhood)


    def slice_suffix(self, pattern: str) -> str:
        """Get the suffix of a string

        :param pattern: The string to slice
        :type pattern: str
        :return: The suffix
        :rtype: str
        """
        return pattern[1:]
    

    def slice_first_nucleotide(self, pattern: str) -> str:
        """Get the first nucleotide of a DNA string

        :param pattern: The DNA string to slice
        :type pattern: str
        :return: The first nucleotide base
        :rtype: str
        """
        return pattern[:1]


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
