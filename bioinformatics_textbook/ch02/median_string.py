import logging

import click

from bioinformatics_textbook.dna import DNA
from bioinformatics_textbook.inout import RosalindDataset


class MedianString:

    def __init__(self) -> None:
        self.median_strings = []
        self.minimum_distance = None
    
    def find_median_strings(self, kmer_length: int, dnas: list[DNA]) -> list[DNA]:
        """Find median string(s), i.e. k-mer(s) that minimize the distance between all k-mers and DNA sequences.

        :param kmer_length: k-mer length
        :type kmer_length: int
        :param dnas: DNA sequences
        :type dnas: list[DNA]
        :return: Median string(s)
        :rtype: list[DNA]
        """
        self.minimum_distance = kmer_length * len(dnas)
        for kmer in DNA.generate_all_possible_kmers(kmer_length=kmer_length):
            distance = self.compute_pattern_strings_distance(pattern=kmer, dnas=dnas)
            if distance == self.minimum_distance:
                self.median_strings.append(kmer)
            elif distance < self.minimum_distance:
                self.median_strings = [kmer]
                self.minimum_distance = distance

        return self.median_strings
        

    def compute_pattern_strings_distance(self, pattern: DNA, dnas: list[DNA]) -> int:
        """Compute the distance between a pattern and a collection of DNA sequences

        :param pattern: k-mer
        :type pattern: DNA
        :param dnas: Collection of DNA sequences
        :type dnas: list[DNA]
        :return: Distance between pattern and DNA sequences
        :rtype: int
        """
        global_distance = 0
        for dna in dnas:
            local_distance = len(pattern)
            for kmer in dna.generate_kmers(kmer_length=len(pattern)):
                hamming_distance = pattern.compute_hamming_distance(dna_q=kmer)
                local_distance = hamming_distance if hamming_distance < local_distance else local_distance
            
            global_distance += local_distance

        return global_distance


class KDNAs(RosalindDataset):
    """Read and represent a Rosalind dataset that contains a k-mer 'pattern' and DNA strings
    """

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with k-mer pattern and DNA strings")

        self.k = int(self._read_first_line())
        self.dnas = [DNA(seq) for seq in self._read_last_lines()]

        self._log_init()

    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("k: %s", self.k)
        self.logger.info("Number DNA strings: %s", len(self.dnas))


class PatternDNAs(RosalindDataset):
    """Read and represent a Rosalind dataset that contains a k-mer 'pattern' and DNA strings
    """

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with k-mer pattern and DNA strings")

        self.pattern = DNA(self._read_first_line())
        self.dnas = [DNA(seq) for seq in self._read_last_line().split(' ')]

        self._log_init()

    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("Pattern: %s", self.pattern)
        self.logger.info("Number DNA strings: %s", len(self.dnas))
