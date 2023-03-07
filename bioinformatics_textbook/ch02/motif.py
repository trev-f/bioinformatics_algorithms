import logging

import click

from bioinformatics_textbook.dna import DNA
from bioinformatics_textbook.inout import RosalindDataset


class Motif:
    def find_k_d_motifs(self, kmer_length: int, num_allowed_mismatches: int, dnas: list) -> set:
        """Find (k,d)-motifs in a collection of DNA sequences. That is, find all k-mers that appear in every string of the collection of DNA sequences with at most d mismatches.

        :param kmer_length: Length of (k,d)-motifs
        :type kmer_length: int
        :param num_allowed_mismatches: Number of allowed mismatches (equivalent to maximum Hamming distance)
        :type num_allowed_mismatches: int
        :param dnas: Collection of DNA sequences
        :type dnas: list
        :return: All unique (k,d)-motifs
        :rtype: set
        """
        patterns = set()
        for dna in dnas:
            candidate_patterns = set()
            for i in range(len(dna) - kmer_length + 1):
                kmer = DNA(dna[i: i + kmer_length])
                candidate_patterns.update(kmer.generate_d_neighborhood(num_allowed_mismatches=num_allowed_mismatches))
            if not patterns:
                patterns = candidate_patterns
            else:
                patterns.intersection_update(candidate_patterns)

        return patterns
        

class KDDNA(RosalindDataset):
    """Read and represent a Rosalind dataset that contains k-mer length (k) and maximum allowed mismatches (d) on the first line (separated by a space) and DNA strings on all subsequent lines.
    """

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with k-mer length, maximum allowed mismatches, and DNA strings")

        k, d = self._read_first_line().split(' ')
        self.kmer_length = int(k)
        self.num_allowed_mismatches = int(d)

        self.dnas = self._read_last_lines()

        self._log_init()

    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("k-mer length: %s", self.kmer_length)
        self.logger.info("Max allowed Hamming distance: %s", self.num_allowed_mismatches)
        self.logger.info("Number DNA strings: %s", len(self.dnas))
