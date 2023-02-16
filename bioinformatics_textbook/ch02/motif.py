import logging

import click

from bioinformatics_textbook.inout import RosalindDataset


class Motif:
    pass

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
