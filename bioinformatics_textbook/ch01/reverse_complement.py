import click

import logging

from bioinformatics_textbook.inout import RosalindDataset


class ReverseComplement:

    def __init__(self, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.logger = logger

    
    def reverse_complement_dna(self, dna: str) -> str:
        """Find the reverse complement of a DNA string

        :param dna: DNA string
        :type dna: str
        :return: Reverse complement of DNA string
        :rtype: str
        """
        complement = self.complement_dna(dna)
        reverse_complement = self.reverse_text(complement)
        
        return reverse_complement


    def complement_dna(self, dna: str) -> str:
        """Find the complement of a DNA string

        :param dna: DNA string
        :type dna: str
        :return: Complement of DNA string
        :rtype: str
        """
        comp_table = self._make_dna_complementation_table()
        complement = dna.translate(comp_table)

        return complement


    def _make_dna_complementation_table(self) -> dict:
        """Construct a translation table for complementing DNA

        :return: DNA complementation table
        :rtype: dict
        """
        comp_table = str.maketrans("ATCG", "TAGC")

        return comp_table
    

    def reverse_text(self, text: str) -> str:
        """Reverse a text string.

        :param text: A string of text (typically a DNA string).
        :type text: str
        :return: The reverse of the text string.
        :rtype: str
        """
        return text[::-1]


class Pattern(RosalindDataset):
    
    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with pattern")

        self.pattern = self._read_all_lines()

        self._log_init()


    def _log_init(self) -> None:
        """Log attributes created during initializtion
        """
        self.logger.info("Pattern: %s", self.pattern)

