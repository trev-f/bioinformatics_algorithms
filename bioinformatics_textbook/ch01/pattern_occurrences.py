import logging

import click

from bioinformatics_textbook.inout import RosalindDataset


class PatternOccurrences:

    def __init__(self, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.logger = logger
        

    def count_pattern(self, text: str, pattern: str) -> int:
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
    

    def find_starting_positions(self, pattern: str, genome: str) -> list:
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


class TextPattern(RosalindDataset):
    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with a DNA string and a pattern (k-mer).")

        self.text = self._read_not_last_line()
        self.pattern = self._read_last_line()

        self._log_init()

    
    def _log_init(self) -> None:
        """Log attributes created during initialization
        """
        self.logger.info("Text: %s+...", self.text[:10])
        self.logger.info("Pattern: %s", self.pattern)


class PatternGenome(RosalindDataset):
    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        super().__init__(input_file=input_file, logger=logger)

        self.logger.info("Initialize object with a pattern (k-mer) and genome.")

        self.pattern = self._read_not_last_line()
        self.genome = self._read_last_line()

        self._log_init()

    
    def _log_init(self) -> None:
        """Log attributes created during initialization
        """
        self.logger.info("Pattern: %s", self.pattern)
        self.logger.info("Pattern: %s+...", self.genome[:10])
