import logging

import click

from bioinformatics_textbook.inout import RosalindDataset


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

