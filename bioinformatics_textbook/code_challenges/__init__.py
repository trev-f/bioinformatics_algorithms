import click

import logging

import bioinformatics_textbook.code_challenges.ch01
import bioinformatics_textbook.code_challenges.frequent_words


class BA1N:

    def __init__(self, input_file: click.File, logger: logging.Logger = logging.getLogger(__name__)) -> None:
        self.logger = logger
        
        self.logger.info("Create instance to solve BA1N")
        
        pattern_hamming_dist = (
            bioinformatics_textbook.code_challenges.frequent_words.PatternHammingDist(
                input_file
            )
        )

        neighborhood = bioinformatics_textbook.code_challenges.frequent_words.FrequentWords().find_neighbors(
            pattern=pattern_hamming_dist.pattern,
            num_allowed_mismatches=pattern_hamming_dist.hamming_dist,
        )
        formatted_neighborhood = (
            bioinformatics_textbook.code_challenges.inout.RosalindSubmission(
                neighborhood
            ).format_rosalind_answer(sep="\n")
        )

        click.echo(formatted_neighborhood)
