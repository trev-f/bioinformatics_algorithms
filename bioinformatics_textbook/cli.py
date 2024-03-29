import logging

import click

import bioinformatics_textbook


class Config(object):
    def __init__(self):
        self.verbose = False
        self.logger = None


pass_config = click.make_pass_decorator(Config, ensure=True)


@click.group()
@click.option("--verbose", "-v", is_flag=True, help="Print more logging messages.")
@pass_config
def cli(config, verbose):
    """Run commands from cfb_rankings_analysis module."""
    config.verbose = verbose
    config.logger = create_root_logger(verbose)

    config.logger.info("Start CLI program.")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1a(config, input_file):
    """
    Program to solve Rosalind problem BA1A: Compute the Number of Times a Pattern Appears in a Text

    https://rosalind.info/problems/ba1a/
    """
    config.logger.info(
        "Run command to solve BA1A: Compute the Number of Times a Pattern Appears in a Text"
    )

    dataset = bioinformatics_textbook.ch01.pattern_occurrences.TextPattern(input_file)
    bioinformatics_textbook.ch01.BA1A(dataset=dataset)

    config.logger.info(
        "Finished command to solve BA1A: Compute the Number of Times a Pattern Appears in a Text"
    )


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1b(config, input_file):
    """
    Program to solve Rosalind problem BA1B: Find the Most Frequent Words in a String

    https://rosalind.info/problems/ba1b/
    """
    config.logger.info(
        "Run command to solve BA1B: Find the most frequent words in a string"
    )

    dataset = bioinformatics_textbook.ch01.frequent_words.TextKmerLength(input_file)
    bioinformatics_textbook.ch01.BA1B(dataset=dataset)

    config.logger.info(
        "Finished command to solve BA1B: Found the most frequent words in a string"
    )


@cli.command()
@click.argument("input_file", type=click.File("r"))
@pass_config
def ba1c(config, input_file):
    """
    Program to solve Rosalind problem BA1C: Find the Reverse Complement of a String

    https://rosalind.info/problems/ba1c/
    """
    config.logger.info("Find the reverse complement of a string")

    dataset = bioinformatics_textbook.ch01.reverse_complement.Pattern(input_file)
    bioinformatics_textbook.ch01.BA1C(dataset=dataset)

    config.logger.info("Found the reverse complement of a string")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1d(config, input_file):
    """
    Program to solve Rosalind problem BA1D: Find All Occurrences of a Pattern in a String

    https://rosalind.info/problems/ba1d/
    """
    config.logger.info("Find all occurrences of a pattern in a string")

    dataset = bioinformatics_textbook.ch01.pattern_occurrences.PatternGenome(input_file)
    bioinformatics_textbook.ch01.BA1D(dataset=dataset)

    config.logger.info("Found all occurrences of a pattern in a string")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1e(config, input_file):
    """
    Program to solve Rosalind problem BA1E: Find Patterns Forming Clumps in a String

    https://rosalind.info/problems/ba1e/
    """
    config.logger.info("Run CLI command to solve BA1E")

    clump_patterns = bioinformatics_textbook.ch01.ch01.ba1e(input_file)
    click.echo(clump_patterns)

    config.logger.info("Finished CLI command to solve BA1E")


@cli.command()
@click.argument("input_file", type=click.File("r"))
@pass_config
def ba1f(config, input_file):
    """Program to solve Rosalind problem BA1F: Find a Position in a Genome Minimizing the Skew.

    https://rosalind.info/problems/ba1f/
    """
    config.logger.info("Run CLI command to solve BA1F")

    min_skew_positions = bioinformatics_textbook.ch01.ch01.ba1f(input_file)
    click.echo(min_skew_positions)

    config.logger.info("Finished CLI command to solve BA1F")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1g(config, input_file):
    """Program to solve Rosalind problem BA1G: Compute the Hamming Distance Between Two Strings

    https://rosalind.info/problems/ba1g/
    """
    config.logger.info("Run CLI command to solve BA1G")

    hamming_distance = bioinformatics_textbook.ch01.ch01.ba1g(input_file)
    click.echo(hamming_distance)

    config.logger.info("Finished CLI command to solve BA1G")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1h(config, input_file):
    """Program to solve Rosalind problem BA1H: Find All Approximate Occurrences of a Pattern in a String

    https://rosalind.info/problems/ba1h/
    """
    config.logger.info("Run CLI command to solve BA1H")

    approx_occurrence_positions = bioinformatics_textbook.ch01.ch01.ba1h(input_file)
    click.echo(approx_occurrence_positions)

    config.logger.info("Finished CLI command to solve BA1H")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1i(config, input_file):
    """Program to solve Rosalind problem BA1I: Find the Most Frequent Words with Mismatches in a String

    https://rosalind.info/problems/ba1i/
    """
    config.logger.info(
        "Run command to solve BA1I: Find the Most Frequent Words with Mismatches in a String"
    )

    dataset = bioinformatics_textbook.ch01.frequent_words.TextKmerLengthHammingDist(
        input_file
    )
    bioinformatics_textbook.ch01.BA1I(dataset=dataset)

    config.logger.info(
        "Finished command to solve BA1I: Find the Most Frequent Words with Mismatches in a String"
    )


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1j(config, input_file):
    """Program to solve Rosalind problem BA1J: Find Frequent Words with Mismatches and Reverse Complements

    https://rosalind.info/problems/ba1j/
    """
    config.logger.info(
        "Run command to solve BA1J: Find Frequent Words with Mismatches and Reverse Complements"
    )

    dataset = bioinformatics_textbook.ch01.frequent_words.TextKmerLengthHammingDist(
        input_file
    )
    bioinformatics_textbook.ch01.BA1J(dataset=dataset)

    config.logger.info(
        "Finished command to solve BA1J: Find Frequent Words with Mismatches and Reverse Complements"
    )


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1n(config, input_file):
    """Program to solve Rosalind problem BA1N: Generate the d-Neighborhood of a String

    https://rosalind.info/problems/ba1n/
    """
    config.logger.info(
        "Run CLI command to solve BA1N: Generate the d-Neighborhood of a String"
    )

    dataset = bioinformatics_textbook.ch01.frequent_words.PatternHammingDist(input_file)
    bioinformatics_textbook.ch01.BA1N(dataset=dataset)

    config.logger.info(
        "Finished CLI command to solve BA1N: Generate the d-Neighborhood of a String"
    )


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba2a(config, input_file):
    """Program to solve Rosalind problem BA2A: Implement MotifEnumeration

    https://rosalind.info/problems/ba2a/
    """
    config.logger.info("Run CLI command to solve BA2A: Implement MotifEnumeration")
    dataset = bioinformatics_textbook.ch02.motif.KDDNA(input_file)
    bioinformatics_textbook.ch02.BA2A(dataset=dataset)


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba2b(config, input_file):
    """Program to solve Rosalind problem BA2B: Find a median string

    https://rosalind.info/problems/ba2h/
    """
    config.logger.info("Run CLI command to solve BA2B: Find a Median String")
    dataset = bioinformatics_textbook.ch02.median_string.KDNAs(input_file)
    bioinformatics_textbook.ch02.BA2B(dataset=dataset)


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba2h(config, input_file):
    """Program to solve Rosalind problem BA2H: Implement DistanceBetweenPatternAndStrings

    https://rosalind.info/problems/ba2h/
    """
    config.logger.info(
        "Run CLI command to solve BA2H: Implement DistanceBetweenPatternAndStrings"
    )
    dataset = bioinformatics_textbook.ch02.median_string.PatternDNAs(input_file)
    bioinformatics_textbook.ch02.BA2H(dataset=dataset)


def create_root_logger(verbose):
    """
    Create a root logger
    """
    # create the logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # create a file handler
    file_handler = logging.FileHandler(".log")
    file_handler.setLevel(logging.DEBUG)

    # create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(set_handler_level(verbose))

    # create formatter and add to handlers
    formatter = logging.Formatter(
        "%(asctime)s :: %(name)s :: %(levelname)s :: %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def set_handler_level(verbose):
    if verbose:
        return logging.DEBUG
    else:
        return logging.ERROR
