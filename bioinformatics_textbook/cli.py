import click
import logging
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
    config.logger.info("Compute the number of times a pattern appears in a text")

    kmer_count = bioinformatics_textbook.code_challenges.ch01.ba1a(input_file)
    click.echo(kmer_count)

    config.logger.info("Computed the number of times a pattern appears in a text")


@cli.command()
@click.argument("input_file", type=click.File("rb"))
@pass_config
def ba1b(config, input_file):
    """
    Program to solve Rosalind problem BA1B: Find the Most Frequent Words in a String

    https://rosalind.info/problems/ba1b/
    """
    config.logger.info("Find the most frequent words in a string")

    most_frequent_words = bioinformatics_textbook.code_challenges.ch01.ba1b(input_file)
    click.echo(most_frequent_words)

    config.logger.info("Found the most frequent words in a string")


@cli.command()
@click.argument("input_file", type=click.File("r"))
@pass_config
def ba1c(config, input_file):
    """
    Program to solve Rosalind problem BA1C: Find the Reverse Complement of a String

    https://rosalind.info/problems/ba1c/
    """
    config.logger.info("Find the reverse complement of a string")

    reverse_complement = bioinformatics_textbook.code_challenges.ch01.ba1c(input_file)
    click.echo(reverse_complement)

    config.logger.info("Found the reverse complement of a string")


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
