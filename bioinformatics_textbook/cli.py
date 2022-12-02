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
@pass_config
def download_data(config):
    """
    Download external data.
    Download external data into `data/external`.
    """
    config.logger.info("Download external data")

    bioinformatics_textbook.data.download_data.main()

    config.logger.info("External data download complete")


@cli.command()
@pass_config
def make_dataset(config):
    """
    Make an interim dataset.
    Transform data from `data/raw` and/or `data/external` into a dataset in `data/interim`.
    """
    config.logger.info("Make dataset")

    bioinformatics_textbook.data.make_dataset.main()

    config.logger.info("Dataset made")


@cli.command()
@pass_config
def process_dataset(config):
    """
    Make a final processed dataset.
    Transform data from `data/interim` into a final processed dataset in `data/processed`.
    """
    config.logger.info("Make final processed dataset")

    bioinformatics_textbook.data.process_dataset.main()

    config.logger.info("Final processed dataset made")


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
