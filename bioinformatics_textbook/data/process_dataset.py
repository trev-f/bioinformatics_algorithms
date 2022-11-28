import click
import logging


def main():
    """
    Main logic for the subcommand finalize-dataset.
    This subcommand creates a finalized dataset.
    Imports data from `data/interim` and
    writes finalized dataset to `data/processed`.
    """
    logger = logging.getLogger(__name__)

    logger.info("Start the main function for processing the dataset")

    click.echo("Process the dataset!")

    logger.info("End of the main function for processing the dataset")


if __name__ == "__main__":
    main()
