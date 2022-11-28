import click
import logging


def main():
    """
    Main logic for the subcommand make-dataset
    This subcommand creates an interim dataset from raw and external data.
    Transform data from `data/raw` and `data/external`
    into an interim dataset in `data/interim`.
    """
    logger = logging.getLogger(__name__)

    logger.info("Start the main function for make dataset")

    click.echo("Make dataset!")

    logger.info("End of the main function for make dataset")


if __name__ == "__main__":
    main()
