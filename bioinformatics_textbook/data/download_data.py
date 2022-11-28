import click
import logging


def main():
    """
    Main logic for the subcommand download-data.
    This subcommand downloads data from an external source to `data/external`.
    """
    logger = logging.getLogger(__name__)

    logger.info("Start the main function for downloading external data")

    click.echo("Download external data!")

    logger.info("End of the main function for downloading external data")


if __name__ == "__main__":
    main()
