from click.testing import CliRunner
from bioinformatics_textbook.cli import cli


def test_download_data():
    runner = CliRunner()
    result = runner.invoke(cli, ["download-data"])
    assert result.exit_code == 0
    assert "Download external data!" in result.output


def test_make_dataset():
    runner = CliRunner()
    result = runner.invoke(cli, ["make-dataset"])
    assert result.exit_code == 0
    assert "Make dataset!" in result.output


def test_process_dataset():
    runner = CliRunner()
    result = runner.invoke(cli, ["process-dataset"])
    assert result.exit_code == 0
    assert "Process the dataset!" in result.output
