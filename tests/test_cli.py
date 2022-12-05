from click.testing import CliRunner
from bioinformatics_textbook.cli import cli


def test_ba1a():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1a", "tests/datasets/ch01/ba1a_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "2"


def test_ba1b():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1b", "tests/datasets/ch01/ba1b_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "GCAT CATG"
