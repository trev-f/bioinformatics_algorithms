from click.testing import CliRunner
from bioinformatics_textbook.cli import cli


def test_ba1f():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1f", "tests/datasets/ch01/ba1f_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "53 97"


def test_ba1e():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1e", "tests/datasets/ch01/ba1e_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "CGACA GAAGA AATGT"


def test_ba1d():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1d", "tests/datasets/ch01/ba1d_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "1 3 9"
    

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


def test_ba1c():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1c", "tests/datasets/ch01/ba1c_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "ACCGGGTTTT"
