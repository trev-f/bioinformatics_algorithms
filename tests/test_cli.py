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


def test_ba1c():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1c", "tests/datasets/ch01/ba1c_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "ACCGGGTTTT"


def test_ba1d():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1d", "tests/datasets/ch01/ba1d_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "1 3 9"


def test_ba1e():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1e", "tests/datasets/ch01/ba1e_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "CGACA GAAGA AATGT"


def test_ba1f():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1f", "tests/datasets/ch01/ba1f_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "53 97"


def test_ba1g():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1g", "tests/datasets/ch01/ba1g_sample_dataset.txt"])
    assert result.exit_code == 0
    assert result.output.rstrip() == "3"


def test_ba1h():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1h", "tests/datasets/ch01/ba1h_sample_dataset.txt"])
    
    assert result.exit_code == 0

    expected_freq_words = set("6 7 26 27 78".split(" "))
    actual_freq_words = set(result.output.rstrip().split(" "))

    assert expected_freq_words == actual_freq_words


def test_ba1i():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1i", "tests/datasets/ch01/ba1i_sample_dataset.txt"])
    
    assert result.exit_code == 0

    expected_freq_words = set("GATG ATGC ATGT".split(" "))
    actual_freq_words = set(result.output.rstrip().split(" "))

    assert expected_freq_words == actual_freq_words
    

def test_ba1j():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1j", "tests/datasets/ch01/ba1j_sample_dataset.txt"])
    
    assert result.exit_code == 0

    expected_freq_words = set("ATGT ACAT".split(" "))
    actual_freq_words = set(result.output.rstrip().split(" "))

    assert actual_freq_words == expected_freq_words


def test_ba1n():
    runner = CliRunner()
    result = runner.invoke(cli, ["ba1n", "tests/datasets/ch01/ba1n_sample_dataset.txt"])
    
    assert result.exit_code == 0

    expected_neighborhood = {
        "CCG", "TCG", "GCG", "AAG", "ATG", "AGG", "ACA", "ACC", "ACT", "ACG"
    }
    actual_neighborhood = set(result.output.rstrip().split("\n"))

    assert expected_neighborhood == actual_neighborhood
