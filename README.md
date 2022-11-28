# bioinformatics_textbook

Notes and solutions for Bioinformatics Algorithms 3rd ed. by Compeau and Pevzner.

## Usage

### Initial setup

Create and activate a virtual environment.

```bash
make virtual_environment
. venv/bin/activate
```

Run an initial check to be sure the template was installed correctly, install the included Python package as an executable (allows for use of CLI), install required external packages,
and make sure intial linting, tests, and formatting runs correctly. All these checks and tests should pass.

```bash
make check
```

### Development

This template was designed with test driven development (TDD) and continuous integration (CI) in mind.

For TDD, add tests to the `tests/` directory.
Throughtout the development cycle, after changes are made, continuously run `make check` to perform:

  1. Installation of packages in `requirements.txt` by `pip`.
  2. Linting of `.py` files with `pylint`.
  3. Testing of `tests/` with `pytest`.
  4. Formatting of `.py` files with `Black`.

For CI, use entrypoints in `Makefile`, such as `make check` to run automated tests with your platform of choice.
I'm a big fan of simply using `make check` in GitHub Actions.

### Running scripts

This template was designed to implement an easy to use CLI implemented by `Click`.
The CLI code is implemented in `bioinformatics_textbook/cli.py`.
All scripts are available as subcommands under a head command specified by the project name.

To list the available subcommands and options:

```bash
bioinformatics-textbook --help
```

For convenience, there are already scripts and subcommands templated for three common data processing operations:

  1. `download-data`: Download external data
  2. `make-dataset`: Transform raw or external data into an interim dataset
  3. `process-dataset`: Transform an interim dataset into a final processed dataset

For added convenience, a `make` command is provided that will perform all of these data commands.
This aids reproducibility as collaborators (or your future self) can simply add the raw data and/or API credentials
and run the following to produce the exact same data:

```bash
make data
```

I highly recommend adding your data processing functionality to the existing data commands
and adding any other functionality as subcommands of the main command.

## Project Organization

```text
├── LICENSE
├── Makefile           <- Makefile with commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
│
├── docs               <- A default Sphinx project; see sphinx-doc.org for details
│
├── models             <- Trained and serialized models, model predictions, or model summaries
│
├── notebooks          <- Jupyter or RMarkdown notebooks.
|                         Naming convention is a number (for ordering), the creator's initials,
|                         and a short `-` delimited description,
|                         e.g. `1.0-jqp-initial-data-exploration`.
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
├── bioinformatics_textbook   <- Source code for python package used in this project.
│   ├── __init__.py    <- Makes bioinformatics_textbook a Python module
│   ├── cli.py         <- Defines a CLI for the package. 
│   │                     Access from the command line with bioinformatics-textbook
│   │
│   ├── data           <- Scripts to download, generate, and/or process data.
│   │   └── __init__.py         <- Exposes the API for the data module.
│   │   └── download_data.py    <- Download data from external sources.
│   │   └── make_dataset.py     <- Make interim dataset from raw and/or external data.
│   │   └── process_dataset.py  <- Process interim dataset into a final processed dataset.
│   │
│   ├── features       <- Scripts to turn raw data into features for modeling
│   │   └── __init__.py         <- Exposes the API for the features module.
│   │   └── build_features.py
│   │
│   ├── models         <- Scripts to train models and then use trained models to make
│   │   └── __init__.py         <- Exposes the API for the models module.
│   │   │                 predictions
│   │   ├── predict_model.py
│   │   └── train_model.py
│   │
│   └── visualization  <- Scripts to create exploratory and results oriented visualizations
│   │   └── __init__.py         <- Exposes the API for the visualization module.
│       └── visualize.py
│
└── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io
```

Project based on the [cookiecutter data science project template][def].

[def]: https://github.com/trev-f/cookiecutter-data-science
