ARG PYTHON_VERSION="3.6.9"
FROM python:${PYTHON_VERSION}

# install vim text editor for git commit messages
RUN apt-get update \
    && apt-get install --no-install-recommends -y vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set working directory
WORKDIR /workspaces/bioinformatics_textbook

# install dependencies
COPY ./requirements.txt .
COPY ./setup.py .
ARG PYTHON_INTERPRETER="python3"
RUN ${PYTHON_INTERPRETER} -m pip install -U pip setuptools wheel
RUN ${PYTHON_INTERPRETER} -m pip install -r requirements.txt
