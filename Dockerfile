ARG PYTHON_VERSION="3.7"
FROM python:${PYTHON_VERSION}

# setup user info
ARG USERNAME=treevooor
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# set user to root for installs
USER root

# create the non-root user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # add sudo support
    && apt-get update \
    && apt-get install --no-install-recommends -y sudo \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# install vim text editor for git commit messages
RUN apt-get update \
    && apt-get install --no-install-recommends -y vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set default user to non-root
USER $USERNAME

# set working directory
WORKDIR /workspaces/bioinformatics_textbook

# install dependencies
COPY ./requirements.txt .
COPY ./setup.py .
ARG PYTHON_INTERPRETER="python3"
RUN ${PYTHON_INTERPRETER} -m pip install --user -U pip setuptools wheel
RUN ${PYTHON_INTERPRETER} -m pip install -r requirements.txt

# default to bash
ENTRYPOINT ["bash"]
