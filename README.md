# Qsar Tutorial

QSAR Tutorial - Obtaining data from ChEMBL and playing around with RDKit

There are two ways you can follow the code in this tutorial:

# Simple Set Up

1. Install Python with [Anaconda](https://www.anaconda.com/distribution/). Make sure you have Python 3.6+ installed.

2. Install RDKit with `conda`:

```console
  $ conda install -y -c rdkit rdkit
```

3. Install the tutorial dependencies with pip:

```console
  $ pip3 install notebooks/requirements.txt
```

4. Run a jupyter server and open the notebooks on the browser

# Docker Set Up

Alternatively, you can run the notebooks with Docker and Docker Compose.
This is more reliable since all dependencies are built within the Docker container.

1. Install [Docker](https://www.docker.com) and [Docker Compose](https://docs.docker.com/compose/overview/).

2. Build the container:

```console
  $ docker-compose build
```

3. Start a Jupyter notebook server:

```console
  $ docker-compose up
```

4. Follow the link to the jupyter server on your browser. It is usually a URL like http://localhost:8888?token=<huge_token>
