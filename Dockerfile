
# Use FEniCS base image
FROM quay.io/fenicsproject/stable:2019.1.0.r3

# Descriptions
LABEL maintainer="Lambert Theisen <lambert.theisen@rwth-aachen.de>"
LABEL description="Linearized R13 Equations Solver Environment"

# Specify software versions
ENV GMSH_VERSION 4.4.0

# Download Install Gmsh SDK with dependecies from Github's dolfinx Dockerfile
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get -qq update && \
    # apt-get -yq --with-new-pkgs -o Dpkg::Options::="--force-confold" upgrade && \ # skip 110 package updates/upgrades
    apt-get -y install \
        libglu1 \
        libxcursor-dev \
        libxinerama1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN cd /usr/local && \
    wget -nc http://gmsh.info/bin/Linux/gmsh-${GMSH_VERSION}-Linux64-sdk.tgz && \
    tar -xf gmsh-${GMSH_VERSION}-Linux64-sdk.tgz
ENV PATH=/usr/local/gmsh-${GMSH_VERSION}-Linux64-sdk/bin:$PATH

# Install additional programs
RUN \
    apt-get update && \
    apt-get install -y \
    numdiff \
    htop \
    imagemagick

# Install any needed packages specified in requirements.txt
# RUN pip install --trusted-host pypi.python.org -r requirements.txt
# COPY requirements.txt /tmp/
RUN pip install \
    # Program related
    pyyaml>=5.1.1 \
    cerberus>=1.3.1 \
    pytest>=5.0.1 \
    pytest-cov>_2.7.1 \
    Sphinx>=2.1.2 \
    sphinxcontrib-napoleon>=0.7 \
    Pygments>=2.4.2 \
    # IDE related (can be skipped)
    pylint>=2.3.1 \
    rope>=0.14.0 \
    doc8>=0.8.0 \
    autopep8>=1.4.4 \
    pydocstyle>=4.0.1

# Install the fenicsR13 package (puts it into the PATH)
COPY . /tmp/
RUN pip install --editable /tmp/.

# Replace default FEniCS Docker WELCOME screen with custom WELCOME screen
COPY WELCOME .
RUN echo "Built: $(date)" >> WELCOME
