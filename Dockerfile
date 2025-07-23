
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
    time \
    imagemagick

# Install any needed packages specified in requirements.txt
# RUN pip install --trusted-host pypi.python.org -r requirements.txt
WORKDIR /fenicsR13
ADD ./requirements.txt /fenicsR13/requirements.txt
RUN pip install -r /fenicsR13/requirements.txt

# Install the fenicsR13 package (puts it into the PATH)
ADD . /fenicsR13
RUN pip install --editable /fenicsR13/.

# Replace default FEniCS Docker WELCOME screen with custom WELCOME screen
ADD WELCOME .
RUN echo "Built: $(date)" >> WELCOME
