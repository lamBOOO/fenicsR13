
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
    apt-get -yq --with-new-pkgs -o Dpkg::Options::="--force-confold" upgrade && \
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

# OLD:

# Download gmsh
# ADD http://gmsh.info/bin/Linux/gmsh-$GMSH_VERSION-Linux64.tgz /home/fenics

# Download and extract gmsh
# RUN apt-get update && apt-get install -y cmake curl g++ gfortran libfltk1.3-dev libfreetype6-dev libgl1-mesa-dev liblapack-dev libxi-dev libxmu-dev mesa-common-dev tcl-dev tk-dev

# RUN apt-get update && apt-get install -y libgl1-mesa-dev
# ADD http://gmsh.info/bin/Linux/gmsh-$GMSH_VERSION-Linux64.tgz /home/fenics/
# RUN tar -xzf /home/fenics/gmsh-$GMSH_VERSION-Linux64.tgz \
#   && rm /home/fenics/gmsh-$GMSH_VERSION-Linux64.tgz

# USER root

# Install any needed packages specified in requirements.txt
# RUN pip install --trusted-host pypi.python.org -r requirements.txt
# COPY requirements.txt /tmp/
# RUN pip install --requirement /tmp/requirements.txt

# RUN apt-get update && apt-get install -y gmsh
