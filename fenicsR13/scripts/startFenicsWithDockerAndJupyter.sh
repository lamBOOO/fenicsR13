#!/bin/bash

docker run -w /home/fenics -v $(pwd):/home/fenics/shared -d -p 127.0.0.1:8888:8888 quay.io/fenicsproject/stable:current 'jupyter-notebook --ip=0.0.0.0'

# docker run --name notebook -w /home/fenics -v $(pwd):/home/fenics/shared -d -p 127.0.0.1:8888:8888 quay.io/fenicsproject/stable:current 'jupyter-notebook --ip=0.0.0.0'

# Start previously created container
# docker start notebook

# Connect to a Docker container
# $ docker exec -it notebook bash

# Attach to stdout of docker container
# $ docker attach notebook